*+QRT_C1NEST Set conical appoximation to Wolter type 1 nest
        SUBROUTINE QRT_C1NEST(AX,AR,FF,IQ,IB,IDF,
     +  NS,PL,PH,HL,HH,RPL,RPH,RHL,RHH,TIN,TJ,TOUT,NW,WRK)
        IMPLICIT NONE
        DOUBLE PRECISION AX(3),AR(3),FF(3)
        INTEGER IQ,IB,IDF
        INTEGER NS,NW
        DOUBLE PRECISION PL,PH,HL,HH,RPL(NS),RPH(NS),RHL(NS),RHH(NS)
        DOUBLE PRECISION TIN(NS),TJ(NS),TOUT(NS),WRK(NW)
*AX        input        optical axis
*AR        input        reference axis
*FF        input        position of focus
*NS        input        number of shells
*PL        input        low axial position of parabola
*PH        input        high axial position of parabola
*HL        input        low axial position of hyperbola
*HH        input        high axial position of hyperbola
*RPL       input        radii parabola near join
*RPH       input         radii parabola at input aperture
*RHL       input        radii hyperbola at exit aperture
*RJH       input         radii hyperbola near join
*TIN       input        thicknesses of shells at input aperture
*TJ        input        thicknesses of shells at join plane (or near join)
*TOUT      input        thicknesses of shells at the output aperture
*NW        input        size of workspace
*WRK       in/out        workspace array
Cf2py  intent(in) AX,AR,FF,IQ,IB,IDF,NS,PL,PH,HL,HH,RPL,RPH,RHL,RHH
Cf2py  intent(in) TIN,TJ,TOUT,NW,WRK
*-Author Dick Willingale 2012-Jun-28
        INTEGER K,J,IDEF(2),N,KSUR,KMIS,JUMP
        DOUBLE PRECISION PP(16),PY(16)
        DOUBLE PRECISION CP(16),CH(16),EA(12)
        INCLUDE 'QR_COM'
C
        IF(ISTAT.NE.0) RETURN
C Check workspace
        IF(NW.LT.(NS-1)*2+12) THEN
           WRITE(*,*) 'qrt_c1nest error - insufficient workspace'
           ISTAT=1
           RETURN
        ENDIF
C
        IDEF(1)=IDF
        IDEF(2)=NS-1
C Set front aperture parameters except for radial limits
C Set exit aperture parameters except for radial limits
        DO J=1,3
                WRK(J)=AX(J)
                EA(J)=AX(J)
        ENDDO
        DO J=1,3
                WRK(J+3)=AR(J)
                EA(J+3)=AR(J)
        ENDDO
        DO J=1,3
                WRK(J+6)=FF(J)+AX(J)*PH
                EA(J+6)=FF(J)+AX(J)*HL
        ENDDO
        WRK(10)=PH
        EA(10)=HL
C Find next surface and surface after nest
        CALL SRT_NSUR(KSUR,ISTAT)
        JUMP=5
        KMIS=KSUR+(NS-1)*JUMP+1
        N=(NS-1)*2+10
        DO J=11,N
                WRK(J)=0.0
        ENDDO
C Set this surface index and jump (number of surfaces per aperture of nest)
        WRK(N+1)=KSUR
        WRK(N+2)=JUMP
        N=N+2
C Set input aperture parameters in common
        CALL SRT_SETF(0,4,N,WRK,IDEF,0,0,-1,ISTAT)
C Loop for shells (last shell is a dummy for stops only)
        DO K=1,NS
C Generate conic coefficients and limits for shell
                CALL SRT_CONE(PL,RPL(K),PH,RPH(K),AX,AR,FF,PP)
                CALL SRT_CONE(HL,RHL(K),HH,RHH(K),AX,AR,FF,PY)
                IF(K.LT.NS) THEN
C Generate coefficients for back of next shell
                        CALL SRT_CONE(PL,RPL(K+1)+TJ(K+1)
     +                        ,PH,RPH(K+1)+TIN(K+1),AX,AR,FF,CP)
                        CALL SRT_CONE(HL,RHL(K+1)+TOUT(K+1)
     +                        ,HH,RHH(K+1)+TJ(K+1),AX,AR,FF,CH)
C Set parameters for parabola
                        IDEF(2)=K
                        CALL SRT_SETF(0,9,16,PP,IDEF,IQ,-1,-1,ISTAT)
C Set parameters in common for back of next parabola
                        IDEF(2)=K+1
                        CALL SRT_SETF(0,9,16,CP,IDEF,IB,
     +                  KSUR+(K-1)*JUMP+1,-1,ISTAT)
C Set parameters in common for hyperbola
                        IDEF(2)=K
                        CALL SRT_SETF(0,9,16,PY,IDEF,IQ,-1,-1,ISTAT)
C Set parameters in common for back of next hyperbola
                        IDEF(2)=K+1
                        CALL SRT_SETF(0,9,16,CH,IDEF,IB,
     +                  KSUR+(K-1)*JUMP+3,-1,ISTAT)
C Set exit aperture
                        EA(11)=CH(14)
                        EA(12)=PY(14)
                        IDEF(2)=K
                        CALL SRT_SETF(0,3,12,EA,IDEF,0,0,KMIS,ISTAT)
C Set maximum radius of input aperture for current shell
                        WRK(12+(K-1)*2)=RPH(K)
                ENDIF
C Set mimumum radius of input aperture for previous shell
                IF(K.GT.1) THEN
                        WRK(11+(K-2)*2)=RPH(K)+TIN(K)
                ENDIF
        ENDDO
C Update input aperture parameters
        IDEF(2)=NS-1
        CALL SRT_SETF(KSUR,4,N,WRK,IDEF,0,0,-1,ISTAT)
C Check next free surface is OK
        CALL SRT_NSUR(KSUR,ISTAT)
        IF(KMIS.NE.KSUR) THEN
                WRITE(*,*) 'qrt_c1nest error - next surface ',KSUR,
     +         ' not ',KMIS
                ISTAT=1
                RETURN
        ENDIF
        END
