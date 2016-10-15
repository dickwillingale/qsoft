*+QRT_W1NEST Set Wolter type 1 nest parameters
        SUBROUTINE QRT_W1NEST(NS,XJ,RJ,RA,PL,PH,HL,HH,TIN,TJ,TOUT,
     +  AX,AR,FF,DEFI,IQ,IB)
*NS        input        number of shells
*XJ        input        axial position of join plane
*RJ        input        radii of shells at join
*RA        input        ratio of grazing angles
*PL        input        low axial position of parabola
*PH        input        high axial position of parabola
*HL        input        low axial position of hyperbola
*HH        input        high axial position of hyperbola
*TIN       input        thicknesses of shells at input aperture
*TJ        input        thicknesses of shells at join plane
*TOUT      input        thicknesses of shells at the output aperture
*AX        input        direction of axis
*AR        input        reference axis in aperture
*FF        input        position of focus
*DEFI      input        deformation index
*IQ        input        reflecting surface quality index
*IB        input        back of shells surface quality index
Cf2py  intent(in) NS,XJ,RJ,RA,PL,PH,HL,HH,TIN,TJ,TOUT
Cf2py  intent(in) AX,AR,FF,DEFI,IQ,IB
        IMPLICIT NONE
        INTEGER NS,DEFI,IQ,IB
        DOUBLE PRECISION XJ,RJ(NS),RA,PL,PH,HL,HH,TIN(NS),TJ(NS)
        DOUBLE PRECISION TOUT(NS),AX(3),AR(3),FF(3)
*-Author Dick Willingale 2012-May-1
        INTEGER MAXN,NWRK
        PARAMETER (MAXN=500,NWRK=(MAXN-1)*2+12)
        DOUBLE PRECISION RIN(MAXN),ROUT(MAXN),WRK(NWRK)
        INTEGER K,J,IDEF(2),N,KSUR,KMIS,JUMP
        DOUBLE PRECISION PP(16),PY(16),PPB(16),PYB(16)
        DOUBLE PRECISION CP(16),CH(16),EA(12)
        DOUBLE PRECISION XMIN,RMIN,XMAX,RMAX
        INCLUDE 'QR_COM'
C
        IF(ISTAT.NE.0) RETURN
C Check workspace
        IF(NS.GT.MAXN) THEN
                WRITE(*,*) 'qrt_w1nest error - insufficient workspace'
                ISTAT=1
                RETURN
        ENDIF
        IDEF(1)=DEFI
        IDEF(2)=NS-1
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
                CALL SRT_WLT1(XJ,RJ(K),RA,PL,PH,HL,HH,AX,AR,FF,PP,PY)
                RIN(K)=PP(16)
                ROUT(K)=PY(14)
                IF(K.LT.NS) THEN
C Generate coefficients for next shell
                        CALL SRT_WLT1(XJ,RJ(K+1),RA,PL,PH,HL,HH,
     +                  AX,AR,FF,PPB,PYB)
C Generate cone coefficients and limits for back of next shell
                        RMAX=PPB(16)+TIN(K+1)
                        RMIN=PPB(14)+TJ(K+1)
                        XMAX=PPB(15)
                        XMIN=PPB(13)
                        CALL SRT_CONE(XMIN,RMIN,XMAX,RMAX,AX,AR,FF,CP)
                        RMAX=PYB(16)+TJ(K+1)
                        RMIN=PYB(14)+TOUT(K+1)
                        XMAX=PYB(15)
                        XMIN=PYB(13)
                        CALL SRT_CONE(XMIN,RMIN,XMAX,RMAX,AX,AR,FF,CH)
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
                        EA(11)=PYB(14)+TOUT(K+1)
                        EA(12)=PY(14)
                        IDEF(2)=K
                        CALL SRT_SETF(0,3,12,EA,IDEF,0,0,KMIS,ISTAT)
C Set maximum radius of input aperture for current shell
                        WRK(12+(K-1)*2)=RIN(K)
                ENDIF
C Set mimumum radius of input aperture for previous shell
                IF(K.GT.1) THEN
                        WRK(11+(K-2)*2)=RIN(K)+TIN(K)
                ENDIF
        ENDDO
C Update input aperture parameters
        IDEF(2)=NS-1
        CALL SRT_SETF(KSUR,4,N,WRK,IDEF,0,0,-1,ISTAT)
C Check next free surface is OK
        CALL SRT_NSUR(KSUR,ISTAT)
        IF(KMIS.NE.KSUR) THEN
                WRITE(*,*) 'qrt_w1nest error - next surface ',
     +          KSUR,' not ',KMIS
                ISTAT=1
                RETURN
        ENDIF
        END
