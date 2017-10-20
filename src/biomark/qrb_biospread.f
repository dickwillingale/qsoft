*+QRB_BIOSPREAD        perform Gaussian spreading of 1-D array
        SUBROUTINE QRB_BIOSPREAD(N,DA,WK,NS,GSIG,RES)
        IMPLICIT NONE
        INTEGER N,NS
        DOUBLE PRECISION DA(N),WK(N),RES(N),GSIG(NS)
Cf2py    intent(in) N,DA,NS,GSIG
Cf2py    intent(inout) WK
Cf2py    intent(out) RES
*N    input     number of data points
*DA   input     data array to be spread
*WK   input     work array
*NS   input     number of sigma values 1 or N
*GSIG input     array of sigma values (in units of samples)
*RES  output    spread array
* qsoft version RW 2017-Jun-26
*-Author Dick Willingale 2014-Mar-30
        INTEGER K,J,KK,JJ,NH
        DOUBLE PRECISION SM        
C Initialize output array
        DO K=1,N
                RES(K)=0.0
        ENDDO
C Loop to spread
        NH=0
        DO K=1,N
                IF(K.EQ.1.AND.NS.EQ.1) THEN
C Set up Gaussian profile once if NS=1
                        SM=0.0
                        NH=0
                        J=1
                        DO WHILE(NH.EQ.0)
                                WK(J)=EXP(-REAL(J-1)**2/
     +                          (2.0*GSIG(K)**2))
                                IF(J.EQ.N.OR.WK(J).LT.0.0001) THEN
                                        NH=J
                                ENDIF
                                IF(J.EQ.1) THEN
                                        SM=SM+WK(J)
                                ELSE
                                        SM=SM+WK(J)*2.0
                                ENDIF
                                J=J+1
                        ENDDO
C Normalise so that sum of Gaussian samples is 1.0
                        DO J=1,NH
                                WK(J)=WK(J)/SM
                        ENDDO
                ELSEIF(NS.GT.1) THEN
C Set up Gaussian profile for each position if NS>1 (NS=N)
                        SM=0.0
                        NH=0
                        J=1
                        DO WHILE(NH.EQ.0)
                                WK(J)=EXP(-REAL(J-1)**2/
     +                          (2.0*GSIG(K)**2))
                                IF(J.EQ.N.OR.WK(J).LT.0.0001) THEN
                                        NH=J
                                ENDIF
                                IF(J.EQ.1) THEN
                                        SM=SM+WK(J)
                                ELSE
                                        SM=SM+WK(J)*2.0
                                ENDIF
                                J=J+1
                        ENDDO
C Normalise so that sum of Gaussian samples is 1.0
                        DO J=1,NH
                                WK(J)=WK(J)/SM
                        ENDDO
                ENDIF
                DO KK=K-NH+1,K+NH-1
                        IF(KK.GT.0.AND.KK.LE.N) THEN
                                JJ=ABS(KK-K)+1
                                RES(KK)=RES(KK)+DA(K)*WK(JJ)
                        ENDIF
                ENDDO
        ENDDO
        END
