*+QRB_BIOMARKERS        Search for biomarkers
        SUBROUTINE QRB_BIOMARKERS(NP,NS,NW,ARRAY,ICLASS,
     +  WK1,WK2,FBOT,FMIN,FDIFF,FPLU,FMAX,PMIN,PDIFF,PMAX,
     +  CHIS,PKS,PCHIS)
        INTEGER NP,NS,NW,ICLASS(NS)
        DOUBLE PRECISION ARRAY(NP,NS),FMIN(NP),FDIFF(NP)
        DOUBLE PRECISION FPLU(NP),FMAX(NP)
        DOUBLE PRECISION PMIN(NP),PDIFF(NP),PMAX(NP),PKS(NP),PCHIS(NP)
        DOUBLE PRECISION WK1(NW),WK2(NW),FBOT(NP),CHIS(NP)
Cf2py    intent(in) NP,NS,NW,ARRAY,ICLASS,WK1,WK2
Cf2py    intent(out) CFBOT,FMIN,FDIFF,FLPU,FMAX
Cf2py    intent(out) PMIN,PDIFF,PMAX,CHIS,PKS,PCHIS
*NP        input        number of peaks
*NS        input        number of individuals
*NW        input        dimension of work arrays max(NP,NS)
*ARRAY     input        flux values of peaks
*ICLASS    input        class of individual 1, 0 (to ignore) or -1
*WK1       input        double precision work array
*WK2       input        double precision work array
*FBOT      output        lowest flux for each peak
*FMIN      output        flux at lower disciminator
*FDIFF     output        flux at maximum fractional difference
*FPLU      output        flux at upper disciminator
*FMAX      output        maximum flux for each peak
*PMIN      output        fraction below lower disciminator
*PDIFF     output        maximum fractional difference
*PMAX      output        fraction above upper disciminator
*CHIS      output        Chi-squared statistic
*PKS       output        Kolmogorov-Smirnov probability for each peak
*PCHIS     output        Chi-squared probability each peak
* QR version RW 2014-Feb-16
* qsoft version RW 2017-Jun-26
*-Author Dick Willingale 2005-Mar-05
        INCLUDE 'QR_COM'
        INTEGER J,K,I,DOF,IOVER,ITAIL
        DOUBLE PRECISION DF,NEFF,KSFUN,FM1,FM2,V1,V2
        DOUBLE PRECISION F(4),FD(4)
        REAL RPS,SR,RS,R(3),S(3),SX_PCHISQ,RCHIS
        EXTERNAL SX_PCHISQ
        EXTERNAL KSFUN
C Check status
        IF(ISTAT.NE.0) RETURN
C Find number of each class
        NM1=0
        NP1=0
        DO J=1,NS
                IF(ICLASS(J).EQ.1) THEN
                        NP1=NP1+1
                ELSEIF(ICLASS(J).EQ.-1) THEN
                        NM1=NM1+1
                ENDIF
        ENDDO
        NEFF=NP1*NM1
        NEFF=NEFF/(NP1+NM1)
C Set size of tail
        ITAIL=NEFF/100
C Loop for all peaks
        DO J=1,NP
C get fluxes and class into work arrays
                DO K=1,NS
                        WK1(K)=ARRAY(J,K)
                        WK2(K)=ICLASS(K)
                        IF(ICLASS(K).EQ.1) THEN
                                WK2(K)=WK2(K)/NP1
                        ELSEIF(ICLASS(K).EQ.-1) THEN
                                WK2(K)=WK2(K)/NM1
                        ENDIF
                ENDDO
C sort fluxes into ascending order
                CALL PDA_DSORT(WK1,WK2,NS,2,ISTAT)
C find disciminator flux levels
                NPD=0
                NND=0
                DF=0.0
                PDIFF(J)=0.0
                FDIFF(J)=0.0
                FM1=0.0
                FM2=0.0
                V1=0.0
                V2=0.0
                DO K=1,NS
                        IF(WK2(K).GT.0) THEN
                                FM1=FM1+WK1(K)
                                V1=V1+WK1(K)**2
                                NPD=NPD+1
                                IF(NPD.EQ.1+ITAIL) THEN
                                        F(1)=WK1(K)
                                ELSEIF(NPD.EQ.NP1-ITAIL) THEN
                                        F(2)=WK1(K)
                                ENDIF
                        ELSE
                                FM2=FM2+WK1(K)
                                V2=V2+WK1(K)**2
                                NND=NND+1
                                IF(NND.EQ.1+ITAIL) THEN
                                        F(3)=WK1(K)
                                ELSEIF(NND.EQ.NM1-ITAIL) THEN
                                        F(4)=WK1(K)
                                ENDIF
                        ENDIF
                        IF(WK2(K).NE.0) THEN
                                DF=DF+WK2(K)
                                IF(ABS(DF).GT.ABS(PDIFF(J))) THEN
                                        FDIFF(J)=WK1(K)
                                        PDIFF(J)=DF
                                ENDIF
                        ENDIF
                ENDDO
C sort levels into order
                FD(1)=1
                FD(2)=1
                FD(3)=-1
                FD(4)=-1
                CALL PDA_DSORT(F,FD,4,2,ISTAT)
                FBOT(J)=F(1)
                FMIN(J)=F(2)
                FPLU(J)=F(3)
                FMAX(J)=F(4)
C Check for overlap
                IOVER=1
                IF(FD(1)+FD(2).NE.0.0) THEN
                        IOVER=0
                ENDIF
C Check that FDIFF is between FMIN and FPLU
                IF(FDIFF(J).LT.FMIN(J).OR.FDIFF(J).GT.FPLU(J)) THEN
                        FDIFF(J)=(FMIN(J)+FPLU(J))*0.5
                ENDIF
C find fractions in wings
                PMIN(J)=0.0
                PMAX(J)=0.0
                DO K=1,NS
                        IF(WK2(K).NE.0) THEN
                            IF(IOVER.EQ.1) THEN
                                IF(WK1(K).LT.FMIN(J)) THEN
                                        PMIN(J)=PMIN(J)+ABS(WK2(K))
                                ELSEIF(WK1(K).GT.FPLU(J)) THEN
                                        PMAX(J)=PMAX(J)+ABS(WK2(K))
                                ENDIF
                            ELSE
                                IF(WK1(K).LE.FMIN(J)) THEN
                                        PMIN(J)=PMIN(J)+ABS(WK2(K))
                                ELSEIF(WK1(K).GE.FPLU(J)) THEN
                                        PMAX(J)=PMAX(J)+ABS(WK2(K))
                                ENDIF
                            ENDIF
                        ENDIF
                ENDDO
                PMIN(J)=PMIN(J)*SIGN(1.0D0,FD(1))
                PMAX(J)=PMAX(J)*SIGN(1.0D0,FD(4))
C Calculate mean fluxes
                FM1=FM1/NP1
                FM2=FM2/NM1
C Calculate variances
                V1=V1/NP1-FM1**2
                V2=V2/NM1-FM2**2
C Calculate Kolomogorov-Smirnov probability using maximum fractional difference
                PKS(J)=KSFUN(NEFF,ABS(PDIFF(J)))
C Calculate Chi-squared using counts outside discriminators
                DOF=2
                IF(PMIN(J).EQ.0) THEN
                        DOF=DOF-1
                ENDIF
                IF(PMAX(J).EQ.0) THEN
                        DOF=DOF-1
                ENDIF
                RS=SQRT(REAL(NP1)/REAL(NM1))
                SR=1.0/RS
                R(1)=0
                S(1)=0
                R(2)=0
                S(2)=0
                R(3)=0
                S(3)=0
                IF(PMIN(J).GT.0) THEN
                        R(1)=PMIN(J)*NP1
                        S(1)=ITAIL
                        R(2)=NP1-R(1)
                        S(2)=NM1-S(1)
                ELSEIF(PMIN(J).LT.0) THEN
                        R(1)=ITAIL
                        S(1)=-PMIN(J)*NM1
                        R(2)=NP1-R(1)
                        S(2)=NM1-S(1)
                ENDIF
                IF(PMAX(J).GT.0) THEN
                        R(3)=PMAX(J)*NP1
                        S(3)=ITAIL
                        R(2)=NP1-R(3)
                        S(2)=NM1-S(3)
                ELSEIF(PMAX(J).LT.0) THEN
                        R(3)=ITAIL
                        S(3)=-PMAX(J)*NM1
                        R(2)=NP1-R(3)
                        S(2)=NM1-S(3)
                ENDIF
                IF(DOF.GT.0) THEN
                        CHIS(J)=0.0
                        DO I=1,3
                            RPS=R(I)+S(I)
                            IF(RPS.GT.0) THEN
                               CHIS(J)=CHIS(J)+(SR*R(I)-RS*S(I))**2/RPS
                            ENDIF
                        ENDDO
                        RCHIS=CHIS(J)
                        PCHIS(J)=SX_PCHISQ(RCHIS,DOF)
                ELSE
                        CHIS(J)=0.0
                        PCHIS(J)=1.0
                ENDIF
C flip sign of PDIFF
                PDIFF(J)=-PDIFF(J)
        ENDDO
        END
*+KSFUN        Kolmogorov-Smirnov probablity function
        FUNCTION KSFUN(ANE,D)
        IMPLICIT NONE
        DOUBLE PRECISION KSFUN,ANE,D
*ANE    effective number of data points n1*n2/(n1+n2)
*D       maximum fractional difference in distributions
*-Author Dick Willingale from Numerical Recipes 1997-May-5
        DOUBLE PRECISION ALAM,EPS1,EPS2,A2,FAC,TERM,TERMBF,AT,SR
        PARAMETER (EPS1=0.001,EPS2=1.D-8)
        INTEGER J
C
        IF(ANE.GT.0.0) THEN
                SR=SQRT(ANE)
        ELSE
                KSFUN=1.0
                RETURN
        ENDIF
        ALAM=D*(SR+0.12+0.11/SR)
        A2=-2.0*ALAM**2
        FAC=2.0
        KSFUN=0.0
        TERMBF=0.0
        DO J=1,100
                TERM=FAC*EXP(A2*J**2)
                KSFUN=KSFUN+TERM
                AT=ABS(TERM)
                IF(AT.LE.EPS1*TERMBF.OR.AT.LE.EPS2*KSFUN) RETURN
                FAC=-FAC
                TERMBF=AT
        ENDDO
C Failed to converge
        KSFUN=1.0
        END
