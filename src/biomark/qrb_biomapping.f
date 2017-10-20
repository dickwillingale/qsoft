*+QRB_BIOMAPPING        Construct biomarker mapping
        SUBROUTINE QRB_BIOMAPPING(NP,NS,ARRAY,ICLASS,FBOT,
     +  FMIN,FDIF,FPLU,FMAX,PMIN,PDIFF,PMAX,PHI,SIGNA,FISH)
        INTEGER NP,NS,ICLASS(NS),SIGNA(NP,2)
        DOUBLE PRECISION ARRAY(NP,NS),FMIN(NP),FDIF(NP),FPLU(NP)
        DOUBLE PRECISION FMAX(NP)
        DOUBLE PRECISION PMIN(NP),PDIFF(NP),PMAX(NP),FBOT(NP),FISH(NP)
        DOUBLE PRECISION PHI(NP,NS)
Cf2py    intent(in) NP,NS,ARRAY,ICLASS,FBOT,FMIN,FDIF
Cf2py    intent(in) FPLU,FMAX,PMIN,PDIFF,PMAX
Cf2py    intent(out) PHI,SIGMA,FISH
*NP        input        number of peaks
*NS        input        number of individuals
*ARRAY     input        flux values of peaks
*ICLASS    input        class -1 1 or 0
*FBOT      input        lowest flux for peak
*FMIN      input        flux at lower disciminator
*FDIF      input        flux at maximum fractional difference
*FPLU      input        flux at upper disciminator
*FMAX      input        maximum flux for peak
*PMIN      input        fraction below lower disciminator
*PDIFF     input        maximum fractional difference
*PMAX      input        fraction above upper disciminator
*PHI       output        mapped marker flux values
*SIGNA     output        signatures of peaks
*FISH      output        Fisher score for peaks
* QR version RW 2014-Feb-16
* qsoft version RW 2017-Jun-26
*-Author Dick Willingale 2005-Mar-05
        INCLUDE 'QR_COM'
        INTEGER J,K,N1,N2
        EXTERNAL SX_CHISQD
        REAL SX_CHISQD
        DOUBLE PRECISION THR,CHIS,M1,M2,V1,V2
C Check status
        IF(ISTAT.NE.0) RETURN
C Set threshold using Chi-squared assuming number of individuals in 
C class is half the total and 99% confidence
        CHIS=SX_CHISQD(0.01,1)
        THR=2.0*CHIS/(CHIS+REAL(NS))
C        write(*,*) 'threshold ',THR,THR*NS*0.5
C Loop for all individuals
        DO J=1,NS
C Loop for all peaks
                DO K=1,NP
                        SIGNA(K,1)=0
                        IF(ABS(PMIN(K)).GT.THR) THEN
                                SIGNA(K,1)=SIGN(1.0D0,PMIN(K))
                        ENDIF
                        SIGNA(K,2)=0
                        IF(ABS(PMAX(K)).GT.THR) THEN
                                SIGNA(K,2)=SIGN(1.0D0,PMAX(K))
                        ENDIF
                        IF(ARRAY(K,J).LE.FMIN(K)) THEN
                                IF(FMIN(K).GT.FBOT(K)) THEN
                                     PHI(K,J)=(ARRAY(K,J)-FBOT(K))
     +                               /(FMIN(K)-FBOT(K))
                                     PHI(K,J)=0.5*PHI(K,J)-1.0
                                ELSE
                                        PHI(K,J)=-1.0
                                ENDIF
                        ELSEIF(ARRAY(K,J).LE.FDIF(K)) THEN
                                PHI(K,J)=(ARRAY(K,J)-FMIN(K))
     +                          /(FDIF(K)-FMIN(K))
                                PHI(K,J)=PHI(K,J)*0.5-0.5
                        ELSEIF(ARRAY(K,J).LE.FPLU(K)) THEN
                                PHI(K,J)=(ARRAY(K,J)-FDIF(K))
     +                          /(FPLU(K)-FDIF(K))
                                PHI(K,J)=PHI(K,J)*0.5
                        ELSE
                                PHI(K,J)=(ARRAY(K,J)-FPLU(K))
     +                          /(FMAX(K)-FPLU(K))
                                PHI(K,J)=PHI(K,J)*0.5+0.5
                        ENDIF
                        IF(PDIFF(K).EQ.0.0) THEN
                                PHI(K,J)=MAX(0.0,PHI(K,J))
                        ENDIF
                ENDDO
        ENDDO
C Calculate Fisher score for each peak
        DO K=1,NP
                M1=0.0
                M2=0.0
                V1=0.0
                V2=0.0
                N1=0
                N2=0
                DO J=1,NS
                        IF(ICLASS(J).EQ.1) THEN
                                M1=M1+PHI(K,J)
                                V1=V1+PHI(K,J)**2
                                N1=N1+1
                        ELSEIF(ICLASS(J).EQ.-1) THEN
                                M2=M2+PHI(K,J)
                                V2=V2+PHI(K,J)**2
                                N2=N2+1
                        ENDIF
                ENDDO
                M1=M1/N1
                M2=M2/N2
                V1=V1/N1-M1**2
                V2=V2/N2-M2**2
                FISH(K)=(M1-M2)**2/(V1+V2)
        ENDDO
        END
