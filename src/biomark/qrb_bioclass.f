*+QRB_BIOCLASS        Classification using biomarkers
        SUBROUTINE QRB_BIOCLASS(NP,NS,NW,ARRAY,ICLASS,NI,IMARK,PDIFF,
     +  WK1,WK2,
     +  RATE,RMSDEV,CMAT,UNAMB,PSEP,ROC)
        INTEGER NP,NS,NW,ICLASS(NS),NI,IMARK(NI),CMAT(2,2),UNAMB
        DOUBLE PRECISION ARRAY(NP,NS),RATE(NS),RMSDEV(NS),PSEP
        DOUBLE PRECISION WK1(NW),WK2(NW)
        DOUBLE PRECISION PDIFF(NP),ROC
Cf2py    intent(in) NP,NS,NW,ARRAY,ICLASS,NI,IMARK,PDIFF
Cf2py    intent(inout) WK1,WK2
Cf2py    intent(out) RATE,RMSDEV,CMAT,UNAMB,PSEP,ROC
*NP        input        number of peaks
*NS        input        number of individuals
*NW        input        size of work array
*ARRAY     input        mapped intensity values of peaks range 0 to 1
*ICLASS    input        class of individual 1, 0 (to ignore) or -1
*NI        input        number of marker peaks in index
*IMARK     input        index of marker peaks
*PDIFF     input        maximum difference between cummulative distributions
*WK1       in/out        working array
*WK2       in/out        working array
*RATE      output        rate of occurance of markers for each individual
*RMSDEV    output        rms deviation of the rate
*CMAT      output        confusion matrix
*UNAMB     output        unambiguous count
*PSEP      output        peak separation
*ROC       output        Relative Operating Characteristic score
* QR version RW 2014-Feb-16
* qsoft version RW 2017-Jun-26
*-Author Dick Willingale 2005-Mar-05
        INCLUDE 'QR_COM'
        INTEGER J,NMARK
        DOUBLE PRECISION WTS
C Check status
        IF(ISTAT.NE.0) RETURN
C initialise rate arrays
        DO J=1,NS
                RATE(J)=0.0
                RMSDEV(J)=0.0
        ENDDO
C Loop for index of marker peaks
        NMARK=0
        DO J=1,NI
                K=IMARK(J)
                NMARK=NMARK+1
C Loop for instances
                DO I=1,NS
                        WTS=ARRAY(K,I)*SIGN(1.0D0,PDIFF(K))
                        RATE(I)=RATE(I)+WTS
                ENDDO
        ENDDO
C Calculate rms deviation of rates
        DO J=1,NS
                RMSDEV(J)=RMSDEV(J)/NMARK-(RATE(J)/NMARK)**2
                IF(RMSDEV(J).GT.0.0) THEN
                        RMSDEV(J)=SQRT(RMSDEV(J))
                ENDIF
        ENDDO
C Analyse distribution of weights
        CALL QRB_BIOSCORE(NS,RATE,ICLASS,WK1,WK2,
     +        CMAT,UNAMB,PSEP,ROC)
        END
