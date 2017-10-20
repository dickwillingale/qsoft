*+QRB_BIOHUNT        Search for optimum set of biomarkers and classify
        SUBROUTINE QRB_BIOHUNT(NP,NS,NW,ARRAY,ICLASS,NI,IMARK,PDIFF,
     +  WK1,WK2,IWK,
     +  NMARK,IOUT,RATE,RMSDEV,CMAT,UNAMB,PSEP,ROC)
        INTEGER NP,NS,NW,ICLASS(NS),NI,IMARK(NI),CMAT(2,2)
        INTEGER UNAMB,IOUT(NI)
        INTEGER NMARK,IWK(NW)
        DOUBLE PRECISION ARRAY(NP,NS),RATE(NS),RMSDEV(NS),PSEP
        DOUBLE PRECISION WK1(NW),WK2(NW)
        DOUBLE PRECISION PDIFF(NP),ROC(NS)
Cf2py    intent(in) NP,NS,NW,ARRAY,ICLASS,NI,IMARK,PDIFF
Cf2py    intent(inout) WK1,WK2,IWK
Cf2py    intent(out) NMARK,IOUT,RATE,RMSDEV,CMAT,UNAMB,PSEP,ROC
*NP        input        number of peaks
*NS        input        number of individuals
*NW        input        size of work array
*ARRAY     input        mapped intensity values of peaks range -1 to 1
*ICLASS    input        class of individual 1, 0 (to ignore) or -1
*NI        input        number of marker peaks in index
*IMARK     input        input index of marker peaks
*PDIFF     input        maximum difference between cummulative distributions
*WK1       in/out        work array
*WK2       in/out        work array
*IWK       in/out        integer work array
*NMARK     output        optimum number of markers
*IOUT      output        output index of markers
*RATE      output        rate of occurance of markers for each individual
*RMSDEV    output        rms deviation of the rate
*CMAT      output        confusion matrix
*UNAMB     output        unambiguous count
*PSEP      output        peak separation
*ROC       output        Relative Operating Characteristic scores
* QR version RW 2014-Feb-16
* qsoft version RW 2017-Jun-26
*-Author Dick Willingale 2005-Sep-05
        INCLUDE 'QR_COM'
        INTEGER I,IBEST
        DOUBLE PRECISION SCOREMAX
C Check status
        IF(ISTAT.NE.0) RETURN
C transfer input marker index to work array
        DO I=1,NI
                IWK(I)=IMARK(I)
        ENDDO
C set initial marker from top of list
        NMARK=1
        IOUT(NMARK)=IWK(1)
        IWK(1)=0
C calculate initial score
        CALL QRB_BIOCLASS(NP,NS,NW,ARRAY,ICLASS,NMARK,IOUT,PDIFF,
     +        WK1,WK2,
     +        RATE,RMSDEV,CMAT,UNAMB,PSEP,ROC(NMARK))
        SCOREMAX=ROC(NMARK)
        IBEST=1
        DO WHILE(IBEST.NE.0)
                IBEST=0
                NMARK=NMARK+1
                DO I=1,NI
                        IF(IWK(I).NE.0) THEN
                                IOUT(NMARK)=IWK(I)
                                CALL QRB_BIOCLASS(NP,NS,NW,ARRAY,
     +                                ICLASS,NMARK,
     +                                IOUT,PDIFF,WK1,WK2,
     +                                RATE,RMSDEV,CMAT,UNAMB,
     +                                PSEP,ROC(NMARK))
                                IF(ROC(NMARK).GT.SCOREMAX) THEN
                                        IBEST=I
                                        SCOREMAX=ROC(NMARK)
                                ENDIF
                        ENDIF
                ENDDO
                IF(IBEST.GT.0) THEN
                        IOUT(NMARK)=IWK(IBEST)
                        IWK(IBEST)=0
                ELSE
                        NMARK=NMARK-1
                ENDIF
        ENDDO
C found optimum set of markers - apply them to return optimum classification
        CALL QRB_BIOCLASS(NP,NS,NW,ARRAY,ICLASS,NMARK,IOUT,PDIFF,
     +        WK1,WK2,
     +        RATE,RMSDEV,CMAT,UNAMB,PSEP,ROC(NMARK))
        END
