*+QRB_BIOSCORE        perform scoring of classification using biomarkers
        SUBROUTINE QRB_BIOSCORE(NS,RATE,ICLASS,WK1,WK2,
     +  CMAT,UNAMB,PSEP,ROC)
        INTEGER NS,ICLASS(NS),CMAT(2,2),UNAMB
        DOUBLE PRECISION RATE(NS),PSEP,WK1(NS),WK2(NS),ROC
Cf2py    intent(in) NS,RATE,ICLASS
Cf2py    intent(inout) WK1,WK2
Cf2py    intent(out) CMAT,UNAMB,PSEP,ROC
*NS        input        number of individuals
*RATE      input        rate of occurance of markers for each individual
*ICLASS    input        class of individual 1, 0 (to ignore) or -1
*WK1       in/out        working array
*WK2       in/out        working array
*CMAT      output        confusion matrix
*UNAMB     output        unambiguous count
*PSEP      output        peak separation
*ROC       output        Relative Operating Characteristic score
* QR version RW 2014-Feb-16
* qsoft version RW 2017-Jun-26
*-Author Dick Willingale 2006-Feb-07
        INCLUDE 'QR_COM'
        DOUBLE PRECISION FMIN,FMAX
        INTEGER YPLUS,YMINUS,NPLUS,NMINUS
C Check status
        IF(ISTAT.NE.0) RETURN
C Analyse distribution of weights
        YPLUS=0
        YMINUS=0
        NPLUS=0
        NMINUS=0
        FMIN=0.0
        FMAX=0.0
C Loop for instances
        DO I=1,NS
                IF(ICLASS(I).GT.0) THEN
                        IF(RATE(I).LE.0) THEN
                                NPLUS=NPLUS+1
                                FMIN=MIN(FMIN,RATE(I))
                        ELSE
                                YPLUS=YPLUS+1
                        ENDIF
                ELSEIF(ICLASS(I).LT.0) THEN
                        IF(RATE(I).GE.0) THEN
                                YMINUS=YMINUS+1
                                FMAX=MAX(FMAX,RATE(I))
                        ELSE
                                NMINUS=NMINUS+1
                        ENDIF
                ENDIF
        ENDDO
C Calculate peak separation
        PSEP=FMIN-FMAX
C Sum rate outside of range FMIN to FMAX
        UNAMB=0
        DO I=1,NS
                IF(ICLASS(I).GT.0) THEN
                        IF(RATE(I).GT.FMAX) THEN
                                UNAMB=UNAMB+1
                        ENDIF
                ELSEIF(ICLASS(I).LT.0) THEN
                        IF(RATE(I).LT.FMIN) THEN
                                UNAMB=UNAMB+1
                        ENDIF
                ENDIF
        ENDDO        
C Set components of confusion matrix
        CMAT(1,1)=YPLUS
        CMAT(1,2)=NPLUS
        CMAT(2,1)=YMINUS
        CMAT(2,2)=NMINUS
C Calculate Relative Operating Characteristic score
        CALL QRB_BIOROC(NS,RATE,ICLASS,WK1,WK2,ROC)
        END
