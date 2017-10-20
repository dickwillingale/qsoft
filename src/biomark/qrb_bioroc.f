*+QRB_BIOROC        Calculation of Relative Operating Characteristic score
        SUBROUTINE QRB_BIOROC(NS,RATE,ICLASS,WK1,WK2,ROC)
        INTEGER NS,ICLASS(NS)
        DOUBLE PRECISION RATE(NS),WK1(NS),WK2(NS),ROC
Cf2py    intent(in) NS,RATE,ICLASS
Cf2py    intent(inout) WK1,WK2
Cf2py    intent(out) ROC
*NS        input        number of individuals
*RATE      input        rate of occurance of markers for each individual
*ICLASS    input        class of individual 1, 0 (to ignore) or -1
*WK1       in/out        working array
*WK2       in/out        working array
*ROC       output        Relative Operating Characteristic score
* QR version RW 2014-Feb-16
* qsoft version RW 2017-Jun-26
*-Author Dick Willingale 2006-Feb-06
        INCLUDE 'QR_COM'
        DOUBLE PRECISION PLUS
        INTEGER I,NPLUS,NMINUS
C Find number in each class
        NPLUS=0
        NMINUS=0
        DO I=1,NS
                IF(ICLASS(I).EQ.1) THEN
                        NPLUS=NPLUS+1
                ELSEIF(ICLASS(I).EQ.-1) THEN
                        NMINUS=NMINUS+1
                ENDIF
        ENDDO
C Sort rates into descending order (classification to follow)
        DO I=1,NS
                WK1(I)=RATE(I)
                WK2(I)=ICLASS(I)
                IF(ICLASS(I).EQ.1) THEN
                        WK2(I)=WK2(I)/NPLUS
                ELSEIF(ICLASS(I).EQ.-1) THEN
                        WK2(I)=WK2(I)/NMINUS
                ENDIF
        ENDDO
        CALL PDA_DSORT(WK1,WK2,NS,-2,ISTAT)
C calculate Relative Operating Characteristic score
        PLUS=0.0
        ROC=0.0
        DO I=1,NS
                IF(WK2(I).GT.0.0) THEN
                        PLUS=PLUS+WK2(I)
                ELSEIF(WK2(I).LT.0.0) THEN
                        ROC=ROC-PLUS*WK2(I)
                ENDIF
        ENDDO
        END
