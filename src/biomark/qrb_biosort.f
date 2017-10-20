*+QRB_BIOSORT        Sort potential biomarkers
        SUBROUTINE QRB_BIOSORT(NP,NW,PDIFF,CHIS,
     +  FISH,WK1,WK2,IRANK,IMARK)
        IMPLICIT NONE
        INTEGER NP,NW,IMARK(NP),IRANK(NP)
        DOUBLE PRECISION PDIFF(NP),WK1(NW),WK2(NW)
        DOUBLE PRECISION CHIS(NP)
        DOUBLE PRECISION FISH(NP)
Cf2py    intent(in) NP,PDIFF,CHIS,SIGMA,FISH
Cf2py    intent(inout) WK1,WK2
Cf2py    intent(out) IRANK,IMARK
*NP        input        number of peaks
*PDIFF     input        maximum fractional difference
*CHIS      input        Chi-squared statistic
*FISH      input        array of Fisher scores for peaks
*WK1       input        work array
*WK2       input        work array
*IRANK     output        ranking index of potential markers (0 if not included)
*IMARK     output        >0 if peak potential marker, 0 if not
* QR version RW 2014-Feb-16
* qsoft version RW 2017-Jun-26
*-Author Dick Willingale 2005-Mar-05
        INCLUDE 'QR_COM'
        INTEGER J
        INTEGER N50
        DOUBLE PRECISION FS50,KS50,CHIS50
C Check status
        IF(ISTAT.NE.0) RETURN
C set 50% index
        N50=INT(NP*50/100)
C write(*,*) 'number of peaks',NP,'half',N50
C sort on Fisher score and pick 50% value
        DO J=1,NP
                WK1(J)=FISH(J)
                WK2(J)=J
        ENDDO
        CALL PDA_DSORT(WK1,WK2,NP,-2,ISTAT)
        FS50=WK1(N50)
C        write(*,*) 'Fisher Score median value',FS50
C sort on K-S statistic and pick 50% value
        DO J=1,NP
                WK1(J)=ABS(PDIFF(J))
                WK2(J)=J
        ENDDO
        CALL PDA_DSORT(WK1,WK2,NP,-2,ISTAT)
        KS50=WK1(N50)
C        write(*,*) 'K-S Statistic median value',KS50
C sort on Chi-Squared and pick 50% value
        DO J=1,NP
                WK1(J)=CHIS(J)
                WK2(J)=J
        ENDDO
        CALL PDA_DSORT(WK1,WK2,NP,-2,ISTAT)
        CHIS50=WK1(N50)
C        write(*,*) 'Chi-Squared Statistic median value',CHIS50
C Calculate figure of merit and sort
        DO J=1,NP
C Combinations of K-S and Fisher
C                WK1(J)=SQRT((FISH(J)/FS50)**2+(PDIFF(J)/KS50)**2)/SQRT(2.0)
C Just K-S statistic
                WK1(J)=ABS(PDIFF(J))/KS50
C Just Fisher score
C                WK1(J)=ABS(FISH(J)/FS50)
C Just Chi-Squared
C                WK1(J)=ABS(CHIS(J)/CHIS50)
                WK2(J)=J
                IF(WK1(J).GT.1.0) THEN
                        IMARK(J)=1
                ELSE
                        IMARK(J)=0
                ENDIF
        ENDDO
        CALL PDA_DSORT(WK1,WK2,NP,-2,ISTAT)
C Pick out index as marker rank
        DO J=1,NP
                IRANK(J)=WK2(J)
        ENDDO
        END
