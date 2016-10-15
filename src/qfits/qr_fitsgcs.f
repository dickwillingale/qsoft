*+QR_FITSGCS        Get FITS column value as string
        SUBROUTINE QR_FITSGCS(ICOL,IROW,IS,STR,SI,IUDEF)
        IMPLICIT NONE
        INTEGER ICOL,IROW,IS,SI,IUDEF
        CHARACTER STR*(IS)
*ICOL        input        column number
Cf2py  intent(in) icol
*IROW        input        row number
Cf2py  intent(in) irow
*IS        input        length of value string
Cf2py  intent(in) is
*STR        output        string returned
Cf2py  intent(out) str
*SI        output        number of characters set in string
Cf2py  intent(out) si
*IUDEF        output        1 if not defined
Cf2py  intent(out) iudef
*-Author: Dick Willingale 2012-Dec-27
        INCLUDE 'QR_COM'
        INTEGER LEN_TRIM
        LOGICAL LVAL
        EXTERNAL LEN_TRIM
C
        IF(ISTAT.NE.0) RETURN
C
        IF(IFITS.EQ.0) THEN
                WRITE(*,*) 'QR_FITSGCS error - file closed'
                ISTAT=1
                RETURN
        ENDIF
C
        CALL FTGCVS(IFITS,ICOL,IROW,1,1,'_',STR(1:IS),LVAL,ISTAT)
        SI=MIN(LEN_TRIM(STR),IS)
        IF(LVAL) THEN
                IUDEF=1
        ELSE
                IUDEF=0
        ENDIF
        END
