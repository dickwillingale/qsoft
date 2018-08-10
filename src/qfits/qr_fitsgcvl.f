*+QR_FITSGCVL        Get FITS column logical
        SUBROUTINE QR_FITSGCVL(ICOL,IROW,NEL,VAL)
        IMPLICIT NONE
        INTEGER ICOL,IROW,NEL
        LOGICAL VAL(NEL)
*ICOL                input        column number
Cf2py  intent(in) icol
*IROW                input        1st row number
Cf2py  intent(in) irow
*NEL                input        number of elements
Cf2py  intent(in) nel
*VAL                output        values returned
Cf2py  intent(out) val
*-Author: Dick Willingale 2012-Dec-27
        INCLUDE 'QR_COM'
C
        IF(ISTAT.NE.0) RETURN
C
        IF(IFITS.EQ.0) THEN
                WRITE(*,*) 'QR_FITSGCVL error - file closed'
                ISTAT=1
                RETURN
        ENDIF
C
        CALL FTGCL(IFITS,ICOL,IROW,1,NEL,VAL,ISTAT)
        END
