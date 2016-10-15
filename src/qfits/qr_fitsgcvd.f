*+QR_FITSGCVD        Get FITS column real (double)
        SUBROUTINE QR_FITSGCVD(NULLVAL,ICOL,IROW,NEL,VAL)
        IMPLICIT NONE
        INTEGER ICOL,IROW,NEL
        DOUBLE PRECISION NULLVAL,VAL(NEL)
*NULLVAL        input        undefined value
Cf2py  intent(in) nullval
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
        LOGICAL LVAL
C
        IF(ISTAT.NE.0) RETURN
C
        IF(IFITS.EQ.0) THEN
                WRITE(*,*) 'QR_FITSGCVD error - file closed'
                ISTAT=1
                RETURN
        ENDIF
C
        CALL FTGCVD(IFITS,ICOL,IROW,1,NEL,NULLVAL,VAL,LVAL,ISTAT)
        END
