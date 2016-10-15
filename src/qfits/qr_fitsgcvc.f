*+QR_FITSGCVC        Get FITS column real complex
        SUBROUTINE QR_FITSGCVC(NULLVAL,ICOL,IROW,NEL,VAL)
        IMPLICIT NONE
        INTEGER ICOL,IROW,NEL
        REAL NULLVAL
        COMPLEX VAL(NEL)
*NULLVAL        input        undefined value
Cf2py  intent(in) nullval
*ICOL                input        column number
Cf2py  intent(in) icol
*IROW                input        1st row number
Cf2py  intent(in) irow
*NEL                input        number of elements
Cf2py  intent(in) nel
*VAL                output        complex values returned
Cf2py  intent(out) val
*-Author: Dick Willingale 2012-Dec-27
        INCLUDE 'QR_COM'
        LOGICAL LVAL
C
        IF(ISTAT.NE.0) RETURN
C
        IF(IFITS.EQ.0) THEN
                WRITE(*,*) 'QR_FITSGCVC error - file closed'
                ISTAT=1
                RETURN
        ENDIF
C
        CALL FTGCVC(IFITS,ICOL,IROW,1,NEL,NULLVAL,VAL,LVAL,ISTAT)
        END
