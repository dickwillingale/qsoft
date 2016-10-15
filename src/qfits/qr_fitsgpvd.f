*+QR_FITSGPVD        Get FITS primary array as real (double)
        SUBROUTINE QR_FITSGPVD(NULLVAL,NEL,VAL)
        IMPLICIT NONE
        INTEGER NEL
        DOUBLE PRECISION NULLVAL,VAL(NEL)
*NULLVAL        input        undefined value
Cf2py  intent(in) nullval
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
                WRITE(*,*) 'QR_FITSGPVD error - file closed'
                ISTAT=1
                RETURN
        ENDIF
C
        CALL FTGPVD(IFITS,0,1,NEL,NULLVAL,VAL,LVAL,ISTAT)
        END
