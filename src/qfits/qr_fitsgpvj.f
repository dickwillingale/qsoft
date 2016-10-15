*+QR_FITSGPVJ        Get FITS primary array as integer
        SUBROUTINE QR_FITSGPVJ(NULLVAL,NEL,IVAL)
        IMPLICIT NONE
        INTEGER NULLVAL,NEL,IVAL(NEL)
*NULLVAL        input        undefined value
Cf2py  intent(in) nullval
*NEL                input        number of elements
Cf2py  intent(in) nel
*IVAL                output        values returned
Cf2py  intent(out) ival
*-Author: Dick Willingale 2012-Dec-27
        INCLUDE 'QR_COM'
        LOGICAL LVAL
C
        IF(ISTAT.NE.0) RETURN
C
        IF(IFITS.EQ.0) THEN
                WRITE(*,*) 'QR_FITSGPVJ error - file closed'
                ISTAT=1
                RETURN
        ENDIF
C
        CALL FTGPVJ(IFITS,0,1,NEL,NULLVAL,IVAL,LVAL,ISTAT)
        END
