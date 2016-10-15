*+QR_FITSPHIS        Put history on current extension
        SUBROUTINE QR_FITSPHIS(NH,HIS)
        IMPLICIT NONE
        INTEGER NH
        CHARACTER HIS*(NH)
*NH        input        length of history string
Cf2py  intent(in) nh
*HIS        input        history string
Cf2py  intent(in) his
*-Author: Dick Willingale 2013-Feb-09
        INCLUDE 'QR_COM'
C
        IF(ISTAT.NE.0) RETURN
C
        IF(IFITS.EQ.0) THEN
                WRITE(*,*) 'QR_FITSPHIS - FITS file not open'
                ISTAT=1
                RETURN
        ENDIF
        CALL FTPHIS(IFITS,HIS,ISTAT)
        END
