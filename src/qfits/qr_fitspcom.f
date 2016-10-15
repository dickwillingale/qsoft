*+QR_FITSPCOM        Put comment on current extension
        SUBROUTINE QR_FITSPCOM(NC,COM)
        IMPLICIT NONE
        INTEGER NC
        CHARACTER COM*(NC)
*NC        input        length of comment
Cf2py  intent(in) nc
*COM        input        comments string
Cf2py  intent(in) com
*-Author: Dick Willingale 2013-Feb-09
        INCLUDE 'QR_COM'
C
        IF(ISTAT.NE.0) RETURN
C
        IF(IFITS.EQ.0) THEN
                WRITE(*,*) 'QR_FITSPCOM - FITS file not open'
                ISTAT=1
                RETURN
        ENDIF
        CALL FTPCOM(IFITS,COM,ISTAT)
        END
