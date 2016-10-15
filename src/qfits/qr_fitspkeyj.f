*+QR_FITSPKEYJ        Put keyword for integer on current extension
        SUBROUTINE QR_FITSPKEYJ(KEY,IVAL,COM)
        IMPLICIT NONE
        INTEGER IVAL
        CHARACTER KEY*(*),COM*(*)
*KEY        input        keyword name
Cf2py  intent(in) key
*IVAL        input        keyword value
Cf2py  intent(in) ival
*COM        input        comments string
Cf2py  intent(in) com
*-Author: Dick Willingale 2013-Feb-09
        INCLUDE 'QR_COM'
        INTEGER ITHERE,IKH,IKA
C
        IF(ISTAT.NE.0) RETURN
C
        IF(IFITS.EQ.0) THEN
                WRITE(*,*) 'QR_FITSPKEYJ - FITS file not open'
                ISTAT=1
                RETURN
        ENDIF
C Get number of existing keys and number that can be added
        CALL FTGHSP(IFITS,IKH,IKA,ISTAT)
C Check if key word already exists in header
        CALL QR_FITSEXISTS(IFITS,1,IKH,KEY,ITHERE,ISTAT)
        IF(ITHERE.EQ.0) THEN
                CALL FTPKYJ(IFITS,KEY,IVAL,COM,ISTAT)
        ENDIF
        END
