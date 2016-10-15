*+QR_FITSPKEYD        Put keyword for double on current extension
        SUBROUTINE QR_FITSPKEYD(KEY,VAL,COM)
        IMPLICIT NONE
        CHARACTER KEY*(*),COM*(*)
        DOUBLE PRECISION VAL
*KEY        input        keyword name
Cf2py  intent(in) key
*VAL        input        keyword value
Cf2py  intent(in) val
*COM        input        comments string
Cf2py  intent(in) com
*-Author: Dick Willingale 2013-Feb-09
        INCLUDE 'QR_COM'
        INTEGER ITHERE,IKH,IKA
C
        IF(ISTAT.NE.0) RETURN
C
        IF(IFITS.EQ.0) THEN
                WRITE(*,*) 'QR_FITSPKEYD - FITS file not open'
                ISTAT=1
                RETURN
        ENDIF
C Get number of existing keys and number that can be added
        CALL FTGHSP(IFITS,IKH,IKA,ISTAT)
C Check if key word already exists in header
        CALL QR_FITSEXISTS(IFITS,1,IKH,KEY,ITHERE,ISTAT)
        IF(ITHERE.EQ.0) THEN
                CALL FTPKYD(IFITS,KEY,VAL,12,COM,ISTAT)
        ENDIF
        END
