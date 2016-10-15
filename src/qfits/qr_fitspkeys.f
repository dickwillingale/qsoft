*+QR_FITSPKEYS        Put keyword for string on current extension
        SUBROUTINE QR_FITSPKEYS(KEY,STR,COM)
        IMPLICIT NONE
        CHARACTER KEY*(*),STR*(*),COM*(*)
*KEY        input        keyword name
Cf2py  intent(in) key
*STL        input        keyword string value
Cf2py  intent(in) stl
*COM        input        comments string
Cf2py  intent(in) com
*-Author: Dick Willingale 2013-Feb-09
              INCLUDE 'QR_COM'
        INTEGER ITHERE,IKH,IKA
C
        IF(ISTAT.NE.0) RETURN
C
        IF(IFITS.EQ.0) THEN
                WRITE(*,*) 'QR_FITSPKEYS - FITS file not open'
                ISTAT=1
                RETURN
        ENDIF
C Get number of existing keys and number that can be added
        CALL FTGHSP(IFITS,IKH,IKA,ISTAT)
C Check if key word already exists in header
        CALL QR_FITSEXISTS(IFITS,1,IKH,KEY,ITHERE,ISTAT)
        IF(ITHERE.EQ.0) THEN
                CALL FTPKYS(IFITS,KEY,STR,COM,ISTAT)
        ENDIF
        END
