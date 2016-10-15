*+QR_FITSPKEYL        Put keyword for logical on current extension
        SUBROUTINE QR_FITSPKEYL(KEY,IVAL,COM)
        IMPLICIT NONE
        INTEGER IVAL
        CHARACTER KEY*(*),COM*(*)
*KEY        input        keyword name
Cf2py  intent(in) key
*IVAL        input        keyword value passed as integer
Cf2py  intent(in) ival
*COM        input        comments string
Cf2py  intent(in) com
*-Author: Dick Willingale 2013-Feb-09
        INCLUDE 'QR_COM'
        INTEGER ITHERE,IKH,IKA
        LOGICAL LVAL
C
        IF(ISTAT.NE.0) RETURN
C
        IF(IFITS.EQ.0) THEN
                WRITE(*,*) 'QR_FITSPKEYL - FITS file not open'
                ISTAT=1
                RETURN
        ENDIF
C Select logical value
        IF(IVAL.EQ.0) THEN
                LVAL=.FALSE.
        ELSE
                LVAL=.TRUE.
        ENDIF
C Get number of existing keys and number that can be added
        CALL FTGHSP(IFITS,IKH,IKA,ISTAT)
C Check if key word already exists in header
        CALL QR_FITSEXISTS(IFITS,1,IKH,KEY,ITHERE,ISTAT)
        IF(ITHERE.EQ.0) THEN
                CALL FTPKYL(IFITS,KEY,LVAL,COM,ISTAT)
        ENDIF
        END
