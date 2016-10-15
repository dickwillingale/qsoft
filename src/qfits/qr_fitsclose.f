*+QR_FITSCLOSE        Close FITS file
        SUBROUTINE QR_FITSCLOSE()
*-Author: Dick Willingale 2012-Jul-25
        IMPLICIT NONE
        INCLUDE 'QR_COM'
C
        IF(ISTAT.NE.0) THEN
                WRITE(*,*) 'QR_FITSCLOSE status error:',ISTAT
                RETURN
        ENDIF
C
        IF(IFITS.EQ.0) THEN
                WRITE(*,*) 'QR_FITSCLOSE error - file not open'
                ISTAT=1
                RETURN
        ENDIF
        CALL FTCLOS(IFITS,ISTAT)
        IFITS=0
        END
