*+QRF_INIT        Initialise/reset the qRfits environment
        SUBROUTINE QRF_INIT()
        IMPLICIT NONE
*ISTAT        in/out        returned status
*-Author Dick Willingale 2016-Sep-28
        INCLUDE 'QR_COM'
C
        ISTAT=0
C HDS 
        ILOC=0
        IPTR=0
        IHDS=0
C FITS
        IFITS=0
        END
