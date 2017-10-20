*+QRB_INIT        Initialise/reset the BIOMARK environment
        SUBROUTINE QRB_INIT()
        IMPLICIT NONE
*ISTAT        in/out        returned status
*-Author Dick Willingale 2017-Jun-26
        CALL QR_INIT()
        END
*+QR_INIT        Initialise/reset the QR environment
        SUBROUTINE QR_INIT()
        IMPLICIT NONE
*ISTAT        in/out        returned status
*-Author Dick Willingale 2012-May-04
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
