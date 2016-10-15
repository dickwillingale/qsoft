*+QR_RSEED        Set random number seed
        SUBROUTINE QR_RSEED(ISEED)
*ISEED        input        integer seed
        IMPLICIT NONE
        INTEGER ISEED
Cf2py  intent(in) ISEED
*-Author Dick Willingale 2012-May-04
        INCLUDE 'QR_COM'
C
        IF(ISTAT.NE.0) RETURN
C
        CALL SYS_SRAND(ISEED)
        END
