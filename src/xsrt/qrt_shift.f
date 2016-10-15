*+QRT_SHIFT shift position of surface element
        SUBROUTINE QRT_SHIFT(IS,PL)
*IS        input        surface element number
*PL        input        vector shift
Cf2py  intent(in) IS,PL
        IMPLICIT NONE
        INTEGER IS
        DOUBLE PRECISION PL(3)
*-Author Dick Willingale 2012-Jun-28
        INCLUDE 'QR_COM'
C
        IF(ISTAT.NE.0) RETURN
C
        CALL SRT_SHIFT(IS,PL,ISTAT)
        END        
