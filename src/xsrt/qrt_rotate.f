*+QRT_ROTATE rotate surface element
        SUBROUTINE QRT_ROTATE(IS,PL,AX,ANGLE)
*IS        input        surface element number
*PL        input        position of rotation centre
*AX        input        rotation axis
*ANGLE     input        rotation angle
Cf2py  intent(in) IS,PL,AX,ANGLE
        IMPLICIT NONE
        INTEGER IS
        DOUBLE PRECISION PL(3),AX(3),ANGLE
*-Author Dick Willingale 2012-Jun-28
        INCLUDE 'QR_COM'
C
        IF(ISTAT.NE.0) RETURN
C
        CALL SRT_ROTATE(IS,PL,AX,ANGLE,ISTAT)
        END        
