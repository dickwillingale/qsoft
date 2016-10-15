*+QRT_WOLTER2 Set Wolter type II
        SUBROUTINE QRT_WOLTER2(RP,GP,RH,GH,RM,FOVR,AX,AR,FF,IDF,IQ)
*RP        input    maximum radius of parabola
*GP        input    grazing angle (degrees) at maximum radius on parabola
*RH        input    maximum radius of hyperbola
*GH        input    grazing angle (degrees) at maximum radius on hyperbola
*RM        input    minimum radius of parabola
*FOVR      input    radius of field of view degrees
*AX        input    direction of axis of telescope
*AR        input    reference direction perpendicular to axis
*FF        input    position of focus of telescope
*IDF       input    deformation index
*IQ        input    surface quality index
Cf2py  intent(in) RP,GP,RH,GH,RM,FOVR,AX,AR,FF,IDF,IQ
        IMPLICIT NONE
        DOUBLE PRECISION RP,GP,RH,GH,RM,FOVR
        DOUBLE PRECISION AX(3),AR(3),FF(3),PI
        INTEGER IDF,IQ
*-Author Dick Willingale 2012-Jun-28
        INTEGER IDEF(2)
        DOUBLE PRECISION PP(16),PY(16)
        INCLUDE 'QR_COM'
C
        IF(ISTAT.NE.0) RETURN
C Generate conic coefficients and limits
        PI=ASIN(1.D0)*2.D0
        CALL SRT_WLT2(RP,GP*PI/180.,RH,GH*PI/180.,RM,FOVR*PI/180.,
     +        AX,AR,FF,PP,PY)
C
        IDEF(1)=IDF
        IDEF(2)=1
        CALL SRT_SETF(0,9,16,PP,IDEF,IQ,-1,-1,ISTAT)
C Set parameters in common for hyperbola
        CALL SRT_SETF(0,9,16,PY,IDEF,IQ,-1,-1,ISTAT)
        END
