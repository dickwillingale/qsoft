*+QRT_BAFFLE  Set cylindrical baffle parameters
        SUBROUTINE QRT_BAFFLE(XMIN,XMAX,RAD,AX,AR,RP,IQ)
*XMIN      input        axial position of base
*XMAX      input        axial position of top
*RAD       input        radius of cylinder
*AX        input        axis direction
*AR        input        reference direction perpendicular to axis
*RP        input        position of vertex
*IQ        input        surface quality index
Cf2py  intent(in) XMIN,XMAX,RAD,AX,AR,RP,IQ
        IMPLICIT NONE
        DOUBLE PRECISION XMIN,XMAX,RAD,AX(3),AR(3),RP(3)
        INTEGER IQ
*-Author Dick Willingale 2012-May-16
        INCLUDE 'SRT_COM'
        INTEGER IDEF(2)
        DOUBLE PRECISION PL(16)
        INCLUDE 'QR_COM'
C
        IF(ISTAT.NE.0) RETURN
        CALL SRT_CONE(XMIN,RAD,XMAX,RAD,AX,AR,RP,PL)
        IDEF(1)=0
        IDEF(2)=0
        CALL SRT_SETF(0,9,16,PL,IDEF,IQ,-2,-1,ISTAT)
        END        
