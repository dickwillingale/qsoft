*+SRT_RCON      Radial distance of ray from axis of conic
        SUBROUTINE SRT_RCON(DIR,XA,ANM,VER,H,COSA,P)
	IMPLICIT NONE
        DOUBLE PRECISION DIR(3),XA(3),ANM(3),VER(3),H,COSA,P(3)
*DIR    input   direction of ray
*XA     input   origin of ray
*ANM    input   axis of conic
*VER    input   reference point on axis (defines a vertex or reference plane)
*H      output  height of ray origin above reference plane of conic
*COSA   output  cosine of angle between axis and ray
*P      output  quadratic coefficients
* The result is expressed as a quadratic r**2=P(1).x**2+P(2).x+P(3) where x 
* is distance along ray from origin of ray.
*Author Dick Willingale 1996-Nov-19
        DOUBLE PRECISION DR(3),S,COSG,T(3),Q(3),COSB,ASMALL
	PARAMETER (ASMALL=1.E-10)
C Calculate cos angle between axis and ray direction
        CALL SRT_VDOT(ANM,DIR,COSA)
C Calculate direction and distance from vertex to ray origin
        CALL SRT_DIDI(VER,XA,DR,S)
C Find cos of angle between conic axis and direction from vertex
        CALL SRT_VDOT(ANM,DR,COSG)
C Set height of ray origin above plane
        H=S*COSG
C Calculate cosa*cosg*cosb
        CALL SRT_VCRS(ANM,DIR,T)
        CALL SRT_VCRS(ANM,DR,Q)
        CALL SRT_VDOT(T,Q,COSB)
C Set coefficients 
        P(1)=1.0-COSA**2
        P(2)=2.0*S*COSB
        P(3)=(1.0-COSG**2)*S**2
	IF(ABS(P(1)).LT.ASMALL) THEN
C Trap rounding errors in calculating angle between ray and axis
		P(1)=0.0
		IF(COSA.LT.0.0) THEN
			COSA=-1.0
		ELSE
			COSA=1.0
		ENDIF
	ENDIF
        END
