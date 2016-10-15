*+AX_DMAT Set up double precision spherical coordinate transformation matrices
	SUBROUTINE AX_DMAT(CEL,ROLL,CTOS,STOC)
	DOUBLE PRECISION CEL(2),ROLL,CTOS(3,3),STOC(3,3)
*CEL	input	celestial pointing direction (radians)
*ROLL	input	roll angle from Cel. North to spacecraft +ve elevation
*		+ve clockwise as look at sky (radians)
*CTOS	output	Celestial to spacecraft transformation matrix
*STOC	output	Spacecraft to Celestial transformation matrix
*-
*Author Dick Willingale 1985-March-29
	DOUBLE PRECISION PI,TWOPI
	PARAMETER ( TWOPI = 6.28318530717957D0 )
	PARAMETER (PI=TWOPI/2.D0)
	DOUBLE PRECISION VPOLE(3),VEQU(3),RLL,DECP,CDUM,RAP
C Force ROLL to range -PI to PI
	RLL=MOD(ROLL,TWOPI)
	IF(RLL.LT.-PI) RLL=RLL+TWOPI
	IF(RLL.GT.PI) RLL=RLL-TWOPI
C Find Celestial position of local north pole (Z-axis)
	DECP=ASIN(COS(RLL)*COS(CEL(2)))
	CDUM=-TAN(DECP)*TAN(CEL(2))
	CDUM=ACOS(MIN(MAX(CDUM,-1.D0),1.D0))
	RAP=CEL(1)+CDUM
	IF(RLL.GT.0.) RAP=CEL(1)-CDUM
C Convert pole and origin to cartesian 3-vectors
	CALL AX_DONA2V(RAP,DECP,VPOLE)
	CALL AX_DONA2V(CEL(1),CEL(2),VEQU)
C Generate transform matrix celestial to space craft
	CALL AX_DONGEN(VPOLE,VEQU,CTOS)
C Find inverse by transposition
	CALL AX_DONMIN(CTOS,STOC)
	RETURN
	END
