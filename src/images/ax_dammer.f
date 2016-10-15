*+AX_DAMMER	Perform Hammer's (Aitoff) projection in double precision
	SUBROUTINE AX_DAMMER(AZ,EL,X,Y)
	DOUBLE PRECISION AZ,EL,X,Y
*AZ	input	azimuth radians
*EL	input	elevation radians
*X	output	radians range -PI to PI
*Y	output	radians range -PI/2 to PI/2
*-Author Dick Willingale 1985-Dec-4
	DOUBLE PRECISION PI,TWOPI,PIBY2
	DOUBLE PRECISION HAZ,T,Z
C
	PIBY2=ASIN(1.D0)
	PI=PIBY2*2.D0
	TWOPI=PI*2.D0
C Find half azimuth in range -PI/2 to +PI/2
	HAZ=(MOD(AZ+PI,TWOPI)-PI)*0.5D0
C Set Z
	Z=SQRT(ABS(1.0D0-COS(HAZ)*COS(EL)))
C Trap origin
	IF(Z.EQ.0.D0) THEN
		T=0.0D0
	ELSE
		T=ATAN2(SIN(HAZ),TAN(EL))
	ENDIF
C Now project
	X=Z*SIN(T)*PI
	Y=Z*COS(T)*PIBY2
	END
