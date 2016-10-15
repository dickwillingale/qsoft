*+AX_LAMBERT	Perform Lambert projection in double precision
*This does so-called equatorial aspect of Azimuthal equal-area projection
	SUBROUTINE AX_LAMBERT(AZ,EL,X,Y)
	DOUBLE PRECISION AZ,EL,X,Y
*AZ	input	azimuth radians
*EL	input	elevation radians
*X	output	radians range -PI/2 to PI/2
*Y	output	radians range -PI/2 to PI/2
*-Author Dick Willingale 2001-Oct-25
	DOUBLE PRECISION RANGE,BEAR
C first transform local azimuth and elevation to a bearing and range
C about the centre. i.e. express with centre as the pole
	RANGE=ABS(ACOS(COS(EL)*COS(AZ)))
	IF(RANGE.NE.0.0) THEN
		BEAR=ATAN2(SIN(EL)/SIN(RANGE),TAN(AZ)/TAN(RANGE))
	ELSE
		BEAR=0.0
	ENDIF
C Project range according to Lambert
	RANGE=2.D0*SIN(RANGE/2.D0)
C Construct new local coodinates using bearing and range
	X=RANGE*COS(BEAR)
	Y=RANGE*SIN(BEAR)
	END
