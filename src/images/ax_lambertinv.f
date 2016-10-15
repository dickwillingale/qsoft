*+AX_LAMBERTINV	Perform inverse Lambert projection in double precision
*inverse of  so-called equatorial aspect of Azimuthal equal-area projection
	SUBROUTINE AX_LAMBERTINV(X,Y,AZ,EL,STATUS)
	DOUBLE PRECISION X,Y,AZ,EL
      	INTEGER STATUS
*X	input	radians range -PI/2 to PI/2
*Y	input	radians range -PI/2 to PI/2
*AZ	output	azimuth radians
*EL	output	elevation radians
*STATUS	output	0 if OK, 1 if X,Y out of range
*-Author Dick Willingale 2012-May-12
	DOUBLE PRECISION RANGE,BEAR,COSEL,PIBY2
C Check 
        PIBY2=ASIN(1.D0)
	IF(ABS(X).GT.PIBY2.OR.ABS(Y).GT.PIBY2) THEN
		STATUS=1
		RETURN
	ENDIF
	STATUS=0
C Construct range and bearing
      	RANGE=SQRT(X**2+Y**2)
	BEAR=ATAN2(Y,X)
C Inverse project according to Lambert
	RANGE=ASIN(RANGE/2.D0)*2.D0
C Construct local azimuth and elevation
	EL=ASIN(SIN(BEAR)*SIN(RANGE))
	COSEL=COS(EL)
	IF(COSEL.NE.0.0) THEN
		AZ=ACOS(COS(RANGE)/COS(EL))
	ELSE
		AZ=0.0
	ENDIF
	END
