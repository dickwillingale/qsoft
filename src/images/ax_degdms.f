*+AX_DEGDMS	Converts decimal degrees to degrees, minutes, seconds of arc.
	SUBROUTINE AX_DEGDMS(DEGS, DMS)
	REAL DEGS		!input: angle in degrees.
	CHARACTER*(*) DMS	!output: character string (>=9 chars normally).
*-Author	Clive Page	1987 MARCH 3
	CHARACTER SIGN
	IF(DEGS .GE. 0.0) THEN
		SIGN = ' '
	ELSE
		SIGN = '-'
	END IF
	NSECS = NINT(ABS(DEGS) * 3600.0)
	ISECS = MOD(NSECS, 60)
	MINS  = NSECS / 60
	IMINS = MOD(MINS, 60)
	IDEGS = MINS / 60
	WRITE(DMS, 15, IOSTAT=K) SIGN, IDEGS, IMINS, ISECS
15	FORMAT(A1, I2.2, 2(':',I2.2))
	END
