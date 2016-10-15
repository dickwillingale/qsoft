*+AX_DEG2DMS	Converts decimal degrees to ddmmss as char string.
	SUBROUTINE AX_DEG2DMS(DEGS, DMS)
	REAL DEGS		!input: angle in degrees.
	CHARACTER*(*) DMS	!output: character string (>=7 chars normally).
*-Author	Clive Page	1986 Nov 7
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
	WRITE(DMS, '(A1, 3I2.2)', IOSTAT=K) SIGN, IDEGS, IMINS, ISECS
	END
