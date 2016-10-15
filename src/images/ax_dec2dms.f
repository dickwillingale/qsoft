*+AX_DEC2DMS	Converts DEC (radians) to string DD MM SS with specified separator
	SUBROUTINE AX_DEC2DMS(DEC, SEP, STRING)
	REAL DEC		!input	Declination, radians
	CHARACTER SEP*2		!input	Chars to go between D M and S.
	CHARACTER STRING*(*)	!output	String DDxMMySS.S with xy separators.
*-Author	Clive Page	1990-Sept-20
	REAL PI, HALFPI, ADEC
	PARAMETER (PI = 3.14159265, HALFPI = PI/2.0)
	INTEGER NSECS
	CHARACTER SIGN*1, TEMP*10
*Compute number of seconds
	IF(DEC .GE. 0.0) THEN
	    SIGN = '+'
	ELSE
	    SIGN = '-'
	END IF
	NSECS = NINT(MIN(99.9, ABS(DEC)*180.0/PI)*3600.0)
	WRITE(TEMP, 5) SIGN, NSECS/3600, SEP(1:1), MOD(NSECS/60,60),
     &   SEP(2:2), MOD(NSECS,60)
5	FORMAT(A1, I2.2, A1, I2.2, A1, I2.2)
	STRING = TEMP
	END
