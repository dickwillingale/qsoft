*+AX_RA2HMS	Converts RA (radians) to string HH MM SS.S with specified separator
	SUBROUTINE AX_RA2HMS(RA, SEP, STRING)
	REAL RA			!input	Right Ascension, radians
	CHARACTER SEP*2		!input	Chars to go between H M and S.
	CHARACTER STRING*(*)	!output	String HHxMMySS.S with xy separators.
*-Author	Clive Page	1990-Sept-20
	REAL PI, TWOPI
	PARAMETER (PI = 3.14159265, TWOPI = 2.0*PI)
	INTEGER NDSECS
	CHARACTER TEMP*10
*Compute number of deciseconds
	NDSECS = NINT(MOD(RA+TWOPI,TWOPI) * 432000.0 / PI)
	WRITE(TEMP, 5) NDSECS/36000, SEP(1:1), MOD(NDSECS/600,60),
     &   SEP(2:2), MOD(NDSECS/10,60), MOD(NDSECS,10)
5	FORMAT(I2.2, A1, I2.2, A1, I2.2, '.', I1.1)
	STRING = TEMP
	END
