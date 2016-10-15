*+AX_DAYHMS	Converts fractions of a day to hours, mins, secs as char string
	SUBROUTINE AX_DAYHMS(DAY,HMS)
	REAL DAY
	CHARACTER*(*) HMS
*DAY	input	Fraction of a day, should be >=0.0 and <86400.0
*HMS	output	Time of day as 8-char string in form "HH:MM:SS".
*-
*Author	Clive Page	1985 Jan 17
	CHARACTER STRING*8
	ISECS = INT(MOD(DAY,1.0)*86400.0)
	IHOUR = ISECS / 3600
	ISECS = ISECS - 3600 * IHOUR
	IMINS = ISECS / 60
	ISECS = ISECS - 60 * IMINS
	WRITE(STRING,11)IHOUR,IMINS,ISECS
11	FORMAT(I2.2,2(':',I2.2))
	HMS = STRING
	END
