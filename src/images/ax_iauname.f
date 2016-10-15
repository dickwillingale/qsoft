*+AX_IAUNAME	Converts source RA, DEC (degs) to standard IAU form 'hhmm+ddd'
	SUBROUTINE AX_IAUNAME(RA,DEC,NAME)
	REAL RA,DEC
	CHARACTER NAME*8
*RA	in	Right ascension, degrees (range 0 to 360)
*DEC	in	Declination, degrees (range -90 to +90)
*NAME	out	Name of source in form 'hhmm+ddd' by truncation of coordinates.
*-
*Author	Clive Page	1982 Aug 3.
	IRA = RA * 4.0
	WRITE(NAME,11)IRA/60, MOD(IRA,60), INT(DEC*10)
11	FORMAT(2I2.2,SP,I4.3)
	END
