*+AX_MJDUTD	Converts MJD to year and day-in-year.
	SUBROUTINE AX_MJDUTD(MJD,IYEAR,IDAY)
	INTEGER MJD, IYEAR, IDAY
*MJD	in	Modified Julian Date.
*IYEAR	out	Year number, range 1901 - 1999.
*IDAY	out	Day number in year, range 1 - 366.
*restriction: only valid for MJDs in this century.
*-
*Author	Clive Page	1984 Aug 28.

	NYEAR  = (MJD - 15019) / 365.2499
	MJDSYR = (NYEAR * 365.25) + 15018.9
	IYEAR  = NYEAR + 1900
	IDAY   = MJD - MJDSYR
	END
