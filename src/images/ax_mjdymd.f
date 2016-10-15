*+AX_MJDYMD	Converts MJD to calendar year, month, date.
	SUBROUTINE AX_MJDYMD(MJD,IYEAR,MONTH,IDATE)
	INTEGER MJD, IYEAR, MONTH,IDATE
	EXTERNAL AX_MJDUTD, AX_UTDYMD
*MJD	in	Modified Julian Date.
*IYEAR	out	Year number, range 1901 - 1999.
*MONTH	out	Month in year, range 1 - 12.
*IDATE	out	Date in month, range 1 - 31.
*-
*Author	Clive Page	1984 Aug 28.
	CALL AX_MJDUTD(MJD,IYEAR,IDAY)
	CALL AX_UTDYMD(IYEAR,IDAY,MONTH,IDATE)
	END
