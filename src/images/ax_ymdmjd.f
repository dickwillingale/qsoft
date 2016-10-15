*+AX_YMDMJD	Converts calendar year+month+date to Modified Julian Date.
	SUBROUTINE AX_YMDMJD(IYEAR,MONTH,IDAY,MJD)
	INTEGER IYEAR, MONTH, IDATE, MJD
*IYEAR	in	Year value, 4 digits, range 1901 to 2099.
*MONTH	in	Month value, range 1 - 12.
*IDAY	in	Date in month, range 1 - 31. (May also be 1-366 if MONTH=1)
*MJD	out	Modified Julian Date.
*-
*Author: C G Page, 1982 Jan 15;  Algorithm from "Almanac for Computers 1982".

	MJD = 367 * IYEAR - (7*(IYEAR + ((MONTH+9)/12))/4) +
     &   ((275 * MONTH)/9) + IDAY - 678987
	END
