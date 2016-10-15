*+AX_UTDYMD	Converts year and day-in-year to month and date.
	SUBROUTINE AX_UTDYMD(IYEAR,IDAY,MONTH,IDATE)
	INTEGER IYEAR,IDAY,MONTH,IDATE
*IYEAR	in	Year number, 4 digits.
*IDAY	in	Day in year, range 1 - 366.
*MONTH	out	Month number, range 1 - 12.
*IDATE	out	Date in month, range 1 - 31.
*-
*Author	Clive Page	1984 Aug 28.

	INTEGER LENGTH(12)
	SAVE LENGTH
	DATA LENGTH /31,28,31, 30,31,30, 31,31,30, 31,30,31/
	IF(MOD(IYEAR,4) .EQ. 0) THEN
		LENGTH(2) = 29
	ELSE
		LENGTH(2) = 28
	END IF
	N = 0
	DO M = 1,12
		N = N + LENGTH(M)
		IF(N .GE. IDAY) GOTO 20
	END DO
20	CONTINUE
	MONTH = M
	IDATE = IDAY - N + LENGTH(M)
	END
