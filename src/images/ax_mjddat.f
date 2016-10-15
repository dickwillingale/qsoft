*+AX_MJDDAT	Converts MJD (D.P.) to calendar date and time (char string).
	SUBROUTINE AX_MJDDAT(DMJD,DAT)
	DOUBLE PRECISION DMJD
	CHARACTER*(*) DAT
*DMJD	input	Modified Julian Date (IAU definition, of course)
*DAT	output	Date and time as 20-char string in form "YYYY-MMM-DD hh:mm:ss"
*-
*Author	Clive Page	1985 Jan 17
*Modified CGP 1985 Sept 2: put month name in 3-letter form.

	CHARACTER HMS*8, STRING*20, MONAME(12)*3
	DATA MONAME/'Jan','Feb','Mar','Apr','May','Jun',
     &   'Jly','Aug','Sep','Oct','Nov','Dec'/
	MJD = INT(DMJD)
	DAY = DMJD - MJD
	CALL AX_MJDYMD(MJD,IYEAR,MONTH,IDATE)
	CALL AX_DAYHMS(DAY,HMS)
	WRITE(STRING,11)IYEAR,MONAME(MONTH),IDATE,HMS
11	FORMAT(I4,'-',A3,'-',I2.2,1X,A8)
	DAT = STRING
	END
