*+AX_MJDNAME	Returns name of day of the week for any Modified Julian Date
	SUBROUTINE AX_MJDNAME(MJD,NAME)
	INTEGER MJD
	CHARACTER*(*) NAME
*JMJD	in	Modified Julian Date
*NAME	out	Name of the day, e.g. 'Wednesday'
*-
*Author	Clive Page	1984 AUG 23.

	CHARACTER*9 DAY(0:6)
	DATA DAY/'Wednesday','Thursday','Friday','Saturday',
     &	 'Sunday','Monday','Tuesday'/
	NAME = DAY(MOD(MJD,7))
	END
