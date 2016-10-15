*+AX_SUNRAD	Computes solar RA and DEC for any MJD in years 1975 to 1999.
	SUBROUTINE AX_SUNRAD(RMJD,SRA,SDEC)
	REAL RMJD, SRA, SDEC
*RMJD	Input	Modified Julian Date
*SRA	Output	Solar R.A. in radians, 1950.0 coords, range 0-2*PI
*SDEC	Output	Solar declination, radians, 1950.0 coords.
*-
* Algorithm from "Almanac for Computers", 1982.
* Epoch changed to 1980 Jan 0 = MJD 44238 to improve accuracy on PDP-11/VAX.
* Results accurate to around 1 arcmin.
* Author: Clive Page, 1982-Sept-24.
	PARAMETER (SINEPS=0.3978045, COSEPS=0.9174702, TWOPI=6.283185308)
	DAYS = RMJD - 44238.0
	RM = 6.21751904 + 0.017201970 * DAYS
	RL = 4.86646193 + 0.017202791 * DAYS + 0.033440 * SIN(RM) +
     &		3.49E-04 * SIN(2*RM)
	SINRL = SIN(RL)
	SRA  = AMOD(ATAN2(COSEPS*SINRL,COS(RL)) + TWOPI, TWOPI)
	SDEC = ASIN(SINEPS * SINRL)
	END
