*+AX_CONGAL	Converts RA,DEC to Galactic longitude and latitude (rads).
	SUBROUTINE AX_CONGAL(RA,DEC,GL,GB)
	REAL RA, DEC, GL, GB
*RA	Input	R.A. in radians, 1950.0 coords.
*DEC	Input	Declination, radians, 1950.0
*GL	Output	Returns new galactic longitude, 0 - 2*PI, radians.
*GB	Output	Returns new galactic latitude.
*-
*Author: Clive Page, 1972.
	REAL CTOG(3,3)
	DATA CTOG /
     &	 -0.06698873942545,-0.87275576585633,-0.48353891462292,
     &	 0.492728466065,-0.4503469580191,0.744584633291,
     &	 -0.86760081115674,-0.18837460170498,0.46019978478119/
	CALL AX_CONVRT(RA,DEC,CTOG,GL,GB)
	END
