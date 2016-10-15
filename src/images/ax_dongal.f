*+AX_DONGAL	Converts RA,DEC to Galactic longitude and latitude (rads).
	SUBROUTINE AX_DONGAL(RA,DEC,GL,GB)
	DOUBLE PRECISION RA, DEC, GL, GB
*RA	Input	R.A. in radians, 1950.0 coords.
*DEC	Input	Declination, radians, 1950.0
*GL	Output	New galactic longitude, 0 - 2*PI, radians.
*GB	Output	New galactic latitude.
*-
*Author: Gavin Eadie, 1972;  Fortran-77 version by Clive Page, 1984 November.
	DOUBLE PRECISION CTOG(3,3)
	DATA CTOG /
     $	 -0.06698873942545,-0.87275576585633,-0.48353891462292,
     +	 0.492728466065,-0.4503469580191,0.744584633291,
     +	 -0.86760081115674,-0.18837460170498,0.46019978478119/
	CALL AX_DONVRT(RA,DEC,CTOG,GL,GB)
	END
