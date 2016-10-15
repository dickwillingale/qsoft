*+AX_CONV2A	Conversion of 3-vector to spherical polars, radians.
	SUBROUTINE AX_CONV2A(V,AZ,EL)
	REAL V(3)
*V	Input 	3-vector
*AZ	Output	Azimuth/RA/longitude, radians 0 - 2*PI
*EL	Output	Elevation/declination/latitude, radians.
*-
*Author: Gavin Eadie, 1972;  Fortran-77 version by Clive Page, 1984 November.
	PARAMETER ( TWOPI = 6.283185308 )
	AZ = ATAN2(V(2),V(1))
	IF(AZ .LT. 0.0) AZ = TWOPI + AZ
	EL = ASIN(V(3))
	END
