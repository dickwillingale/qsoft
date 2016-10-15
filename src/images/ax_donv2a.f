*+AX_DONV2A	Conversion of 3-vector to spherical polars, radians.
	SUBROUTINE AX_DONV2A(V,AZ,EL)
	DOUBLE PRECISION V(3), AZ, EL
*V	Input 	3-vector
*AZ	Output	Azimuth/RA/longitude, radians 0-2*PI
*EL	Output	Elevation/declination/latitude, radians.
*-
*Author: Gavin Eadie, 1972;  Fortran-77 version by Clive Page, 1984 November.
	DOUBLE PRECISION TWOPI
	PARAMETER ( TWOPI = 6.28318530717957D0 )
	AZ = ATAN2(V(2),V(1))
	IF(AZ .LT. 0.D0) AZ = TWOPI + AZ
	EL = ASIN(V(3))
	END
