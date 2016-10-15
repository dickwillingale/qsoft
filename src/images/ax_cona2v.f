*+AX_CONA2V	Converts spherical-polar coordinate to 3-vector
	SUBROUTINE AX_CONA2V(AZ,EL,V)
	REAL AZ, EL, V(3)
*AZ	Input	R.A. (or longitude or azimuth), radians.
*EL	Input	Dec (or latitude or elevation), radians.
*V	Output	Returns 3-vector of unit length.
*-
*Author: Gavin Eadie, 1972;  Fortran-77 version by Clive Page, 1984 November.
	CEL = COS(EL)
	V(1) = COS(AZ) * CEL
	V(2) = SIN(AZ) * CEL
	V(3) = SIN(EL)
	END
