*+AX_DONA2V	Converts spherical-polar coordinate to 3-vector
	SUBROUTINE AX_DONA2V(AZ,EL,V)
	DOUBLE PRECISION AZ, EL, V(3)
*AZ	Input	R.A. or longitude or azimuth, radians.
*EL	Input	Dec or latitude or elevation, radians.
*V	Output	Unit 3-vector.
*-
*Author: Gavin Eadie, 1972;  Fortran-77 version by Clive Page, 1984 November.
	DOUBLE PRECISION CEL
	CEL = COS(EL)
	V(1) = COS(AZ) * CEL
	V(2) = SIN(AZ) * CEL
	V(3) = SIN(EL)
	END
