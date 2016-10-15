*+AX_CONGEN	Generates coordinate Conversion matrix from two 3-vectors
	SUBROUTINE AX_CONGEN(VZ,VX,CTOS)
	REAL VZ(3), VX(3), CTOS(3,3)
*VZ	Input	Vector of pole of new coords in old frame
*VX	Input	Vector  of reference point on equator (ditto)
*CTOS	Output	Unit matrix for Conversion from old to new coordinates.
*-
*Author: Gavin Eadie, 1972;  Fortran-77 version by Clive Page, 1984 November.
	DO I = 1,3
	  CTOS(I,3) = VZ(I)
	END DO
	CALL AX_CONVCP(VZ,VX,CTOS(1,2))
	CALL AX_CONVCP(CTOS(1,2),VZ,CTOS(1,1))
	END
