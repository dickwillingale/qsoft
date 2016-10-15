*+AX_DONGEN	Generates coordinate Conversion matrix from two 3-vectors
	SUBROUTINE AX_DONGEN(VZ,VX,CTOS)
	DOUBLE PRECISION VZ(3), VX(3), CTOS(3,3)
*VZ	Input	Vector of pole of new coords in old frame
*VX	Input	Vector  of reference point on equator (ditto)
*CTOS	Output	Unit matrix for conversion to new coordinates.
*-
*Author: Gavin Eadie, 1972;  Fortran-77 version by Clive Page, 1984 November.
	DO I = 1,3
	  CTOS(I,3) = VZ(I)
	END DO
	CALL AX_DONVCP(VZ,VX,CTOS(1,2))
	CALL AX_DONVCP(CTOS(1,2),VZ,CTOS(1,1))
	END
