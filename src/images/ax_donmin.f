*+AX_DONMIN	Conversion matrix inversion (by transpose)
	SUBROUTINE AX_DONMIN(CTOS,STOC)
	DOUBLE PRECISION CTOS(3,3),STOC(3,3)
*CTOS	Input 	Conversion matrix
*STOC	Output	Conversion matrix for inverse Conversion.
*-
*Author: Gavin Eadie, 1972;  Fortran-77 version by Clive Page, 1984 November.
	DO I = 1,3
	  DO J = 1,3
	   STOC(J,I) = CTOS(I,J)
	  END DO
	END DO
	END
