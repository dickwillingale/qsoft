*+SRT_CPSD	Chebyshev Polynomials of the 2nd kind
	FUNCTION SRT_CPSD(N,X)
	IMPLICIT NONE
	INTEGER N
	COMPLEX SRT_CPSD,X
* Generate Chebyshev Polynomial using the recurrence relation.
*-Author Dick Willingale 1997
	COMPLEX X2,UJ,UJ1,UJ2
	INTEGER J
C
	IF(N.EQ.0) THEN
		SRT_CPSD=CMPLX(1.0,0.0)
	ELSEIF(N.EQ.1) THEN
		SRT_CPSD=X*CMPLX(2.0,0.0)
	ELSE
		X2=X*CMPLX(2.0,0.0)
		UJ2=CMPLX(1.0,0.0)
		UJ1=X2
		DO J=2,N
			UJ=UJ1*X2-UJ2
			UJ2=UJ1
			UJ1=UJ
		ENDDO
		SRT_CPSD=UJ
	ENDIF
	END
