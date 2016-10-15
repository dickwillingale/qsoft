*+SRT_CPII	Chebyshev Polynomials of the 2nd kind
	FUNCTION SRT_CPII(N,X)
	IMPLICIT NONE
	INTEGER N
	COMPLEX X,SRT_CPII
*N	input	order of polynomial
*X	input	argument of polynomial
*-Author Dick Willingale 1997-Feb-28
	EXTERNAL SRT_ACOS
	COMPLEX SRT_ACOS,TT,AC,ACN
C
	TT=SQRT(1.0-X**2)
	IF(ABS(TT).GT.0.0) THEN
		AC=SRT_ACOS(X)
		ACN=AC*(N+1)
		IF(ABS(AIMAG(ACN)).GT.80.0) THEN
			ACN=ACN*80.0/AIMAG(ACN)
		ENDIF
		SRT_CPII=SIN(ACN)/TT
	ELSE
		SRT_CPII=N+1
	ENDIF
	END
