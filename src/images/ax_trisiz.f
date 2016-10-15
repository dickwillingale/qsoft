
*+AX_TRISIZ	Computes area (steradians) of spherical triangle given its sides.
	REAL FUNCTION AX_TRISIZ(A,B,C)
*A,B,C		Input	Sides of the triangle, radians.
*AX_TRISIZ	Returns	Area of the spherical triangle, steradians.
*-
*Author: Clive Page, 1983.
	S = 0.5 * (A + B + C)
	AX_TRISIZ = 2.0 * ASIN(SQRT(SIN(S) * SIN(S-A) * SIN(S-B) * SIN(S-C))
     $   / (2.0 * COS(0.5*A) * COS(0.5*B) * COS(0.5*C)))
	END
