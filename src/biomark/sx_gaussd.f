*+SX_GAUSSD	Returns normal deviate from upper tail area of Gaussian.
	REAL FUNCTION SX_GAUSSD(PROB)
	REAL PROB
*PROB	   input	Upper tail probability, range 0.0 to 0.9999997
*SX_GAUSSD output	Returns normal deviate, +ve if P < 0.5, else -ve.
*-
* Author	Clive Page	1983 July 5
*Standard Fortran-77 version 1984 Sept 25, CGP.
* Algorithm AS111 from J.D.Beasley and S.G.Springer, Appl Stats 26,118 (1977)

	PARAMETER (A0 =   2.50622823884E0, A1 = -18.61500062529E0,
     &		  A2 =  41.39119773534E0, A3 = -25.44106049637E0,
     &		  B1 =  -8.47351093090E0, B2 =  23.08336743743E0,
     &		  B3 = -21.06224101826E0, B4 =   3.13082909833E0,
     &		  C0 =  -2.78718931138E0, C1 =  -2.29796479134E0,
     &		  C2 =   4.85014127135E0, C3 =   2.32121276858E0,
     &		  D1 =   3.54388924762E0, D2 =   1.63706781897E0)

	Q = PROB - 0.5
	IF(ABS(Q) .LE. 0.42) THEN
	    R = Q**2
	    SX_GAUSSD = -Q * (((A3 * R + A2) * R + A1) * R + A0) /
     &             ((((B4 * R + B3) * R + B2) * R + B1) * R + 1.0)
	ELSE
	    R = PROB
	    IF(Q .GT. 0.0) R = 1.0 - PROB
	    R = SQRT(-LOG(R))
	    SX_GAUSSD = (((C3 * R + C2) * R + C1) * R + C0) /
     &        ((D2 * R + D1) * R + 1.0)
	    IF(Q .GT. 0.0) SX_GAUSSD = -SX_GAUSSD
	END IF
	END
