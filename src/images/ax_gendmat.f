*+AX_GENDMAT Generates DP transformation matrices given RA, DEC, rollangle.
	SUBROUTINE AX_GENDMAT(DRA, DEC, DROLL, DCTOS, DSTOC)
	DOUBLE PRECISION DRA		!input	Right Ascension (radians).
	DOUBLE PRECISION DEC		!input	Declination (radians).
	DOUBLE PRECISION DROLL		!input	Position angle of y-axis
					!	anticlockwise N thro E (rads)
	DOUBLE PRECISION DCTOS(3,3)	!output	Matrix celestial to spacecraft
	DOUBLE PRECISION DSTOC(3,3)	!output	Matrix spacecraft to celestial.
*-Author	Clive Page	1988 Sept 20
*After AX_DMAT of Dick Willingale.
	DOUBLE PRECISION PI, TWOPI
	PARAMETER (PI = 3.141592653589793D0, TWOPI = 2.0D0 * PI)
	DOUBLE PRECISION VPOL(3), VEQU(3), ROLL, DECPOL, CDUM, RAPOL
* Force ROLL to range -PI to PI
	ROLL = MOD(DROLL+TWOPI+PI,TWOPI) - PI
* Find Celestial position of local north pole (Z-axis)
	DECPOL = ASIN(COS(ROLL) * COS(DEC))
	CDUM   = ACOS(MIN(MAX(-1.0D0,-TAN(DECPOL) * TAN(DEC)),1.0D0))
	RAPOL  = DRA + CDUM
	IF(ROLL .GT. 0.0D0) RAPOL = DRA - CDUM
* Convert pole and origin to cartesian 3-vectors
	CALL AX_DONA2V(RAPOL, DECPOL, VPOL)
	CALL AX_DONA2V(DRA,   DEC,    VEQU)
* Generate transform matrix celestial to space craft and inverse
	CALL AX_DONGEN(VPOL, VEQU, DCTOS)
	CALL AX_DONMIN(DCTOS, DSTOC)
	END
