*+AX_CONVRT	Converts spherical-polar coordinates to another frame
	SUBROUTINE AX_CONVRT(RA, DEC, CTOE, AZ, EL)
	REAL RA		!input	RA-like coordinate, radians.
	REAL DEC	!input	DEC-like coordinate, radians.
	REAL CTOE(3,3)	!input	3x3 converrsion e.g. from AX_CONGEN
	REAL AZ		!output	Azimuth in new frame, radians 0 to 2pi
	REAL EL		!output	Elevation in new frame, radians range
			!	-pi/2 to +pi/2.
*-Author	Clive Page	1989 Feb 24
	PARAMETER (PI = 3.14159265, TWOPI = 2.0 * PI)
*	CALL MTH$SINCOS(RA, SRA, CRA)
*	CALL MTH$SINCOS(DEC, V3, CDEC)
        SRA=SIN(RA)
        CRA=COS(RA)
        V3=SIN(DEC)
        CDEC=COS(DEC)
	V1 = CRA * CDEC
	V2 = SRA * CDEC
	X1 = V1 * CTOE(1,1) + V2 * CTOE(2,1) + V3 * CTOE(3,1)
	X2 = V1 * CTOE(1,2) + V2 * CTOE(2,2) + V3 * CTOE(3,2)
	X3 = V1 * CTOE(1,3) + V2 * CTOE(2,3) + V3 * CTOE(3,3)
	AZ = ATAN2(X2, X1)
	IF(AZ .LT. 0.0) AZ = AZ + TWOPI
	EL = ASIN(X3)
	END
