*+AX_GTCIRC	Computes great-circle angle between two sperical-polar positions.
	REAL FUNCTION AX_GTCIRC(RA1,DEC1,RA2,DEC2)
	REAL RA1, DEC1, RA2, DEC2
*RA1		Input	R.A. (or longitude) of first point, radians.
*DEC1		Input	Declination (or latitude) of first point, radians.
*RA2		Input	R.A. (or longitude) of second point, radians.
*DEC2		Input	Declination (or latitude) of second point, radians.
*AZ_GTCIRC	Output	Angular distance along great circle, radians.
*-
*Author: Clive Page, 1976.
	AX_GTCIRC = ACOS(COS(DEC1) * COS(DEC2) * COS(RA1-RA2) +
     &		SIN(DEC1) * SIN(DEC2))
	END
