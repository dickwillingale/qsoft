*+AX_MJDGMST	Computes Greenwich Mean Sidereal Time for any date 1950 - 1999.
	SUBROUTINE AX_MJDGMST(ZMJD,GMST)
	IMPLICIT DOUBLE PRECISION (Z)
	DOUBLE PRECISION ZMJD
	REAL GMST
*ZMJD	in	Modified Julian Date.
*GMST	out	Greenwich Mean Sidereal Time
*		i.e. R.A. of Greenwich meridian, radians.
*-
*Author	Clive Page	1984-JUN-17.

	ZDAYS = ZMJD - 33282D0
	ZINT  = AINT(ZDAYS)
	ZMST  = 1.7466477033819D0 + 1.72027914861503D-2 * ZINT +
     &		6.30038809866574D0 * (ZDAYS - ZINT)
	GMST  = MOD(ZMST,6.28318530717957D0)
	END
