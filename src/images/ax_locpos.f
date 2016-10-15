*+AX_LOCPOS calculates 'local' az.and el. of object for given geographical posn
	SUBROUTINE AX_LOCPOS(DMJD,RAS,DECS,GLONG,GLAT,AZ,EL)
	DOUBLE PRECISION DMJD
	REAL RAS,DECS
	REAL GLAT,GLONG
	REAL AZ,EL
*DMJD	in	MJD of observation
*RAS	in	R.A. of object (radians)
*DECS	in	Dec. of object (radians)
*GLONG	in	geographical longitude (radians)
*		[-ve if W.of Greenwich]
*GLAT	in	geographical latitude (radians)
*AZ	out	local azimuth of object (radians)
*		measured from meridian
*EL	out 	local elevation of object (radians)
*		measured from horizon
*
*- author MGW 	12-12-88 	calls: AX_MJDGMST etc.

	PARAMETER (PI=3.141592654, TWOPI=2.0*PI, DTOR=PI/180.0)
	REAL STG, STL
	REAL VP(3),VH(3),CTOAA(3,3)

C sidereal time (greenwich + local)
	CALL AX_MJDGMST(DMJD,STG)
	STL = STG + GLONG
	IF(STL.LT.0.0) STL = STL + TWOPI
	IF(STL.GT.TWOPI) STL = STL - TWOPI
C dec of zenith and horizon
	ZDEC = GLAT
	ZHOR = PI/2.0 - GLAT
C compute matrix- pole at zenith, equator at horizon (0=meridian)
	CALL AX_CONA2V(STL,ZDEC,VP)
	CALL AX_CONA2V(STL,ZHOR,VH)
	CALL AX_CONGEN(VP,VH,CTOAA)
C
	CALL AX_CONVRT(RAS,DECS,CTOAA,AZ,EL)

	END
