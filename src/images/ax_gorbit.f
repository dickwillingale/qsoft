*+AX_GORBIT Predicts geocentric position/velocity for satellite in Earth orbit.
	SUBROUTINE AX_GORBIT(ORBEL,DMJD,POSN,VELY)
	IMPLICIT DOUBLE PRECISION (A-H, O-Z)
	DOUBLE PRECISION ORBEL(6), POSN(3), VELY(3)
*
*ORBEL	Input	Array of orbital elements as follows:
*  ORBEL(1) = semi-major axis, km.
*  ORBEL(2) = eccentricity.
*  ORBEL(3) = orbital inclination to celestial equator, radians.
*  ORBEL(4) = Right Ascension of ascending node, radians.
*  ORBEL(5) = Longitude of perigee, radians.
*  ORBEL(6) = Epoch of elements and time of perigee passage, days [e.g. MJD]
*DMJD	Input	Date for which predicion required, days
*		[e.g. MJD]   Note: must be later than ORBEL(6).
*POSN	Output	Geocentric X,Y,Z coordinates, km.
*			X-axis towards first point of Aries,
*			Z-axis towards celestial pole,
*			Y-axis at right-angles to both, right-handed set.
*VELY	Output	Velocity components along (X,Y,Z), km/sec.
*-
*Author: C G Page	1984-MAY-24
*Algorithm from A.E.Roy, "Orbital Motion" pp80-87,  285, 471.
*

*GMU is constant of Earth's Gravitation,
*RE  is Radius of the Earth,
*CJ2 is first coefficient of Earth's oblateness.
*EPS is convergence criterion for iterative solution of Kepler's equation

	PARAMETER (GMU = 398603.2D0, RE = 6378.165D0,
     $   CJ2 = 1.08263D-3, EPS = 1D-12, PI = 3.14159265358979D0,
     $   TWOPI = 2.0D0 * PI, PIBY2 = 0.5D0 * PI, PI3BY2 = 1.5D0 * PI)

*local storage
	DIMENSION VONE(3), VTWO(3), CTOS(3,3), STOC(3,3)

*RATEZ = mean motion, rads/sec
*RATE  = mean motion corrected for oblateness, rads/sec
*DELTAT = time after perigee (and epoch) secs.

	RATEZ = SQRT(GMU/ORBEL(1)**3)
	E2F   = 1.0 - ORBEL(2)**2
	OBF   = 1.5 * CJ2 * (RE / (ORBEL(1) * E2F))**2
	S2I   = SIN(ORBEL(3))**2
	RATE  = RATEZ * (1D0 + OBF * (1D0 - 1.5D0 * S2I) * SQRT(E2F))
	DELTAT = (DMJD - ORBEL(6)) * 86400D0

*AMEAN = mean anomaly, radians
*AECC  = eccentric anomaly, radians
*ATRUE = true anomaly, radians.
*Solve Keplers equation by iteration, seems to take 3 - 7 cycles

	AMEAN = MOD(RATE * DELTAT, TWOPI)
	IF(AMEAN .LT. 0D0) AMEAN = AMEAN + TWOPI
	AECC  = AMEAN + ORBEL(2) * SIN(AMEAN)
	DO I = 1,10
	  DE = (AMEAN - AECC + ORBEL(2) * SIN(AECC)) /
     $                              (1D0 - ORBEL(2) * COS(AECC))
	  AECC = AECC + DE
	  IF(ABS(DE) .LT. EPS) GOTO 20
15	END DO
20	CONTINUE
	EROOT = SQRT( (1D0 + ORBEL(2)) / (1D0 - ORBEL(2)) )
	ATRUE = 2D0 * ATAN(EROOT * TAN(0.5D0 * AECC) )
	IF(ATRUE .LT. 0D0) ATRUE = ATRUE + TWOPI

*R = radius vector, km
*V = total velocity, km/sec

	R = ORBEL(1) * (1D0 - ORBEL(2) * COS(AECC))
	V = SQRT(GMU * ( (2D0/R) - (1D0/ORBEL(1))) )

*ASCNOD = precessed value of RA of ascending node
*PERIGE = precessed value of longitude of perigee

	OBDELT = OBF * RATE * DELTAT
	ASCNOD = ORBEL(4) - OBDELT * COS(ORBEL(3))
	PERIGE = ORBEL(5) + OBDELT * (2D0 - 2.5D0 * S2I)

*PHI = angle between velocity and radius vectors,
*PSI = angle between velocity vector and ascending node

	PHI = ACOS(-ORBEL(2) * SIN(AECC) /
     $                        SQRT(1D0 - (ORBEL(2) * COS(AECC))**2))
	PSI = ATRUE + (PI - PHI) + PERIGE

*Pole of orbit is at RA = ASCNOD + 1.5*PI, DEC = PI/2 - inc
*Zero point of reference is RA = ASCNOD, DEC = 0
*Generate coordinate conversion matrices between orbit plane and celestials.

	CALL AX_DONA2V(ASCNOD+PI3BY2, PIBY2-ORBEL(3), VONE)
	CALL AX_DONA2V(ASCNOD, 0D0, VTWO)
	CALL AX_DONGEN(VONE,VTWO,CTOS)
	CALL AX_DONMIN(CTOS,STOC)
	CALL AX_DONVRT(PERIGE+ATRUE, 0D0, STOC, RA, DEC)
	CALL AX_DONA2V(RA, DEC, VONE)
	CALL AX_DONVRT(PSI, 0D0, STOC, VRA, VDEC)
	CALL AX_DONA2V(VRA, VDEC, VTWO)
	DO 35,I = 1,3
	  POSN(I) = VONE(I) * R
	  VELY(I) = VTWO(I) * V
35	CONTINUE
	END
