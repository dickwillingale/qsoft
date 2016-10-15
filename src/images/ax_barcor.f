*+AX_BARCOR	Computes barycentric/heliocentric coordinates of Earth.
	SUBROUTINE AX_BARCOR(DCORH,DCORB)
	IMPLICIT DOUBLE PRECISION (D)
	DIMENSION DCORH(3),DCORB(3)
*	Input	Common block /BARXYZ/ computed by prior call of AX_BARVEL.
*DCORH	Output	Heliocentric coordinates of Earth
*		in x,y,z directions in A.U.
*DCORB	Output	Barycentric coordinates as above.
*-
*Accuracy: The largest deviations from the JPL-DE96 are
*     0.000011 A.U. (5 millisecs) for the heliocentric,
*     and 0.000046 A.U. (23 millisecs) for the barycentric coordinates.
*Must by preceded by CALL AX_BARVEL to set values in common block.
*
* Written by P.Stumpff	1979 July 15.
* Modified for Standard Fortran by C.G. Page 1983 Sept 27.
*
	DIMENSION CCPAM(4)
	COMMON /BARXYZ/ DPREMA(3,3), DPSI, D1PDRO, DSINLS, DCOSLS,
     &   DSINEP, DCOSEP, FORBEL(7), SORBEL(17), SINLP(4), COSLP(4),
     &   SINLM, COSLM, SIGMA, IDEQ

*   CCPAM(K) = A*M(PLANETS),CCIM = INCLINATION(MOON),DC1MME = 1-MASS(EMB)

	PARAMETER (CCIM = 8.978749E-2, DC1MME = 0.99999696D0)
	DATA CCPAM/4.960906E-3,2.727436E-3,8.392311E-4,1.556861E-3/

*     HELIOCENTRIC COORDINATES OF THE EARTH (BARVL197-232)

	DR    = DPSI*D1PDRO
	FLATM = CCIM*SIN(FORBEL(3))
	A     = SIGMA*COS(FLATM)
	DXH   = DR*DCOSLS - A*COSLM
	DYH   = DR*DSINLS - A*SINLM
	DZH   = -SIGMA*SIN(FLATM)

*     BARYCENTRIC COORDINATES OF THE EARTH (BARVL234-248)

	DXB = DXH*DC1MME
	DYB = DYH*DC1MME
	DZB = DZH*DC1MME
	DO 10, K = 1,4
	    FLAT = SORBEL(K+13)*SIN(FORBEL(K+3)-SORBEL(K+5))
	    A = CCPAM(K)*(1.0-SORBEL(K+9)*COS(FORBEL(K+3)-
     &		SORBEL(K+1)))
	    B = A*COS(FLAT)
	    DXB = DXB-B*COSLP(K)
	    DYB = DYB-B*SINLP(K)
	    DZB = DZB-A*SIN(FLAT)
10  	CONTINUE

*    TRANSITION TO MEAN EQUATOR OF DATE (BARVL250-256)

	DYAH = DCOSEP*DYH-DSINEP*DZH
	DZAH = DSINEP*DYH+DCOSEP*DZH
	DYAB = DCOSEP*DYB-DSINEP*DZB
	DZAB = DSINEP*DYB+DCOSEP*DZB

	IF(IDEQ .EQ. 0) THEN
	    DCORH(1) = DXH
	    DCORH(2) = DYAH
	    DCORH(3) = DZAH
	    DCORB(1) = DXB
	    DCORB(2) = DYAB
	    DCORB(3) = DZAB
	ELSE

*     GENERAL PRECESSION FROM EPOCH DJE TO DEQ (BARVL267-275)

	    DO 30, N = 1,3
	 	DCORH(N) = DXH*DPREMA(N,1)+DYAH*DPREMA(N,2)+
     &			DZAH*DPREMA(N,3)
	 	DCORB(N) = DXB*DPREMA(N,1)+DYAB*DPREMA(N,2)+
     &			DZAB*DPREMA(N,3)
30  	    CONTINUE
	END IF
	END
