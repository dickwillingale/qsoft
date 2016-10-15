*+AX_DTRANS	Transforms celestial to other frame or vice-versa
	SUBROUTINE AX_DTRANS(RA, DEC, CTOS, FWD, AZ, EL)
	DOUBLE PRECISION RA		!input	Radians
	DOUBLE PRECISION DEC		!input	Radians
	DOUBLE PRECISION CTOS(3,3)	!input	Celestial - other matrix
	LOGICAL FWD			!input	.true. for cel->other
	DOUBLE PRECISION AZ		!output	Radians
	DOUBLE PRECISION EL		!output	Radians
*-Author	Clive Page	1988 Sept 27
	DOUBLE PRECISION V1, V2, V3, V(3), CRA, CDEC, SRA, TWOPI
	PARAMETER (TWOPI = 6.28318530717957D0)

*	CALL MTH$DSINCOS(RA, SRA, CRA)
*	CALL MTH$DSINCOS(DEC, V3, CDEC)
        SRA=SIN(RA)
        CRA=COS(RA)
        V3=SIN(DEC)
        CDEC=COS(DEC)
	V1 = CRA * CDEC
	V2 = SRA * CDEC
	IF(FWD) THEN
	    DO J = 1,3
		V(J) = V1 * CTOS(1,J) + V2 * CTOS(2,J) + V3 * CTOS(3,J)
	    END DO
	ELSE
	    DO J = 1,3
		V(J) = V1 * CTOS(J,1) + V2 * CTOS(J,2) + V3 * CTOS(J,3)
	    END DO
	END IF
	AZ = ATAN2(V(2), V(1))
	IF(AZ .LT. 0.0D0) AZ = AZ + TWOPI
	EL = ASIN(V(3))
	END
