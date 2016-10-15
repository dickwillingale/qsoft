*+AX_BOXSIZ	Computes area of error-box composed of spherical quadrilateral.
	REAL FUNCTION AX_BOXSIZ(EBOX)
	REAL EBOX(2,4)
*EBOX	Input	RA,DEC of 4 corners taken cyclically, degrees.
*BOXSIZ	Output	Area of error-box, square degrees.
*-
*
* Author	Clive Page	1983-JULY-20.
*
	REAL CD(4),SD(4),RA(4)
	PARAMETER (DTOR = 1.745329252E-2, SRTOSD = 3.282806350E3)
	DIST(I,J) = ACOS(SD(I) * SD(J) + CD(I) * CD(J) *
     &    COS(RA(I) - RA(J)))
* CD,SD hold cosine and sine of the four declinations, RA holds value in rads.
	DO I = 1,4
	  CD(I) = COS(EBOX(2,I) * DTOR)
	  SD(I) = SIN(EBOX(2,I) * DTOR)
	  RA(I) = EBOX(1,I) * DTOR
	END DO
	D24 = DIST(2,4)
	AX_BOXSIZ = (AX_TRISIZ(DIST(1,2),D24,DIST(4,1)) +
     &            AX_TRISIZ(DIST(2,3),D24,DIST(3,4))) * SRTOSD
	END
