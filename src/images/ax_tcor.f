*+AX_TCOR provides barycentric and heliocentric corrections
	SUBROUTINE AX_TCOR(DMJD,RA,DEC,SECH,SECB)
	IMPLICIT DOUBLE PRECISION (D)
	DOUBLE PRECISION DMJD
	REAL RA,DEC,SECH,SECB
*DMJD      Input  MJD
*RA        Input  source position (radians)
*DEC       Input  source position (radians)
*SECH      Output heliocentric time difference (sec) ADD to DMJD
*SECB	   Output barycentric time difference (sec)
*-
C Calls: AX_BARCOR, AX_BARVEL, AX_DONA2V
C
C Author:   MGW 1987 SEP 22
C
	DIMENSION DVELH(3),DVELB(3),DCORH(3),DCORB(3),DVA(3),DVB(3)

	DATA DEQ/1950.0D0/
	DATA AUSEC/499.01265/

	CALL AX_BARVEL(DMJD,DEQ,DVELH,DVELB)
	CALL AX_BARCOR(DCORH,DCORB)
* sourcce position vector
	CALL AX_DONA2V(DBLE(RA),DBLE(DEC),DVB)
* heliocentric
	DO I=1,3
		DVA(I) = DCORH(I)
	ENDDO
	SECH = AUSEC*(DVA(1)*DVB(1)+DVA(2)*DVB(2)+DVA(3)*DVB(3))
* barycentric
	DO I=1,3
		DVA(I) = DCORB(I)
	ENDDO
	SECB = AUSEC*(DVA(1)*DVB(1)+DVA(2)*DVB(2)+DVA(3)*DVB(3))
	END
