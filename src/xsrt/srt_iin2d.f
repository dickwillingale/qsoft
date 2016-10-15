*-SRT_IIN2D	Indexing of 2-d deformation grid
	SUBROUTINE SRT_IIN2D(IX,IY,NS,NSMAX,NX,NY,DELZ,GRADX,GRADY,
     +	Z,DZDX,DZDY,ISTAT)
	IMPLICIT NONE
	INTEGER IX,IY,NS,NSMAX,NX,NY
	DOUBLE PRECISION DELZ(NX,NY,NSMAX),GRADX(NX,NY,NSMAX)
	DOUBLE PRECISION GRADY(NX,NY,NSMAX),Z,DZDX,DZDY
	INTEGER ISTAT
*IX	input	index of position on surface
*IY	input	index of position on surface
*NS	input	sub-surface number
*NSMAX	input	maximum number of sub-surfaces held
*NX	input	number of x samples
*NY	input	number of y samples
*DELZ	input	deformation samples
*GRADX	input	x gradient samples
*GRADY	input	y gradiant samples
*Z	output	deformation at X,Y
*DZDX	output	gradient of deformation wrt x at X,Y
*DZDY	output	gradient of deformation wrt y at X,Y
*ISTAT	in/out	returned status
*-Author Dick Willingale 2012-Jun-18
	INTEGER KK,JJ
C
	IF(ISTAT.NE.0) RETURN
C
	IF(NS.GT.0.AND.NS.LE.NSMAX) THEN
      		KK=IX
		IF(KK.LT.1.OR.KK.GT.NX) THEN
			KK=0
		ENDIF
		JJ=IY
		IF(JJ.LT.1.OR.JJ.GT.NY) THEN
			JJ=0
		ENDIF
	ELSE
		WRITE(*,*) 'SRT_IN2D error - internal index',NS,'  out of range'
		ISTAT=1
		RETURN
	ENDIF
	IF(KK.GT.0.AND.JJ.GT.0) THEN
		Z=DELZ(KK,JJ,NS)
		DZDX=GRADX(KK,JJ,NS)
		DZDY=GRADY(KK,JJ,NS)
	ELSE
C Missed grid so return zeros
		Z=0.0
		DZDX=0.0
		DZDY=0.0
	ENDIF
	END
