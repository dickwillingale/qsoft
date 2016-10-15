*-SRT_IN2D	Interpolation of 2-d deformation grid
	SUBROUTINE SRT_IN2D(X,Y,NS,NSMAX,NX,NY,XSAM,YSAM,DELZ,GRADX,GRADY,
     +	Z,DZDX,DZDY,ISTAT)
	IMPLICIT NONE
	INTEGER NS,NSMAX,NX,NY
	DOUBLE PRECISION X,Y,XSAM(NX,NSMAX),YSAM(NY,NSMAX)
	DOUBLE PRECISION DELZ(NX,NY,NSMAX),GRADX(NX,NY,NSMAX)
	DOUBLE PRECISION GRADY(NX,NY,NSMAX),Z,DZDX,DZDY
	INTEGER ISTAT
*X	input	X position on surface (axial for conic)
*Y	input	Y position in surface (azimuth for conic)
*NS	input	sub-surface number (shell for Wolter I nest)
*NSMAX	input	maximum number of sub-surfaces held
*NX	input	number of x samples
*NY	input	number of y samples
*XSAM	input	x samples
*YSAM	input	y samples
*DELZ	input	deformation samples
*GRADX	input	x gradient samples
*GRADY	input	y gradiant samples
*Z	output	deformation at X,Y
*DZDX	output	gradient of deformation wrt x at X,Y
*DZDY	output	gradient of deformation wrt y at X,Y
*ISTAT	in/out	returned status
*-Author Dick Willingale 1996-Nov-21
	DOUBLE PRECISION DX,DY
	INTEGER KK,JJ,KPOS,JPOS,IK,IJ
C
	IF(ISTAT.NE.0) RETURN
C
	IF(NS.GT.0.AND.NS.LE.NSMAX) THEN
C Find position in grid
		KK=NINT(NX*(X-XSAM(1,NS))/(XSAM(NX,NS)-XSAM(1,NS)))
		KK=MIN(MAX(KK,1),NX-1)
		IF(X.GE.XSAM(KK,NS)) THEN
			IK=1
		ELSE
			IK=-1
		ENDIF
		KPOS=0
		DO WHILE(KPOS.EQ.0)
		    IF(KK.LT.NX.AND.KK.GT.0) THEN
			IF(X.GE.XSAM(KK,NS).AND.X.LE.XSAM(KK+1,NS)) THEN
				KPOS=KK
			ENDIF
		    ELSEIF(KK.GE.NX) THEN
			KPOS=-1
		    ELSEIF(KK.LE.0) THEN
			KPOS=-1
		    ENDIF
		    KK=KK+IK
		ENDDO			
		JJ=NINT(NY*(Y-YSAM(1,NS))/(YSAM(NY,NS)-YSAM(1,NS)))
		JJ=MIN(MAX(JJ,1),NY-1)
		IF(Y.GE.YSAM(JJ,NS)) THEN
			IJ=1
		ELSE
			IJ=-1
		ENDIF
		JPOS=0
		DO WHILE(JPOS.EQ.0)
		    IF(JJ.LT.NY.AND.JJ.GT.0) THEN
			IF(Y.GE.YSAM(JJ,NS).AND.Y.LE.YSAM(JJ+1,NS)) THEN
				JPOS=JJ
			ENDIF
		    ELSEIF(JJ.GE.NY) THEN
			JPOS=-1
		    ELSEIF(JJ.LE.0) THEN
			JPOS=-1
		    ENDIF
		    JJ=JJ+IJ
		ENDDO
		KK=KPOS
		JJ=JPOS
	ELSE
		WRITE(*,*) 'SRT_IN2D error - internal index',NS,'  out of range'
		ISTAT=1
		RETURN
	ENDIF
	IF(KK.GT.0.AND.JJ.GT.0) THEN
C Find y gradient by linear interpolation between profiles
		DY=YSAM(JJ+1,NS)-YSAM(JJ,NS)
		IF(DY.NE.0.0) THEN
			DY=(Y-YSAM(JJ,NS))/DY
		ELSE
			DY=0.0
		ENDIF
		DZDY=GRADY(KK,JJ,NS)+(GRADY(KK,JJ+1,NS)-GRADY(KK,JJ,NS))*DY
C Interpolate deformation between profiles
		Z=DELZ(KK,JJ,NS)+(DELZ(KK,JJ+1,NS)-DELZ(KK,JJ,NS))*DY
C Interpolate x gradient between profiles
		DZDX=GRADX(KK,JJ,NS)+(GRADX(KK,JJ+1,NS)-GRADX(KK,JJ,NS))*DY
	ELSE
C Missed grid so return zeros
		Z=0.0
		DZDX=0.0
		DZDY=0.0
	ENDIF
	END
