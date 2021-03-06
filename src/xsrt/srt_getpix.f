*-SRT_GETPIX	Get position and pixel value from 2-d deformation grid
	SUBROUTINE SRT_GETPIX(IX,IY,NS,NSMAX,NX,NY,XX,YY,DELZ,X,Y,Z,ISTAT)
	IMPLICIT NONE
	INTEGER IX,IY,NS,NSMAX,NX,NY
	DOUBLE PRECISION XX(NX),YY(NY),DELZ(NX,NY,NSMAX),X,Y,Z
	INTEGER ISTAT
*IX,IY	input	pixel indices
*NS	input	sub-surface number (shell for Wolter I nest)
*NSMAX	input	maximum number of sub-surfaces held
*NX	input	number of x samples
*NY	input	number of y samples
*XX	input	x samples
*YY	input	y samples
*DELZ	input	deformation samples
*X	output	x position at IX
*Y	output	y position at IY
*Z	output	deformation IX,IY
*ISTAT	in/out	returned status
*-Author Dick Willingale 2004-Mar-18
C
	IF(ISTAT.NE.0) RETURN
      write(*,*) "srt_getpix",ix,iy,ns,nsmax,nx,ny
C
	IF(NS.GT.0.AND.NS.LE.NSMAX) THEN
		IF(IX.GT.0.AND.IX.LE.NX.AND.IY.GT.0.AND.IY.LE.NY) THEN
			X=XX(IX)
			Y=YY(IY)
			Z=DELZ(IX,IY,NS)
		ELSE
			WRITE(*,*) 'SRT_GETPIX - pixel indices out of range'
			ISTAT=1
			RETURN
		ENDIF
	ELSE
		WRITE(*,*) 'SRT_GETPIX - sub-matrix index out of range'
		ISTAT=1
		RETURN
	ENDIF
	END
