*+SRT_NXTSHFT	Get source shifts for next ray from deformation matrix
	SUBROUTINE SRT_NXTSHFT(IDEF,DX,DY,ISTAT)
	IMPLICIT NONE
	INTEGER IDEF,ISTAT
	DOUBLE PRECISION DX,DY
*IDEF	input	deformation index (used to convert point source into pixels)
*DX	output	shift along source aperture reference axis
*DY	output	shift along source aperture other axis
*ISTAT	in/out	returned status
*-Author Dick Willingale 2004-Mar-18
	INCLUDE 'SRT_COM'
	INTEGER NDONE,IAT,IHERE,NHERE,IX,IY,NPIX
	DOUBLE PRECISION SCALE,X,Y,Z
	SAVE NDONE,IAT,IHERE,NHERE,IX,IY,NPIX
	SAVE SCALE,X,Y
	DATA NDONE/0/
C	
	IF(ISTAT.NE.0) RETURN
C Check deformation index
	IF(IDEF.EQ.0) THEN
		DX=0.0
		DY=0.0
		RETURN
	ENDIF
	IF(NDONE.EQ.0) THEN
C If first ray then calculate flux scaling
		SCALE=0.0
		DO IY=1,IDFM(3,IDEF)
			DO IX=1,IDFM(2,IDEF)
				CALL SRT_GETPIX(IX,IY,1,IDFM(1,IDEF),
     +				IDFM(2,IDEF),IDFM(3,IDEF),
     +				%val(IDFP(1,IDEF)),%val(IDFP(2,IDEF)),
     +				%val(IDFP(3,IDEF)),X,Y,Z,ISTAT)
				SCALE=SCALE+Z
			ENDDO
		ENDDO
		SCALE=DBLE(NRAYS)/SCALE
C Set initial indices etc.
		NHERE=0
		IHERE=0
		IAT=-1	
		NPIX=IDFM(2,IDEF)*IDFM(3,IDEF)
	ENDIF
C Check to see if need to move to next pixel
	IF(IHERE.EQ.NHERE) THEN
		IHERE=0
		NHERE=0	
		DO WHILE(NHERE.EQ.0)
			IAT=IAT+1
			IF(IAT.EQ.NPIX) IAT=0
			IX=MOD(IAT,IDFM(2,IDEF))+1
			IY=IAT/IDFM(2,IDEF)+1
			CALL SRT_GETPIX(IX,IY,1,IDFM(1,IDEF),
     +			IDFM(2,IDEF),IDFM(3,IDEF),
     +			%val(IDFP(1,IDEF)),%val(IDFP(2,IDEF)),
     +			%val(IDFP(3,IDEF)),X,Y,Z,ISTAT)
			NHERE=NINT(Z*SCALE)
		ENDDO
	ENDIF
C Set current offsets
	DX=X
	DY=Y
C Increment number of rays done at current pixel
	IHERE=IHERE+1
	NDONE=NDONE+1
C If last ray or finished last pixel reset 
	IF(NDONE.EQ.NRAYS) THEN
		NDONE=0
	ENDIF
	END
