*+SRT_SETD	Set deformation parameters
	SUBROUTINE SRT_SETD(ID,NM,NX,NY,IX,IY,IMZ,IMX,IMY,ISTAT)
	IMPLICIT NONE
	INTEGER ID,NM,NX,NY,ISTAT
        INTEGER IX,IY,IMZ,IMX,IMY
*ID	input	deformation index
*NM	input	number of sub-matrices
*NX	input	number of x samples
*NY	input	number of y samples
*IX	input	index to x samples
*IY	input	index to y samples
*IMZ	input	index to z matrix
*IMX	input	index to x gradient matrix
*IMY	input	index to y gradient matrix
*ISTAT	in/out	returned status
*-Author Dick Willingale 1996-Dec-6
	INCLUDE 'SRT_COM'
	IF(ISTAT.NE.0) RETURN
C Check deformation index
	IF(ID.LT.0.OR.ID.GT.MAXDF) THEN
		WRITE(*,*) 'SRT_SETD error - deformation index out of range'
		ISTAT=1
		RETURN
	ENDIF
C Set parameters
	IDFM(1,ID)=NM
	IDFM(2,ID)=NX
	IDFM(3,ID)=NY
	IDFP(1,ID)=IX
	IDFP(2,ID)=IY
	IDFP(3,ID)=IMZ
	IDFP(4,ID)=IMX
	IDFP(5,ID)=IMY
	END
