*-SRT_IDFRM      Deformation of surfaces
	SUBROUTINE SRT_IDFRM(IDEF,IX,IY,Z,DZDX,DZDY,ISTAT)
	IMPLICIT NONE
	INTEGER IX,IY,IDEF(2),ISTAT
	DOUBLE PRECISION Z,DZDX,DZDY
*IDEF   input   IDEF(1)=main deformation index, IDEF(2)=internal index
*IX	input	index of position on surface
*IY	input	index of position on surface
*Z      output  deformation normal to surface 
*DZDX   output  gradient of deformation wrt X
*DZDY   output  gradient of deformation wrt Y
*ISTAT  in/out  returned status
*Note for a surface of revolution X is axial direction and Y is azimuth.
*-Author Dick Willingale 2012-Jun-18
	INCLUDE 'SRT_COM'
	INTEGER ID
C
	IF(ISTAT.NE.0) RETURN
C Format of deformation data in common IDFM(3,MAXDF) and IDFP(5,MAXDF)
C Deformation parameter values held in array DSAM
C IDFM(1,IDEF(1))       maximum number of sub-matrices
C IDFM(2,IDEF(1))       number of X sample positions
C IDFM(3,IDEF(1))       number of Y sample positions
C IDFP(1,IDEF(1))       index to X samples
C IDFP(2,IDEF(1))       index to Y samples
C IDFP(3,IDEF(1))       index to Z matrices
C IDFP(4,IDEF(1))       index to X gradient matrices
C IDFP(5,IDEF(1))       index to Y gradient matcices
	ID=IDEF(1)
	IF(ID.LT.1.OR.ID.GT.MAXDF) THEN
		WRITE(*,*) 'SRT_DFRM error - index ',ID,' out of range'
		ISTAT=1
		RETURN
	ENDIF
        CALL SRT_IIN2D(IX,IY,IDEF(2),IDFM(1,ID),IDFM(2,ID),IDFM(3,ID),
     +  DSAM(IDFP(3,ID)),DSAM(IDFP(4,ID)),DSAM(IDFP(5,ID)),
     +	Z,DZDX,DZDY,ISTAT)
	END
