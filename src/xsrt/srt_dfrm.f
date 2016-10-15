*-SRT_DFRM      Deformation of surfaces
	SUBROUTINE SRT_DFRM(IDEF,X,Y,Z,DZDX,DZDY,ISTAT)
	IMPLICIT NONE
	INTEGER IDEF(2),ISTAT
	DOUBLE PRECISION X,Y,Z,DZDX,DZDY
*IDEF   input   IDEF(1)=main deformation index, IDEF(2)=internal index
*X,Y    input   position on surface (in surface coordinates)
*Z      output  deformation normal to surface 
*DZDX   output  gradient of deformation wrt X
*DZDY   output  gradient of deformation wrt Y
*ISTAT  in/out  returned status
*Note for a surface of revolution X is axial direction and Y is azimuth.
*-Author Dick Willingale 1996-Nov-21
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
        CALL SRT_IN2D(X,Y,IDEF(2),IDFM(1,ID),IDFM(2,ID),IDFM(3,ID),
     +  DSAM(IDFP(1,ID)),DSAM(IDFP(2,ID)),DSAM(IDFP(3,ID)),
     +  DSAM(IDFP(4,ID)),DSAM(IDFP(5,ID)),Z,DZDX,DZDY,ISTAT)
	END
