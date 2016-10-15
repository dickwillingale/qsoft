*+SRT_TDEF	Tabulate deformation
	SUBROUTINE SRT_TDEF(IDEF,NX,NY,X,Y,Z,DZDX,DZDY,ISTAT)
	IMPLICIT NONE
	INTEGER IDEF(2),NX,NY,ISTAT
	DOUBLE PRECISION X(NX),Y(NY),Z(NX,NY),DZDX(NX,NY),DZDY(NX,NY)
*IDEF	input	deformation index and submatrix
*NX	input	size in x direction
*NY	input	size in y direction
*X	input	array of x values
*Y	input	array of y values
*Z	output	array of z deformation samples
*DZDX	output	array of gradient samples
*DZDY	output	array of gradient samples
*ISTAT	in/out	returned status
*-Author Dick Willingale 1997-Jan-16
	INTEGER J,K
C
	IF(ISTAT.NE.0) RETURN
C	
	DO J=1,NY
		DO K=1,NX
			CALL SRT_DFRM(IDEF,X(K),Y(J),Z(K,J),DZDX(K,J),
     +			DZDY(K,J),ISTAT)
		ENDDO
	ENDDO
	END
