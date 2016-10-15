*+SRT_MATS	Set deformation arrays
	SUBROUTINE SRT_MATS(NSM,NS,NX,NY,X,Y,Z,NN,G1,G2,XD,YD,ZD,
     +	GXD,GYD,ISTAT)
	IMPLICIT NONE
	INTEGER NSM,NS,NX,NY,ISTAT,NN
	DOUBLE PRECISION X(NX),Y(NY),Z(NX,NY),G1(NN),G2(NN)
	DOUBLE PRECISION XD(NX,NSM),YD(NY,NSM),ZD(NX,NY,NSM)
	DOUBLE PRECISION GXD(NX,NY,NSM),GYD(NX,NY,NSM)
*NSM	input	total number of sub-matrices
*NS	input	index of sub-matrix
*NX	input	number of x samples
*NY	input	number of y samples
*X	input	x values
*Y	input	y values
*Z	input	z (deformation) values
*NN	input	size of gradient workspace arrays
*G1	input	work area for gradient values
*G1	input	work area for gradient values
*XD	in/out	x values
*YD	in/out	y values
*ZD	in/out	z values
*GXD	in/out	x gradient
*GYD	in/out	y gradient
*ISTAT	in/out	returned status
*-Author Dick Willingale 1996-Dec-11
	INTEGER J,K
C
	IF(ISTAT.NE.0) RETURN
C
	DO J=1,NX
		XD(J,NS)=X(J)
		IF(X(J).LT.-1.D5.OR.X(J).GT.1.D5) THEN
			WRITE(*,*) 'srt_mats error - dodgy x',X(J),' in ',NS
			ISTAT=1
		ENDIF
	ENDDO
	DO J=1,NY
		YD(J,NS)=Y(J)
		IF(Y(J).LT.-1.D5.OR.Y(J).GT.1.D5) THEN
			WRITE(*,*) 'srt_mats error - dodgy y',Y(J),' in ',NS
			ISTAT=1
		ENDIF
	ENDDO
	DO J=1,NY
		DO K=1,NX
			G1(K)=Z(K,J)
			ZD(K,J,NS)=Z(K,J)
			IF(Z(K,J).LT.-1.D5.OR.Z(K,J).GT.1.D5) THEN
			WRITE(*,*) 'srt_mats error - dodgy z',Z(K,J),' in ',NS
			ISTAT=1
			ENDIF
		ENDDO
		CALL SRT_GRAD(NX,X,G1,0,G2)
		DO K=1,NX
			GXD(K,J,NS)=G2(K)
		ENDDO
	ENDDO
	DO K=1,NX
		DO J=1,NY
			G1(J)=Z(K,J)
		ENDDO
		CALL SRT_GRAD(NY,Y,G1,0,G2)
		DO J=1,NY
			GYD(K,J,NS)=G2(J)
		ENDDO
	ENDDO
	END
