*+SRT_PLNA         Intersection of a ray with a plane aperture
	SUBROUTINE SRT_PLNA(DIR,XA,PPT,PLN,XRF,P,IDEF,IKON,IH,ANS,RN,
     +	ISTAT)
	IMPLICIT NONE
	DOUBLE PRECISION DIR(3),XA(3),PPT(3),PLN(3),XRF(3),P(9)
	DOUBLE PRECISION ANS(3),RN(3)
	INTEGER IDEF(2),IKON,ISTAT,IH
*DIR    input   direction of ray
*XA     input   origin of ray
*PPT    input   point on plane
*PLN    input   normal to plane
*XRF	input	reference axis in plane
*P	input	limits
*IDEF	input	deformation indices
*IKON	input 	configuration
*IH	output  aperture number hit (0 if not hit)
*ANS    output  intersection point on plane
*RN	output	normal to plane at intersection point (including deformation)
* The configuration is given by:
* IKON		deformation	limits
*  1		normal		cartesian (hole)
*  2		normal		radial
*  3		radial		radial
*  -ve		radial		nested radial
*  4		normal		azimuthal
*  5		normal		cartesian (block)
*  6		normal		cartesian grid
*  7		normal		cartesian slats
*  8		normal		radial/azimuthal sector
*  9		normal		radial+pores
* 10		normal		parallelogram
* If limits parallelogram P(1)=xmin, P(2)=ymin, P(3)=xmax, P(4)=ymax, P(5)=dx
* If limits cartesian then P(1)=xmin, P(2)=ymin, P(3)=xmax, P(4)=ymax
* If limits cartesian slats then P(1)=xmin, P(2)=ymin, P(3)=xmax, P(4)=ymax
* P(5)=pitchx, P(6)=ribx
* If limits cartesian grid P(1)=pitchx, P(2)=pitchy, P(3)=ribx, P(4)=riby
* If limits radial then P(1)=aref, P(2)=rmin, P(3)=rmax, P(4)=rmin, P(5)=rmax..
* If limits azimuthal then P(1)=nsec, P(2)=cw, P(3)=aw where arm
* width= cw+aw*r
* If limits radial/azimuthal P(1)=rmin, P(2)=rmax, P(3)=amin, P(4)=amax
* If deformation is normal then this is used during search
* If deformation is radial then we perturb the radial limits
* If IKON -ve then ABS(IKON) is the number of nested apertures. In this
* case search for aperture which contains the ray. The aperture radii must
* start with the largest (outer) first. 2 calls are required to SRT_DFRM,
*	Use IDEF(ISHELL) for outer radius
*	Use IDEF(ISHELL)+1 for inner radius
*ISTAT	in/out	returned status
*-Author Dick Willingale 1996-Nov-28
	INCLUDE 'SRT_COM'
	DOUBLE PRECISION DR(3),S,A,G,VP(3),R,X,Y,VY(3),DN,DM
	DOUBLE PRECISION DHDX,DHDY,DH,DL,AZ,WS,PS,HW,XX,YY,GR
	INTEGER J,JJ,IDL(3),NPASS
	LOGICAL SEARCH
	DOUBLE PRECISION SMALL
	PARAMETER (SMALL=1.E-4)
C
	IF(ISTAT.NE.0.0) RETURN
C Calculate cos of angle between ray and normal
	CALL SRT_VDOT(PLN,DIR,A)
	IF(A.EQ.0.0) THEN
		IH=0
	ELSE
C Calculate direction and distance to point on plane
		CALL SRT_DIDI(XA,PPT,DR,S)
C Calculate cos of angle between line to point on plane and normal
		CALL SRT_VDOT(PLN,DR,G)
C Find y reference axis in plane
		CALL SRT_VCRS(PLN,XRF,VY)
C Start with no deformation in the normal direction
		DHDX=0.0
		DHDY=0.0
		DH=0.0
C Search for intersection
		SEARCH=.TRUE.
		NPASS=0
		DM=0.0
		DO WHILE(SEARCH)
C Calculate intersection position
			DN=S*G+DH
			DO J=1,3
				ANS(J)=DIR(J)*DN/A+XA(J)
			ENDDO
C Find vector from origin on plane and intersection
			CALL SRT_DIDI(PPT,ANS,VP,R)
C Calculate local x and y on plane
			CALL SRT_VDOT(VP,XRF,X)
			X=X*R
			CALL SRT_VDOT(VP,VY,Y)
			Y=Y*R
			IF((IKON.EQ.1.OR.IKON.EQ.2.OR.IKON.EQ.4.OR.IKON.EQ.5
     +			.OR.IKON.EQ.6.or.ikon.eq.7.or.ikon.eq.8)
     +			.AND.IDEF(1).GT.0
     +			.AND.ABS(DN-DM).GT.SMALL) THEN
				DM=DN
C Get deformation
				CALL SRT_DFRM(IDEF,X,Y,DH,DHDX,DHDY,ISTAT)
			ELSE
C Converged
				SEARCH=.FALSE.
			ENDIF
			NPASS=NPASS+1
			IF(NPASS.GT.10) THEN
			  WRITE(*,*)  'SRT_PLNA error - failed to converge'
			  SEARCH=.FALSE.
			ENDIF
		ENDDO
C Check limits
		IF(IKON.EQ.1) THEN
			IH=1
			IF(X.LT.P(1).OR.X.GT.P(3)
     +			.OR.Y.LT.P(2).OR.Y.GT.P(4)) THEN
				IH=0
			ENDIF
		ELSEIF(IKON.EQ.2) THEN
			IH=1
			IF(R.LT.P(1).OR.R.GT.P(2)) THEN
				IH=0
			ENDIF
		ELSEIF(IKON.EQ.3) THEN
			IH=0
			IF(IDEF(1).GT.0) THEN
				IDL(1)=IDEF(1)
				IDL(2)=IDEF(2)
				AZ=ATAN2(Y,X)
				CALL SRT_DFRM(IDL,P(1),AZ,DH,DHDX,DHDY,ISTAT)
C Deformation of low radius is obtained using higher index
				IDL(2)=IDEF(2)+1
				CALL SRT_DFRM(IDL,P(1),AZ,DL,DHDX,DHDY,ISTAT)
			ELSE
				DL=0.0
				DH=0.0
			ENDIF
			IF(R.GT.P(2)+DL.AND.R.LE.P(3)+DH) THEN
				IH=1
			ENDIF
			DHDX=0.0
			DHDY=0.0
		ELSEIF(IKON.EQ.4) THEN
			IH=0
			AZ=ATAN2(Y,X)
			IF(AZ.LT.0.0) THEN	
				AZ=2.0*PI+AZ
			ENDIF
			WS=2.0*PI/P(1)
			PS=MOD(AZ,WS)*R
			HW=(P(2)+P(3)*R)*0.5
			IF(PS.GT.HW.AND.PS.LE.WS*R-HW) THEN
				IH=1
			ENDIF
		ELSEIF(IKON.EQ.5) THEN
			IH=0
			IF(X.LT.P(1).OR.X.GT.P(3)
     +			.OR.Y.LT.P(2).OR.Y.GT.P(4)) THEN
				IH=1
			ENDIF
		ELSEIF(IKON.EQ.6) THEN
			IH=0
			XX=MOD(X+(P(3)+P(1))*0.5,P(1))
			IF(XX.LT.0.0) XX=P(1)+XX
			YY=MOD(Y+(P(4)+P(2))*0.5,P(2))
			IF(YY.LT.0.0) YY=P(2)+YY
			IF(XX.GT.P(3).AND.YY.GT.P(4)) THEN
				IH=1
			ENDIF
		ELSEIF(IKON.EQ.7) THEN
			IH=1
			IF(X.LT.P(1).OR.X.GT.P(3)
     +			.OR.Y.LT.P(2).OR.Y.GT.P(4)) THEN
				IH=0
			ELSE
				XX=MOD(X-P(1),P(5))
				IF(XX.LT.P(6)) THEN
					IH=0
				ENDIF
			ENDIF
		ELSEIF(IKON.EQ.8) THEN
			IH=0
			AZ=ATAN2(Y,X)
			IF(AZ.LT.0.0) THEN	
				AZ=2.0*PI+AZ
			ENDIF
			IF(R.GT.P(1).AND.R.LE.P(2).AND.AZ.GT.P(3)
     +			.AND.AZ.LT.P(4)) THEN
				IH=1
			ENDIF
			DHDX=0.0
			DHDY=0.0
		ELSEIF(IKON.EQ.9) THEN
			IH=0
			IF(R.LT.P(1).OR.R.GT.P(2)) THEN
				IH=1
			ENDIF
			P(8)=X
			P(9)=Y
		ELSEIF(IKON.EQ.10) THEN
			IH=1
			GR=Y*P(5)/(P(4)-P(2))
			IF(X.LT.(P(1)+GR).OR.X.GT.(P(3)+GR)
     +			.OR.Y.LT.P(2).OR.Y.GT.P(4)) THEN
				IH=0
			ENDIF
		ELSEIF(IKON.LT.0) THEN
C Nested set of aperture
			IH=0
			IDL(1)=IDEF(1)
			AZ=ATAN2(Y,X)
			J=0
			DO WHILE(IH.EQ.0.AND.J.LT.-IKON)
				J=J+1
				JJ=J+J
				IF(IDEF(1).GT.0) THEN
					IDL(2)=J
					CALL SRT_DFRM(IDL,P(1),AZ,
     +					DH,DHDX,DHDY,ISTAT)
C Deformation of low radius is obtained using higher index
					IDL(2)=J+1
					CALL SRT_DFRM(IDL,P(1),AZ,
     +					DL,DHDX,DHDY,ISTAT)
				ELSE
					DL=0.0
					DH=0.0
				ENDIF
				IF(R.GT.P(JJ)+DL.AND.R.LE.P(JJ+1)+DH) THEN
					IH=J
				ENDIF
			ENDDO	
			DHDX=0.0
			DHDY=0.0
		ENDIF
		IF(IH.GT.0) THEN
C Calculate normal including deformation
			DO J=1,3
				RN(J)=PLN(J)-XRF(J)*DHDX-VY(J)*DHDY
			ENDDO
			CALL SRT_VNRM(RN,ISTAT)
		ENDIF
	ENDIF
	END
