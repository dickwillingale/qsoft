*+SRT_SPHR      Intersection of a ray with a sphere
	SUBROUTINE SRT_SPHR(DIR,XA,CAX,CAR,CVR,PC,IDEF,IKON,HIT,PH,CNL,
     +	ISTAT)
	IMPLICIT NONE
	DOUBLE PRECISION DIR(3),XA(3),CAX(3),CAR(3),CVR(3),PC(9)
	DOUBLE PRECISION PH(3),CNL(3)
	INTEGER IDEF(2),IKON,ISTAT
	LOGICAL HIT
*DIR    input   direction of ray
*XA     input   origin of ray
*CAX    input   axis of sphere, from centre to local origin on surface
*CAR    input   reference direction
*CVR    input   centre of sphere
*PC     input   radius of sphere and limits
*IDEF   input   deformation indices
*IKON   input   configuration
*HIT    output  .TRUE. if hit
*PH     output  intersection point
*CNL    output  normal to surface at intersection point
*ISTAT  in/out  returned status
* The configuration is given by:
* IKON          incidence       deformation     limits
*  1            grazing         radial          cartesian
*  2		grazing		radial		radial
*  3            normal   	radial          cartesian
*  4            normal          radial          radial
*  5		normal		radial		radial+grid
*  6		normal		radial		cartesian+grid
* Local x on the surface is given by CAR
* Local y on the surface is perpendicular to CAX and CAR
* x and y are positions in the tangential plane
* PC(1) is the radius of the sphere
* If cartesian limits PC(2)=xmin, PC(3)=ymin, PC(4)=xmax, PC(5)=ymax
* If limits radial or cartesian +grid then
*			PC(2)=aperture radius or aperture half width
*			PC(3)=pitch of grid in x
*			PC(4)=pitch of grid in y
*			PC(5)=rib width of grid in x
*			PC(6)=rib width of grid in y
*			PC(7)=pore length
*			PC(8)=intersection x (returned)
*			PC(9)=intersection y (returned)
*			PC(10)=surface index
* radial/cartesian grid added 2004-Feb for MCP ray tracing RW
*-Author Dick Willingale 1996-Nov-28
	DOUBLE PRECISION P(3),D1,D2,DM,A,B,C,VP(3),VR(3)
	DOUBLE PRECISION VY(3),X,Y,R,RXY,XX,YY
	DOUBLE PRECISION DR,DRDX,DRDY,DN
	DOUBLE PRECISION OPOS(3),PH1,PH2
	INTEGER IHIT,NPASS,J
	LOGICAL SEARCH
	DOUBLE PRECISION PI,SMALL
	PARAMETER (PI=3.14159265359,SMALL=1.D-6)
C	
	IF(ISTAT.NE.0) RETURN
C Find quadratic for distance of ray from centre. r**2=P(1).d**2+p(2).d+P(3)
C where d is the distance along the ray from the ray origin
	CALL SRT_VDOT(DIR,DIR,P(1))
	CALL SRT_DIDI(XA,CVR,VP,P(3))
	CALL SRT_VDOT(VP,DIR,P(2))
	P(2)=-P(2)*2.0*P(3)
	P(3)=P(3)**2
C Start with no deformation
	DRDX=0.0
	DR=0.0
	DRDY=0.0
	X=0.0
	Y=0.0
	R=0.0
C Search for intersection
	DM=0.0
	SEARCH=.TRUE.
	NPASS=0
	DO WHILE(SEARCH)
C Calculate coefficients of quadratic of separation of ray and sphere
		A=P(1)
		B=P(2)
		C=P(3)-PC(1)**2
C Include radial deformation, ignore deformation gradients
		C=C-(DR**2)
C Solve quadratic for distance of intersection points along ray
		CALL SRT_QURT(A,B,C,IHIT,D1,D2)
C Find origin on surface
		DO J=1,3
			OPOS(J)=CVR(J)+CAX(J)
		ENDDO
C Look for points at positive distance along ray
		IF(IHIT.EQ.1.AND.D1.GT.0.0) THEN
			HIT=.TRUE.
			DN=D1
		ELSEIF(IHIT.EQ.2) THEN
			IF(D1.GT.0.0.AND.D2.GT.0.0) THEN
C Both intersection points in +ve ray direction
C Find distance of both intersection points from surface origin
				PH1=0.0
				PH2=0.0
				DO J=1,3
					PH1=PH1+(XA(J)+DIR(J)*D1-OPOS(J))**2
					PH2=PH2+(XA(J)+DIR(J)*D2-OPOS(J))**2
				ENDDO
				IF(PH1.LT.PH2) THEN
					DN=D1
				ELSE
					DN=D2
				ENDIF
				HIT=.TRUE.
			ELSE
				DN=MAX(D1,D2)
				IF(DN.GT.0.0) THEN
					HIT=.TRUE.
				ELSE
					HIT=.FALSE.
				ENDIF
			ENDIF
		ELSE
			HIT=.FALSE.
		ENDIF
		IF(HIT) THEN
C Find intersection point
			DO J=1,3
				PH(J)=XA(J)+DIR(J)*DN
			ENDDO
C Find vector from centre to intersection and radial distance
			CALL SRT_DIDI(CVR,PH,VR,R)
C Calculate local x of intersection point
			CALL SRT_VDOT(VR,CAR,X)
			X=X*R
C Calculate local y of intersection point
			CALL SRT_VCRS(CAX,CAR,VY)
			CALL SRT_VDOT(VR,VY,Y)
			Y=Y*R
			IF(IDEF(1).GT.0.AND.ABS(DN-DM).GT.SMALL) THEN
C Still not converged 
				DM=DN
C Perturb intersection point by deformation
				CALL SRT_DFRM(IDEF,X,Y,DR,DRDX,DRDY,ISTAT)
			ELSE
C Converged
				SEARCH=.FALSE.
			ENDIF
		ELSE
			SEARCH=.FALSE.
		ENDIF
		NPASS=NPASS+1
		IF(NPASS.GT.10) THEN
			WRITE(*,*) 'SRT_SPHR error - failed to converge'
			SEARCH=.FALSE.  
		ENDIF
	ENDDO
	IF(HIT) THEN
C Check limits
                IF(IKON.EQ.5) THEN
C Return local impact coordinates and check aperture radius
                        PC(8)=X
			PC(9)=Y
			RXY=SQRT(X**2+Y**2)
        		IF(RXY.LT.PC(2)) THEN
                		HIT=.FALSE.
			ENDIF
                ELSEIF(IKON.EQ.6) THEN
C Return local impact coordinates and check aperture half length
                        PC(8)=X
			PC(9)=Y
        		IF(ABS(X).LT.PC(2).AND.ABS(Y).LT.PC(2)) THEN
                		HIT=.FALSE.
			ENDIF
		ELSEIF(IKON.EQ.2.OR.IKON.EQ.4) THEN
			RXY=SQRT(X**2+Y**2)
			IF(RXY.LT.PC(2).OR.RXY.GE.PC(3)) THEN
				HIT=.FALSE.
			ENDIF
		ELSE
			IF(X.LT.PC(2).OR.X.GE.PC(4).OR.
     +                  Y.LT.PC(3).OR.Y.GE.PC(5)) THEN
				HIT=.FALSE.
			ENDIF
		ENDIF
C Calculate normal to surface including deformation
		DO J=1,3
			CNL(J)= VR(J)-CAR(J)*DRDX-VY(J)*DRDY
		ENDDO
		CALL SRT_VNRM(CNL,ISTAT)
	ENDIF
	END
