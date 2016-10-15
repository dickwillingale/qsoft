*+SRT_CNIN      Intersection of a ray with a conic at normal incidence
	SUBROUTINE SRT_CNIN(DIR,XA,CAX,CAR,CVR,PC,IDEF,IKON,HIT,PH,CNL,
     +	ISTAT)
	IMPLICIT NONE
	DOUBLE PRECISION DIR(3),XA(3),CAX(3),CAR(3),CVR(3),PC(7)
	DOUBLE PRECISION PH(3),CNL(3)
	INTEGER IDEF(2),IKON,ISTAT
	LOGICAL HIT
*DIR    input   direction of ray
*XA     input   origin of ray
*CAX    input   axis of cone
*CAR    input   azimuthal reference direction
*CVR    input   vertex of cone
*PC     input   quadratic coefficients of conic and limits
*IDEF   input   deformation indices
*IKON   input   configuration
*HIT    output  .TRUE. if hit
*PH     output  intersection point
*CNL    output  normal to surface at intersection point
*ISTAT  in/out  returned status
* The configuration is given by:
* IKON          incidence       deformation     limits
*  1            normal          axial           radial
*  2            normal          axial           cartesian
*  3		normal		axial		azimuthal
* The conic equation is of the form r**2=PC(1).t**2+PC(2).t+PC(3) where
* r is the radius from the axis and t is the distance along axis from the
* vertex.
*-Author Dick Willingale 1996-Dec-2
	INCLUDE 'SRT_COM'
	DOUBLE PRECISION T1,COSA,P(3),D1,D2,DM,A,B,C,VP(3),VY(3),VR(3)
	DOUBLE PRECISION T,Q,SRT_QUFD,DTDX,DTDY,DN
	DOUBLE PRECISION R,X,Y,DT,TDR,WS,PS,HW
	INTEGER IH,NPASS,J
	LOGICAL SEARCH
	DOUBLE PRECISION SMALL
	PARAMETER (SMALL=1.D-3)
C Calculate local y axis
	CALL SRT_VCRS(CAX,CAR,VY)
C Find coefficients of radial distance of ray from axis
C r**2=P(1).x**2+P(2).x+P(3) where t=x.COSA+T1, 
C x is distance along ray from XA
	CALL SRT_RCON(DIR,XA,CAX,CVR,T1,COSA,P)
C Start with no deformation
	DT=0.0
C Search for intersection
	DM=0.0
	SEARCH=.TRUE.
	NPASS=0
	DO WHILE(SEARCH)
C Calculate coefficients of quadratic of separation of ray and conic
C Include axial deformation
		A=P(1)-(PC(1)*COSA**2)
		B=P(2)-(PC(2)*COSA+PC(1)*T1*COSA*2.0)
		C=P(3)-(PC(3)+PC(1)*T1**2+PC(2)*T1)
		B=B-(2.0*PC(1)*COSA*DT)
		C=C-(2.0*PC(1)*T1*DT+PC(1)*DT**2+PC(2)*DT)
C Solve quadratic for distance of intersection points along ray
		CALL SRT_QURT(A,B,C,IH,D1,D2)
C Look for point at minimum positive distance along ray
		IF(IH.EQ.1.AND.D1.GT.0.0) THEN
			HIT=.TRUE.
			DN=D1
		ELSEIF(IH.EQ.2) THEN
			DN=MIN(D1,D2)
			IF(DN.GT.0.0) THEN
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
C Calculate axial position of intersection point
			T=T1+DN*COSA
C Find intersection point, radial vector and radius
			DO J=1,3
				PH(J)=XA(J)+DIR(J)*DN
				VP(J)=CVR(J)+CAX(J)*T                        
			ENDDO
			CALL SRT_DIDI(VP,PH,VR,R)
C Find local x, y and azimuth
			CALL SRT_VDOT(VR,CAR,X)
			X=X*R
			CALL SRT_VDOT(VR,VY,Y)
			Y=Y*R
			Q=ATAN2(Y,X)
			IF(IDEF(1).GT.0) THEN
C Axial deformation
				CALL SRT_DFRM(IDEF,X,Y,DT,DTDX,DTDY,ISTAT)
				IF(ABS(DN-DM).GT.SMALL) THEN
					DM=DN
				ELSE
					SEARCH=.FALSE.
				ENDIF
			ELSE
				SEARCH=.FALSE.
			ENDIF
		ELSE
			SEARCH=.FALSE.
		ENDIF
		NPASS=NPASS+1
		IF(NPASS.GT.10) THEN
			SEARCH=.FALSE.  
		ENDIF
	ENDDO
	IF(HIT) THEN
C Check limits
		IF(IKON.EQ.1) THEN
			IF(R.LT.PC(5).OR.R.GE.PC(7)) THEN
				HIT=.FALSE.
			ENDIF
		ELSEIF(IKON.EQ.2) THEN
			IF(X.LT.PC(4).OR.X.GE.PC(6).OR.
     +                  Y.LT.PC(5).OR.Y.GE.PC(7)) THEN
				HIT=.FALSE.
			ENDIF
		ELSEIF(IKON.EQ.3) THEN
			IF(Q.LT.0.0) THEN	
				Q=2.0*PI+Q
			ENDIF
			WS=2.0*PI/PC(4)
			PS=MOD(Q,WS)*R
			HW=(PC(5)+PC(6)*R)*0.5
			IF(PS.GT.HW.AND.PS.LE.WS*R-HW) THEN
				HIT=.FALSE.
			ENDIF
		ENDIF
	ENDIF
	IF(HIT) THEN
		TDR=SRT_QUFD(PC,T+DT)
		IF(TDR.NE.0.0) THEN
			TDR=2.0*R/TDR
C Calculate normal to surface including deformation
			DO J=1,3
		  		CNL(J)=-CAX(J)+VR(J)*TDR+CAR(J)*DTDX+
     +                  	VY(J)*DTDY
			ENDDO
			CALL SRT_VNRM(CNL,ISTAT)
		ELSE
			CNL(2)=VR(1)
			CNL(2)=VR(2)
			CNL(3)=VR(3)
		ENDIF
		IF(NPASS.GT.10) THEN
			WRITE(*,*) 'SRT_CNIN warning - failed to converge'
		ENDIF
	ENDIF
	END
