*+SRT_CNIC      Intersection of a ray with a conic grazing incidence
	SUBROUTINE SRT_CNIC(DIR,XA,CAX,CAR,CVR,PC,IDEF,IKON,HIT,PH,CNL,
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
*  1            grazing         radial          axial
*  2            grazing         radial          axial+azimuthal
*  3            normal          radial          axial+circumferencial
* The conic equation is of the form r**2=PC(1).t**2+PC(2).t+PC(3) where
* r is the radius from the axis and t is the distance along axis from the
* vertex.
*-Author Dick Willingale 1996-Nov-16
	DOUBLE PRECISION T1,COSA,P(3),D1,D2,A,B,C,VP(3),VY(3),VR(3),VA(3)
	DOUBLE PRECISION T,Q,SRT_QUFU,SRT_QUFD,DR,DRDT,DRDQ,DN
	DOUBLE PRECISION R,X,Y,RDT,RM,RRDD
	INTEGER IHIT,NPASS,J
	LOGICAL SEARCH
	DOUBLE PRECISION SMALL
	PARAMETER (SMALL=1.D-4)
C Calculate local y axis
	CALL SRT_VCRS(CAX,CAR,VY)
C Find coefficients of radial distance of ray from axis
C r**2=P(1).x**2+P(2).x+P(3) where t=x.COSA+T1, 
C x is distance along ray from XA
	CALL SRT_RCON(DIR,XA,CAX,CVR,T1,COSA,P)
C Calculate coefficients of quadratic of separation of ray and conic
	A=P(1)-(PC(1)*COSA**2)
	B=P(2)-(PC(2)*COSA+PC(1)*T1*COSA*2.0)
	C=P(3)-(PC(3)+PC(1)*T1**2+PC(2)*T1)
C Solve quadratic for distance of intersection points along ray
	CALL SRT_QURT(A,B,C,IHIT,D1,D2)
C Look for point at minimum positive distance along ray
	IF(IHIT.EQ.1.AND.D1.GT.SMALL) THEN
		HIT=.TRUE.
		DN=D1
	ELSEIF(IHIT.EQ.2) THEN
		IF(D1.LE.D2) THEN
			T=T1+D1*COSA
			IF(D1.GT.SMALL.AND.T.GE.PC(4).AND.T.LT.PC(6)) THEN
				HIT=.TRUE.
				DN=D1
			ELSEIF(D2.GT.SMALL) THEN
				HIT=.TRUE.
				DN=D2
			ELSE
				HIT=.FALSE.
			ENDIF
		ELSE
			T=T1+D2*COSA
			IF(D2.GT.SMALL.AND.T.GE.PC(4).AND.T.LT.PC(6)) THEN
				HIT=.TRUE.
				DN=D2
			ELSEIF(D1.GT.SMALL) THEN
				HIT=.TRUE.
				DN=D1
			ELSE
				HIT=.FALSE.
			ENDIF
		ENDIF
	ELSE
		HIT=.FALSE.
	ENDIF
	IF(HIT) THEN
C Search for intersection including deformation
		SEARCH=.TRUE.
		NPASS=0
		DO WHILE(SEARCH)
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
C Find radius and gradient of conic generating function
			RM=SQRT(SRT_QUFU(PC,T))
			RDT=SRT_QUFD(PC,T)/(2.0*RM)
			DRDQ=0.0
			IF(IDEF(1).GT.0) THEN
C Include radial deformation
				CALL SRT_DFRM(IDEF,T,Q,DR,DRDT,DRDQ,ISTAT)
				RM=RM+DR
				RDT=RDT+DRDT
C Find gradient of separation wrt distance alone ray
				RRDD=SRT_QUFD(P,DN)/(2.0*R)
				IF(COSA.NE.0.0) THEN
					RRDD=RRDD-RDT/COSA
				ELSE
					RRDD=1.0
				ENDIF
				IF(ABS(R-RM).GT.SMALL.AND.RRDD.NE.0.0) THEN
C Calculate search step along ray by Newton-Raphson
					DN=DN-(R-RM)/RRDD
				ELSE
					SEARCH=.FALSE.
				ENDIF
			ELSE
C If no deformation then found intersection
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
				IF(T.LT.PC(4).OR.T.GE.PC(6)) THEN
					HIT=.FALSE.
				ENDIF
			ELSEIF(IKON.EQ.2) THEN
				IF(T.LT.PC(4).OR.T.GE.PC(6)) THEN
                                	HIT=.FALSE.
                        	ENDIF
				IF(Q.LT.PC(5).OR.Q.GE.PC(7)) THEN
                                	HIT=.FALSE.
                        	ENDIF
			ELSEIF(IKON.EQ.3) THEN
				IF(T.LT.PC(4).OR.T.GE.PC(6)) THEN
                                	HIT=.FALSE.
                        	ENDIF
				IF(X.LT.PC(5).OR.X.GE.PC(7)) THEN
                                	HIT=.FALSE.
                        	ENDIF
			ENDIF
		ENDIF
		IF(HIT) THEN	
		IF(R.GT.0.0) THEN
C Calculate normal to surface including deformation
			CALL SRT_VCRS(CAX,VR,VA)
			DO J=1,3
			  CNL(J)= VR(J)-CAX(J)*RDT-VA(J)*DRDQ/R
			ENDDO
			CALL SRT_VNRM(CNL,ISTAT)
		ELSE
			CNL(1)=CAX(1)
			CNL(2)=CAX(2)
			CNL(3)=CAX(3)
			IF(NPASS.GT.10) THEN
			    WRITE(*,*) 'SRT_CNIC warning',
     +		            ' - failed to converge, last delta r ',R-RM
			ENDIF
		ENDIF
		ENDIF
	ENDIF
	END
