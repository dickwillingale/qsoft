*+SRT_SIPORE     Dynamically set up a pair of square Si pores for HPO
	SUBROUTINE SRT_SIPORE(DIR,HPOS,CVR,CAX,CAR,PC,IDEF,ISTAT)
	IMPLICIT NONE
	DOUBLE PRECISION DIR(3),HPOS(3),CAX(3),CAR(3),CVR(3),PC(500)
	INTEGER IDEF(2),ISTAT
*DIR	input	direction of ray
*HPOS	input	position of hit on aperture surface 
*CVR    input   centre of aperture plane (above join plane of Wolter I)
*CAX    input   axis of optic
*CAR    input   reference direction in aperture
*PC     input   remaining parameters
*IDEF   input   deformation indices
*ISTAT  in/out  returned status
* Parameters for pores are:
* PC(1) inner radius of aperture
* PC(2) outer radius of aperture
* PC(3) focal length
* PC(4) radial pitch of pores
* PC(5) azimuthal pitch of pores
* PC(6) wall thickness of pores
* PC(7) surface quality
* PC(8) local x of ray intersection
* PC(9) local y of ray intersection
* PC(10) this surface index
* PC(11) axial curvature 0 conical, 1 Wolter, 2 curved-plane, 3 constant
*	  4 spherical principal surface with Wolter curvature
*	  5 spherical principal surface with conical approximation
*	  6 spherical principal surface with Wolter curvature and
*	  aperture grid
*	  7 spherical principal surface with curved-plane curvature
*	  8 spherical principal surface with plane-curved curvature
*	  9 spherical principal surface with -1/2, 5/2 curvature
*	 10 spherical principal surface with 1/2, 3/2 curvature 
*	 11 spherical principal surface with -1, 3 curvature 
*	 12 spherical principal surface with -3/2, 7/2 curvature 
*	 13 spherical principal surface with -2, 4 curvature 
* PC(12) module frame width
* PC(13) axial distance from aperture (CVR above) to join plane
* Modified to use SRT_FINDMOD RW 2012-Jul-10
* Modified to include grid in front of module RW 2013-May-22
* Axial curvature configuration changed from PC(11) to CMOD from the module
* aperture search 29/04/2015 RW
* Added GMOD as grazing angle ratio for module 29/04/2015 RW
* Added axial variation in figure errors 19/07/2017 RW
*-Author Dick Willingale 2008-Oct-2
	DOUBLE PRECISION PI
	PARAMETER (PI=3.14159265359)
	DOUBLE PRECISION VY(3),X,Y,Z,XP(3),YP(3),ZP(3),RXY,POS(3)
	DOUBLE PRECISION RNM(3),RPOS,CVRN(3)
	DOUBLE PRECISION PP(16),HX,HY,HL,THETA,DTHETA,RCURP,RCURH
	DOUBLE PRECISION DX,DY,DZ,DDX,DDY,DDZ,THETAJ,THETAP,THETAH,DV
	DOUBLE PRECISION FL,PL,RP,RAD,RMIN,DR,DA,DC1,DB1,DC2,DB2,WFR
        DOUBLE PRECISION DBS,DBC,SCA
	DOUBLE PRECISION A2J,WALL,DELTA,HGRID,GX,GY,RAP,ROFF
	DOUBLE PRECISION RMOD,PMOD,TMOD,WMOD,HMOD,XX,YY,RAN(4),TT,DC,DB
	DOUBLE PRECISION CMOD,GMOD,YYY
	DOUBLE PRECISION TFOV,GDEP,GBIT
	INTEGER IX,IY,I,KSUR,IQ,IP,IMOD,ICURV
	LOGICAL HIT
C
	IF(ISTAT.NE.0) RETURN
C Get focal length
	FL=PC(3)
C Find module corresponding to impact position in aperture - module index IMOD
	CALL SRT_FINDMOD(PC(8),PC(9),IMOD,RMOD,PMOD,TMOD,WMOD,HMOD,
     +	PL,CMOD,GMOD,XX,YY)
C module frame width
	WFR=PC(12)
C tracing aperture to join distance
	A2J=PC(13)
C wall thickness
	WALL=PC(6)
C radial pitch
	RAP=PC(4)
C Radial and azimuthal width of pore aperture
	DR=PC(4)-PC(6)
	DA=PC(5)-PC(6)
C
	IF(IMOD.GT.0) THEN
C Get axial curvature 
C	ICURV=INT(PC(11))
		ICURV=INT(CMOD)
C Note GMOD is the grazing angle ratio for this module
		IF(ICURV.EQ.4.or.ICURV.eq.5.or.ICURV.eq.6
     +		.or.ICURV.eq.7.or.ICURV.eq.8
     +		.or.ICURV.eq.9.or.ICURV.eq.10
     +		.or.ICURV.eq.11.or.ICURV.eq.12.or.ICURV.eq.13) THEN
C Spherical principal surface
			DELTA=FL-SQRT(FL**2-RMOD**2)
			IF(ICURV.EQ.6) THEN
C HGRID is axial distance of top of grid from entrance to pore
C DR is the radial width of the pore
C DR/2+GBIT is the half width of the grid aperture
C GDEP is the axial depth of the grid bars
C WALL is the radial thickness of the pore walls
C
C for a free standing grid
C				GBIT=0.02
				GBIT=0.07
				GDEP=0.775
				HGRID=PL*WALL/DR*2.0
C for extension of module plates
C				GBIT=WALL*0.5-0.025
C				HGRID=PL*GBIT/DR
C				GDEP=HGRID-0.025
			ELSE
				HGRID=0.0
				GBIT=0.0
				GDEP=0.0
			ENDIF
		ELSE
C Planar principal surface
			DELTA=0.0
			HGRID=0.0
		ENDIF
C Move along normal to just above grid aperture at top of module.
		DO I=1,3
			CVRN(I)=CVR(I)-CAX(I)*(A2J+DELTA-HGRID-PL-WALL)
		ENDDO
C Get hit position at entrance to module (if grid present then entrance
C to grid)
		CALL SRT_PLNA(DIR,HPOS,CVRN,CAX,CAR,PC,0,9,HIT,POS,RNM,ISTAT)
C Find module again to get new impact position
		CALL SRT_FINDMOD(PC(8),PC(9),IMOD,RMOD,PMOD,TMOD,WMOD,HMOD,
     +		PL,CMOD,GMOD,XX,YY)
C Check inside frame
		WMOD=WMOD-2.0*WFR
		HMOD=HMOD-2.0*WFR
		IF((ABS(XX).GT.WMOD*0.5).OR.(ABS(YY).GT.HMOD*0.5)) THEN
			IMOD=0
		ENDIF
	ENDIF
	IF(IMOD.GT.0) THEN
		ICURV=INT(CMOD)
C Set radial offset of pore wrt the centre of the module aperture
		ROFF=XX
C Find new local hit position radius and angle in aperture
		RXY=SQRT(PC(8)**2+PC(9)**2)
		THETA=ATAN2(PC(9),PC(8))
		DTHETA=PC(5)/RXY
C Radial packing of pores - find centre of pore for ray using hit position
		XX=RXY-(RMOD-WMOD*0.5)
		RXY=(INT(XX/PC(4))+0.5)*PC(4)+(RMOD-WMOD*0.5)
C
		IF(THETA-PMOD.GT.PI) THEN
			PMOD=PMOD+PI*2.0D0
		ELSEIF(THETA-PMOD.LT.-PI) THEN
			PMOD=PMOD-PI*2.0D0
		ENDIF
		YYY=THETA-(PMOD-HMOD*0.5/RMOD)
		THETA=(INT(YYY/DTHETA)+0.5)*DTHETA+(PMOD-HMOD*0.5/RMOD)
		X=RXY*COS(THETA)
		Y=RXY*SIN(THETA)
	ELSE
C Move pore so that ray misses
		X=PC(2)*10.0
		Y=PC(2)*10.0
	ENDIF
C Estimate grazing angle at join plane for this pore
	THETAJ=ATAN2(RXY,(FL-DELTA))/2.0/(1.0+GMOD)
C The end of the pore on principal plane is at a slightly smaller
C radius that the entrance to the grid or pore
	X=X*(1.0-ABS(PL*TAN(THETAJ)*COS(THETA)/X))
	Y=Y*(1.0-ABS(PL*TAN(THETAJ)*SIN(THETA)/Y))
	RXY=SQRT(X**2+Y**2)
C Radius of curvature of 1st surface
	RCURP=RXY+DR
C Radius of curvature of 2nd surface
	RCURH=RXY
C Grazing angle at join for centre of pore
	THETAJ=ATAN2(RXY,(FL-DELTA))/2.0/(1.0+GMOD)
C Grazing angle of parabola (1st surface)
	THETAP=ATAN2(RXY+PL*TAN(THETAJ)/2,FL-DELTA)/2.0/(1.0+GMOD)
C Grazing angle of hyperbola (2nd surface)
	THETAH=THETAP*GMOD
C Find Y reference axis at local origin on aperture
	CALL SRT_VCRS(CAX,CAR,VY)
	CALL SRT_VNRM(VY,ISTAT)
C Find position of pore centre on join plane
	POS(1)=CVR(1)-CAX(1)*(A2J+DELTA)+X*CAR(1)+Y*VY(1)
	POS(2)=CVR(2)-CAX(2)*(A2J+DELTA)+X*CAR(2)+Y*VY(2)
	POS(3)=CVR(3)-CAX(3)*(A2J+DELTA)+X*CAR(3)+Y*VY(3)
C Deformations
	IF(IDEF(1).GT.0) THEN
C DX is effective displacement of module focus in direction CAR
C DY is effective displacement of module focus in direction CAX X CAR
C DZ is change in focal length of module
C DB mean value for Gaussian in-plane figure errors radians
C DC mean value for Gaussian out-of-plane figure errors radians
C Note TMOD provides the module rotation about the axis which gives
C a contribution to the out-of-plane error
C Also note that the derivatives are not used
		IDEF(2)=1
		CALL SRT_IDFRM(IDEF,IMOD,1,DX,DDX,DDY,ISTAT)
		IDEF(2)=2
		CALL SRT_IDFRM(IDEF,IMOD,1,DY,DDX,DDY,ISTAT)
		IDEF(2)=3
		CALL SRT_IDFRM(IDEF,IMOD,1,DZ,DDX,DDY,ISTAT)
		IDEF(2)=4
		CALL SRT_IDFRM(IDEF,IMOD,1,DB,DDX,DDY,ISTAT)
		IDEF(2)=5
		CALL SRT_IDFRM(IDEF,IMOD,1,DC,DDX,DDY,ISTAT)
		CALL SYS_GAUSS(4,RAN,0.0D0,1.0D0,ISTAT)
C Scale in-plane and out-of-plane figure errors according to distance
C from centre of plate towards the axial edges of the plate
                IF(HMOD.GT.0.0) THEN
                  SCA=1.0+(YY/HMOD*2.0)**2
		ELSE
		  SCA=1.0
		ENDIF
                DB=DB*SCA
                DC=DC*SCA
		DB1=RAN(1)*DB
		DB2=RAN(2)*DB
		DC1=RAN(3)*DC+TMOD
		DC2=RAN(4)*DC+TMOD
	ELSE
		DX=0.0
		DY=0.0
		DZ=0.0
		DB1=0.0
		DB2=0.0
		DC1=TMOD
		DC2=TMOD
	ENDIF
C Find distance to vertex of 1st surface
	Z=RCURP/TAN(THETAP)
C Convert DZ focal length error to a shift of the vertex in-plane
	DV=DZ*ROFF/(FL+DZ)
C Find direction of pore axis including deformations
	DBC=DB1*COS(THETA)*Z
	DBS=DB1*SIN(THETA)*Z
	ZP(1)=Z*CAX(1)+(X-DX+DBC+DV)*CAR(1)+(Y-DY+DBS)*VY(1)
	ZP(2)=Z*CAX(2)+(X-DX+DBC+DV)*CAR(2)+(Y-DY+DBS)*VY(2)
	ZP(3)=Z*CAX(3)+(X-DX+DBC+DV)*CAR(3)+(Y-DY+DBS)*VY(3)
	CALL SRT_VNRM(ZP,ISTAT)
C Now set up axes tangential to surface at pore centre
C Note that there is no pore centred at the origin so can use the cross product
	CALL SRT_VCRS(CAX,ZP,XP)
	CALL SRT_VNRM(XP,ISTAT)
	CALL SRT_VCRS(ZP,XP,YP)
C Apply out-of-plane figure error and module rotation by rotating pore
	IF(DC1.NE.0.0) THEN
		DO I=1,3
			TT=XP(I)*COS(DC1)+YP(I)*SIN(DC1)
			YP(I)=-XP(I)*SIN(DC1)+YP(I)*COS(DC1)
			XP(I)=TT
		ENDDO
	ENDIF
C Dimensions of grid aperture
	GX=DA*0.5+GBIT
	GY=DR*0.5+GBIT
C Set vector parameters for grid aperture
	PP(1)=ZP(1)
	PP(2)=ZP(2)
	PP(3)=ZP(3)
	PP(4)=XP(1)
	PP(5)=XP(2)
	PP(6)=XP(3)
	IF(HGRID.GT.0.0) THEN
		PP(7)=POS(1)+ZP(1)*PL+CAX(1)*HGRID
		PP(8)=POS(2)+ZP(2)*PL+CAX(2)*HGRID
		PP(9)=POS(3)+ZP(3)*PL+CAX(3)*HGRID
	ELSE
		PP(7)=POS(1)+ZP(1)*PL+CAX(1)*WALL
		PP(8)=POS(2)+ZP(2)*PL+CAX(2)*WALL
		PP(9)=POS(3)+ZP(3)*PL+CAX(3)*WALL
	ENDIF
        PP(10)=-GX
        PP(11)=-GY
        PP(12)=GX
        PP(13)=GY
C Find surface index of grid aperture and set that surface element
	KSUR=PC(10)+1
	CALL SRT_SETF(KSUR,1,13,PP,0,0,0,KSUR+1,ISTAT)
C Set vector parameters for grid exit aperture
	IF(HGRID.GT.0.0) THEN
		PP(7)=POS(1)+ZP(1)*PL+CAX(1)*(HGRID-GDEP)
		PP(8)=POS(2)+ZP(2)*PL+CAX(2)*(HGRID-GDEP)
		PP(9)=POS(3)+ZP(3)*PL+CAX(3)*(HGRID-GDEP)
C set exit aperture for a plate
	ELSE
		PP(7)=POS(1)+ZP(1)*PL+CAX(1)*WALL/2.0
		PP(8)=POS(2)+ZP(2)*PL+CAX(2)*WALL/2.0
		PP(9)=POS(3)+ZP(3)*PL+CAX(3)*WALL/2.0
	ENDIF
C Set grid exit surface element
	KSUR=PC(10)+2
	CALL SRT_SETF(KSUR,1,13,PP,0,0,0,KSUR+1,ISTAT)
C Dimensions of pore
	HX=DA*0.5
	HY=DR*0.5
	HL=PL*0.5
C Get surface quality of pore reflecting wall
	IQ=INT(PC(7))
C Set vector parameters for pore aperture
	PP(1)=ZP(1)
	PP(2)=ZP(2)
	PP(3)=ZP(3)
	PP(4)=XP(1)
	PP(5)=XP(2)
	PP(6)=XP(3)
	PP(7)=POS(1)+ZP(1)*PL
	PP(8)=POS(2)+ZP(2)*PL
	PP(9)=POS(3)+ZP(3)*PL
        PP(10)=-HX
        PP(11)=-HY
        PP(12)=HX
        PP(13)=HY
C Set pore aperture and set that surface element
	KSUR=PC(10)+3
	CALL SRT_SETF(KSUR,1,13,PP,0,0,0,KSUR+1,ISTAT)
C Set sides of pore - local x is along pore axis
C This is an absorbing wall
	PP(1)=XP(1)
	PP(2)=XP(2)
	PP(3)=XP(3)
	PP(4)=ZP(1)
	PP(5)=ZP(2)
	PP(6)=ZP(3)
	PP(7)=POS(1)+ZP(1)*HL-XP(1)*HX
	PP(8)=POS(2)+ZP(2)*HL-XP(2)*HX
	PP(9)=POS(3)+ZP(3)*HL-XP(3)*HX
	PP(10)=-HL
        PP(11)=-HY
        PP(12)=HL
        PP(13)=HY
	CALL SRT_SETF(KSUR+1,5,13,PP,0,0,KSUR+1,-1,ISTAT)
C This is the reflecting surface modelled as a cylinder
	PP(1)=ZP(1)
	PP(2)=ZP(2)
	PP(3)=ZP(3)
	PP(4)=XP(1)
	PP(5)=XP(2)
	PP(6)=XP(3)
	PP(7)=POS(1)+ZP(1)*HL-YP(1)*HY+YP(1)*RCURP
	PP(8)=POS(2)+ZP(2)*HL-YP(2)*HY+YP(2)*RCURP
	PP(9)=POS(3)+ZP(3)*HL-YP(3)*HY+YP(3)*RCURP
	IF(ICURV.EQ.0.or.ICURV.eq.5.or.ICURV.eq.8) THEN
C Axial profile linear (conical approximation)
		PP(10)=0.0
	ELSEIF(ICURV.EQ.1.OR.ICURV.EQ.4.OR.ICURV.EQ.6) THEN
C Axial profile curved (Wolter I)
		PP(10)=-TAN(THETAP)**2
	ELSEIF(ICURV.EQ.2.or.ICURV.eq.7) THEN
C Axial profiles curved-plane
		PP(10)=-2.0D0*TAN(THETAP)**2
	ELSEIF(ICURV.EQ.9) THEN
C Axial profiles -1/2, +5/2
		PP(10)=0.5*TAN(THETAP)**2
	ELSEIF(ICURV.EQ.10) THEN
C Axial profiles 1/2, 3/2
		PP(10)=-0.5*TAN(THETAP)**2
	ELSEIF(ICURV.EQ.11) THEN
C Axial profiles -1, 3 
		PP(10)=TAN(THETAP)**2
	ELSEIF(ICURV.EQ.12) THEN
C Axial profiles -1.5, 3.5
		PP(10)=1.5D0*TAN(THETAP)**2
	ELSEIF(ICURV.EQ.13) THEN
C Axial profiles -2, 4
		PP(10)=2.0D0*TAN(THETAP)**2
	ELSE
C Axial profile constant curvature
		PP(10)=-4.3D-4*RCURP
	ENDIF
	PP(11)=0.0
	PP(12)=RCURP**2
	PP(13)=-HL
        PP(14)=-HX
        PP(15)=HL
        PP(16)=HX
	CALL SRT_SETF(KSUR+4,26,16,PP,0,IQ,KSUR+1,-1,ISTAT)
C This is an absorbing wall
	PP(1)=-XP(1)
	PP(2)=-XP(2)
	PP(3)=-XP(3)
	PP(4)=ZP(1)
	PP(5)=ZP(2)
	PP(6)=ZP(3)
	PP(7)=POS(1)+ZP(1)*HL+XP(1)*HX
	PP(8)=POS(2)+ZP(2)*HL+XP(2)*HX
	PP(9)=POS(3)+ZP(3)*HL+XP(3)*HX
	PP(10)=-HL
        PP(11)=-HY
        PP(12)=HL
        PP(13)=HY
	CALL SRT_SETF(KSUR+5,5,13,PP,0,0,KSUR+1,-1,ISTAT)
C This is an absorbing wall
	PP(1)=-YP(1)
	PP(2)=-YP(2)
	PP(3)=-YP(3)
	PP(4)=ZP(1)
	PP(5)=ZP(2)
	PP(6)=ZP(3)
	PP(7)=POS(1)+ZP(1)*HL+YP(1)*HY
	PP(8)=POS(2)+ZP(2)*HL+YP(2)*HY
	PP(9)=POS(3)+ZP(3)*HL+YP(3)*HY
	PP(10)=-HL
        PP(11)=-HX
        PP(12)=HL
        PP(13)=HX
	CALL SRT_SETF(KSUR+6,5,13,PP,0,0,KSUR+1,-1,ISTAT)
C Now 2nd surface part of pore
	HL=PL*0.5/GMOD
C Find distance for vertex of 2nd surface from Principal plane
	Z=RCURH/TAN(THETAP*2.0+THETAH)
C Find direction of 2nd pore axis
	DBC=DB2*COS(THETA)*Z
	DBS=DB2*SIN(THETA)*Z
	ZP(1)=Z*CAX(1)+(X-DX+DBC+DV)*CAR(1)+(Y-DY+DBS)*VY(1)
	ZP(2)=Z*CAX(2)+(X-DX+DBC+DV)*CAR(2)+(Y-DY+DBS)*VY(2)
	ZP(3)=Z*CAX(3)+(X-DX+DBC+DV)*CAR(3)+(Y-DY+DBS)*VY(3)
	CALL SRT_VNRM(ZP,ISTAT)
C Shift pore position to give a small gap between 1st and 2nd pore
	POS(1)=POS(1)-PC(6)*CAX(1)
	POS(2)=POS(2)-PC(6)*CAX(2)
	POS(3)=POS(3)-PC(6)*CAX(3)
C Now set up axes tangential to surface at pore centre
C Note that there is no pore centred at the origin so can use the cross product
	CALL SRT_VCRS(CAX,ZP,XP)
	CALL SRT_VNRM(XP,ISTAT)
	CALL SRT_VCRS(ZP,XP,YP)
C Apply out-of-plane figure error by rotating pore
	IF(DC2.NE.0.0) THEN
		DO I=1,3
			TT=XP(I)*COS(DC2)+YP(I)*SIN(DC2)
			YP(I)=-XP(I)*SIN(DC2)+YP(I)*COS(DC2)
			XP(I)=TT
		ENDDO
	ENDIF
C Set vector parameters for pore aperture
	PP(1)=ZP(1)
	PP(2)=ZP(2)
	PP(3)=ZP(3)
	PP(4)=XP(1)
	PP(5)=XP(2)
	PP(6)=XP(3)
	PP(7)=POS(1)
	PP(8)=POS(2)
	PP(9)=POS(3)
        PP(10)=-HX
        PP(11)=-HY
        PP(12)=HX
        PP(13)=HY
C Find surface index of 2nd pore aperture and set that surface element
	KSUR=PC(10)+8
	CALL SRT_SETF(KSUR,1,13,PP,0,0,0,KSUR+1,ISTAT)
C Set sides of pore - local x is along pore axis
C This is an absorbing wall
	PP(1)=XP(1)
	PP(2)=XP(2)
	PP(3)=XP(3)
	PP(4)=ZP(1)
	PP(5)=ZP(2)
	PP(6)=ZP(3)
	PP(7)=POS(1)-ZP(1)*HL-XP(1)*HX
	PP(8)=POS(2)-ZP(2)*HL-XP(2)*HX
	PP(9)=POS(3)-ZP(3)*HL-XP(3)*HX
	PP(10)=-HL
        PP(11)=-HY
        PP(12)=HL
        PP(13)=HY
	CALL SRT_SETF(KSUR+1,5,13,PP,0,0,KSUR+1,-1,ISTAT)
C This is reflecting surface modelled as a cylinder
	PP(1)=ZP(1)
	PP(2)=ZP(2)
	PP(3)=ZP(3)
	PP(4)=XP(1)
	PP(5)=XP(2)
	PP(6)=XP(3)
	PP(7)=POS(1)-ZP(1)*HL-YP(1)*HY+YP(1)*RCURH
	PP(8)=POS(2)-ZP(2)*HL-YP(2)*HY+YP(2)*RCURH
	PP(9)=POS(3)-ZP(3)*HL-YP(3)*HY+YP(3)*RCURH
	IF(ICURV.EQ.0.OR.ICURV.EQ.2.or.ICURV.eq.5.or.ICURV.eq.7) THEN
C Axial profile linear (conical approximation or curved-plane)
		PP(10)=0.0
	ELSEIF(ICURV.EQ.1.OR.ICURV.EQ.4.OR.ICURV.EQ.6) THEN
C Axial profile Wolter I
		PP(10)=-TAN(THETAH)**2
	ELSEIF(ICURV.EQ.8) THEN
		PP(10)=-TAN(THETAP)**2-TAN(THETAH)**2
	ELSEIF(ICURV.EQ.9) THEN
		PP(10)=-2.5D0*TAN(THETAP)**2
	ELSEIF(ICURV.EQ.10) THEN
		PP(10)=-1.5D0*TAN(THETAP)**2
	ELSEIF(ICURV.EQ.11) THEN
		PP(10)=-3.0D0*TAN(THETAP)**2
	ELSEIF(ICURV.EQ.12) THEN
		PP(10)=-3.5D0*TAN(THETAP)**2
	ELSEIF(ICURV.EQ.13) THEN
		PP(10)=-4.0D0*TAN(THETAP)**2
	ELSE
C Axial profile constant curvature
		PP(10)=-1.56D-7*RCURH
	ENDIF
	PP(11)=0.0
	PP(12)=RCURH**2
	PP(13)=-HL
        PP(14)=-HX
        PP(15)=HL
        PP(16)=HX
	CALL SRT_SETF(KSUR+2,26,16,PP,0,IQ,KSUR+1,-1,ISTAT)
C This is an absorbing wall
	PP(1)=-XP(1)
	PP(2)=-XP(2)
	PP(3)=-XP(3)
	PP(4)=ZP(1)
	PP(5)=ZP(2)
	PP(6)=ZP(3)
	PP(7)=POS(1)-ZP(1)*HL+XP(1)*HX
	PP(8)=POS(2)-ZP(2)*HL+XP(2)*HX
	PP(9)=POS(3)-ZP(3)*HL+XP(3)*HX
	PP(10)=-HL
        PP(11)=-HY
        PP(12)=HL
        PP(13)=HY
	CALL SRT_SETF(KSUR+3,5,13,PP,0,0,KSUR+1,-1,ISTAT)
C This is an absorbing wall
	PP(1)=-YP(1)
	PP(2)=-YP(2)
	PP(3)=-YP(3)
	PP(4)=ZP(1)
	PP(5)=ZP(2)
	PP(6)=ZP(3)
	PP(7)=POS(1)-ZP(1)*HL+YP(1)*HY
	PP(8)=POS(2)-ZP(2)*HL+YP(2)*HY
	PP(9)=POS(3)-ZP(3)*HL+YP(3)*HY
	PP(10)=-HL
        PP(11)=-HX
        PP(12)=HL
        PP(13)=HX
	CALL SRT_SETF(KSUR+4,5,13,PP,0,0,KSUR+1,-1,ISTAT)
	END
