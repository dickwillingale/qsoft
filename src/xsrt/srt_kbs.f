*+SRT_KBS    Dynamically set up a rectangular slot for K-B stack
	SUBROUTINE SRT_KBS(DIR,HPOS,CAX,CAR,CVR,PC,IDEF,ISTAT)
	IMPLICIT NONE
	DOUBLE PRECISION DIR(3),HPOS(3),CAX(3),CAR(3),CVR(3),PC(50)
	INTEGER IDEF(2),ISTAT
*DIR	input	direction of ray
*HPOS	input	position of hit on aperture plane reference surface
*CAX    input   axis of optic
*CAR    input   reference direction in aperture
*CVR    input   centre of curvature of front aperture 
*PC     input   remaining parameters
*IDEF   input   deformation indices
*ISTAT  in/out  returned status
* Parameters for slots are:
* PC(1) radius of spherical aperture surface if 1st stack 2*FLEN+PLMAX
*					     if 2nd stack 2*FLEN
* PC(2) outer radius of aperture
* PC(3) focal length FLEN (-ve for 2nd stack)
* PC(4) pitch of slots
* PC(5) wall thickness of slots
* PC(6) minimum axial slot length (PLMIN)
* PC(7) maximum axial slot length (PLMAX)
* PC(8) local x of ray intersection
* PC(9) local y of ray intersection
* PC(10) this surface index
* PC(11) surface quality
* Note: -ve focal length indicates 2nd stack
*-Author Dick Willingale 2008-Nov-20
	DOUBLE PRECISION PI
	PARAMETER (PI=3.14159265359)
	DOUBLE PRECISION VY(3),X,Y,XP(3),YP(3),ZP(3),POS(3)
	DOUBLE PRECISION XC,YC,PLMAX,MCAX(3),MCVR(3),DELT,RXY
	DOUBLE PRECISION PP(16),HX,HY,HL,THETA,PHI,RSPHR
	DOUBLE PRECISION DX,DY,DZ,DB,DC,DDX,DDY,RAN(2)
	DOUBLE PRECISION XX,FL,YY,PL,ELH,ELW,ZS,SHWID,SHLEN,THETAG
	INTEGER I,KSUR,IQ,NSLOT,ISTACK,ITRY
      	DOUBLE PRECISION RM,PM,TM,WM,HM,CM,GM
	INTEGER MODHIT,ICURV
	LOGICAL HIT
C Get radius of curvature of spherial surface
	RSPHR=PC(1)
C Get focal length and set stack number using sign
	FL=PC(3)
	if(FL.GT.0.0) THEN
		ISTACK=1
	else
		ISTACK=2
		FL=-FL
	endif
C Find module index
	CALL SRT_SFINDMOD(RSPHR,PC(8),PC(9),MODHIT,RM,PM,TM,WM,HM,PL,
     +	CM,GM,XX,YY)
     	if(MODHIT.GT.0) then
		THETA=PM+TM
		ELH=HM*0.5
		ELW=WM*0.5
C Set slot width
		SHWID=(PC(4)-PC(5))*0.5
C Centre of slot reflecting wall in local module aperture coordinates XX,YY
		IF(ISTACK.EQ.1) THEN
			NSLOT=INT((YY+ELH)/PC(4))
			YY=(NSLOT+0.5)*PC(4)-ELH+SHWID
			XX=0.0
			SHLEN=ELW
		ELSE
			NSLOT=INT((XX+ELW)/PC(4))
			XX=(NSLOT+0.5)*PC(4)-ELW+SHWID
			YY=0.0
			SHLEN=ELH
		ENDIF
C Calculate position of module centre in full aperture coordinates
		XC=RM*COS(PM)
		YC=RM*SIN(PM)
C Transform back to full aperture coordinates
		X=XX*COS(THETA)-YY*SIN(THETA)+XC
		Y=XX*SIN(THETA)+YY*COS(THETA)+YC
C Set slot rotation for 2nd stack
		IF(ISTACK.NE.1) THEN
			THETA=THETA-PI*0.5
		ENDIF
		PLMAX=PC(7)
		THETAG=SHWID*2.0/PLMAX
	ELSE
C Ray outside module so set slot to be daft position such that ray misses
		X=1.e5
		Y=1.e5
		THETA=0.0
		THETAG=0.0
	endif
C X,Y are centre of slot reflecting wall in local full aperture coordinates
C Find Y reference axis at local origin on aperture
	CALL SRT_VCRS(CAX,CAR,VY)
	CALL SRT_VNRM(VY,ISTAT)
C Calculate position of centre of slot reflecting wall
		DO I=1,3
		POS(I)=CVR(I)+CAX(I)*SQRT(RSPHR**2-X**2-Y**2)
		POS(I)=POS(I)+CAR(I)*X+VY(I)*Y
	ENDDO
C Find module axis
C	DO I=1,3
C		MCAX(I)=CAX(I)*SQRT(RSPHR**2-XC**2-YC**2)
C		MCAX(I)=MCAX(I)+CAR(I)*XC+VY(I)*YC
C	ENDDO
C	CALL SRT_VNRM(MCAX,ISTAT)
C Find module vertex
C	PLMAX=PC(7)
C	IF(ISTACK.EQ.1) THEN
C		DELT=-PLMAX*0.5
C	ELSE
C		DELT=-PLMAX*0.5
C	ENDIF
C 
	DO I=1,3
C		MCVR(I)=CVR(I)+MCAX(I)*DELT
		MCVR(I)=CVR(I)
	ENDDO
C Deformations
	IF(IDEF(1).GT.0.AND.MODHIT.GT.0) THEN
C DX is effective displacement of module vertex in direction CAR
C DY is effective displacement of module vertex in direction CAX X CAR
C DZ is effective displacement of module vertex in axial direction CAX
C DB mean value for Gaussian in-plane figure errors radians
C DC mean value for Gaussian out-of-plane figure errors radians
C Note that the derivatives are not used
		IDEF(2)=1
		CALL SRT_IDFRM(IDEF,MODHIT,1,DX,DDX,DDY,ISTAT)
		IDEF(2)=2
		CALL SRT_IDFRM(IDEF,MODHIT,1,DY,DDX,DDY,ISTAT)
		IDEF(2)=3
		CALL SRT_IDFRM(IDEF,MODHIT,1,DZ,DDX,DDY,ISTAT)
		IDEF(2)=4
		CALL SRT_IDFRM(IDEF,MODHIT,1,DB,DDX,DDY,ISTAT)
		IDEF(2)=5
		CALL SRT_IDFRM(IDEF,MODHIT,1,DC,DDX,DDY,ISTAT)
		CALL SYS_GAUSS(2,RAN,0.0D0,1.0D0,ISTAT)
		DB=RAN(1)*DB
		DC=RAN(2)*DC
	ELSE
		DX=0.0
		DY=0.0
		DB=0.0
		DC=0.0
	ENDIF
C Find direction of slot axis
	DO I=1,3
		ZP(I)=POS(I)-MCVR(I)-DX*CAR(I)-DY*VY(I)-DZ*CAX(I)
	ENDDO
	CALL SRT_VNRM(ZP,ISTAT)
C Now set up axes tangential to spherical surface at slot centre
        CALL SRT_VCRS(ZP,CAR,YP)
        CALL SRT_VNRM(YP,ISTAT)
	CALL SRT_VCRS(YP,ZP,XP)
C Rotate reference axes to align with module
C Include out-of-plane figure errors by including slot rotation error
	THETA=THETA+DC
	DO I=1,3
		XP(I)=COS(THETA)*XP(I)+SIN(THETA)*YP(I)
	ENDDO
        CALL SRT_VCRS(ZP,XP,YP)
C Include in-plane figure errors by shifting axis in YP direction
	DO I=1,3
		ZP(I)=ZP(I)*SQRT(1.0-DB**2)+YP(I)*DB
	ENDDO
        CALL SRT_VNRM(ZP,ISTAT)
        CALL SRT_VCRS(ZP,XP,YP)
C Set dimensions of slot aperture
	HX=SHLEN
	HY=SHWID
C Shift position from reflecting wall to slot centre
	DO I=1,3
		POS(I)=POS(I)-YP(I)*HY
	ENDDO
C Set vector parameters for slot aperture
	PP(1)=ZP(1)
	PP(2)=ZP(2)
	PP(3)=ZP(3)
	PP(4)=XP(1)
	PP(5)=XP(2)
	PP(6)=XP(3)
	PP(7)=POS(1)
	PP(8)=POS(2)
	PP(9)=POS(3)
C
        PP(10)=-HX
        PP(11)=-HY
        PP(12)=HX
        PP(13)=HY
C Find surface index of slot aperture and set that surface element
	KSUR=PC(10)+1
	CALL SRT_SETF(KSUR,1,13,PP,0,0,0,KSUR+1,ISTAT)
C set half axial length of slot
	HL=PL*0.5
C Get surface quality of slot walls
	IQ=INT(PC(11))
C local x is along slot axis
C absorbing wall of slot
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
C absorbing wall of slot
	PP(1)=YP(1)
	PP(2)=YP(2)
	PP(3)=YP(3)
	PP(4)=ZP(1)
	PP(5)=ZP(2)
	PP(6)=ZP(3)
	PP(7)=POS(1)-ZP(1)*HL-YP(1)*HY
	PP(8)=POS(2)-ZP(2)*HL-YP(2)*HY
	PP(9)=POS(3)-ZP(3)*HL-YP(3)*HY
	PP(10)=-HL
        PP(11)=-HX
        PP(12)=HL
        PP(13)=HX
	CALL SRT_SETF(KSUR+2,5,13,PP,0,0,KSUR+1,-1,ISTAT)
C absorbing wall of slot
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
C This is the reflecting wall of slot modelled as a cylinder with very
C large radius of curvature
	RXY=1.e6
	PP(1)=ZP(1)
	PP(2)=ZP(2)
	PP(3)=ZP(3)
	PP(4)=XP(1)
	PP(5)=XP(2)
	PP(6)=XP(3)
	PP(7)=POS(1)-ZP(1)*HL+YP(1)*HY-YP(1)*RXY
	PP(8)=POS(2)-ZP(2)*HL+YP(2)*HY-YP(2)*RXY
	PP(9)=POS(3)-ZP(3)*HL+YP(3)*HY-YP(3)*RXY
	ICURV=0
	IF(ICURV.EQ.0) THEN
C Axial profile linear (flat plate)
		PP(10)=0.0
	ELSE IF(ICURV.EQ.1) THEN
C Axial profile curved 
		PP(10)=2.0D0*TAN(THETAG)**2
	ENDIF
	PP(11)=0.0
	PP(12)=RXY**2
	PP(13)=-HL
        PP(14)=-HX
        PP(15)=HL
        PP(16)=HX
	CALL SRT_SETF(KSUR+4,26,16,PP,0,IQ,KSUR+1,-1,ISTAT)
	END
