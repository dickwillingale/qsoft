*+SRT_SQMPOARR    Dynamically set up pore in a square pore MPO array
	SUBROUTINE SRT_SQMPOARR(DIR,HPOS,CAX,CAR,CVR,PC,IDEF,ISTAT)
	IMPLICIT NONE
	DOUBLE PRECISION DIR(3),HPOS(3),CAX(3),CAR(3),CVR(3),PC(50)
	INTEGER IDEF(2),ISTAT
*DIR	input	direction of ray
*HPOS	input	position of hit on spherical reference surface
*CAX    input   axis of sphere, from centre to local origin on surface
*CAR    input   reference direction
*CVR    input   centre of sphere
*PC     input   radius of sphere and other paramters (see below)
*IDEF   input   deformation indices
*ISTAT  in/out  returned status
* Parameters are:
* PC(1) is the radius of curvature of the spherical MPO array
* PC(2) is the full aperture half width of the MPO array
* PC(3) void (spare)
* PC(4) void (spare)
* PC(5) void (spare)
* PC(6) void (spare)
* PC(7) void (spare)
* PC(8) local x of ray intersection
* PC(9) local y of ray intersection
* PC(10) this surface index
* All pore parameters are set for the individual modules in the MPO array 
*-Author Dick Willingale 2017-Oct-23
	DOUBLE PRECISION PI,SMALL
	PARAMETER (PI=3.14159265359,SMALL=1.D-4)
	DOUBLE PRECISION VY(3),X,Y,Z,XP(3),YP(3),ZP(3),RXY,CTH,STH
	DOUBLE PRECISION T,POS(3)
	DOUBLE PRECISION MPC(3),RPOS,DA,DC,XPP(3),YPP(3),DB,RAN(3)
	DOUBLE PRECISION PP(14),HX,HY,HL,THETA,DTHETA,DX,DY,DT,DDX,DDY
        DOUBLE PRECISION DE,DS
	DOUBLE PRECISION PWAFFLE,X0,SS,Y0,TPORE,RP,RPS,RLO,RHI
	DOUBLE PRECISION GAP,FIBRE,XFIB,YFIB,RFIB,SFIB
	DOUBLE PRECISION ATHETA,T1,T2,T3,T4,PIBY4,ABIT,S45,PL,COST,TGAP
        DOUBLE PRECISION SSHR,DSHR,TFIB,CRFACT,BAND,WFR,DIQ,SP
	DOUBLE PRECISION XMOD,YMOD,TMOD,WMOD,HMOD,LMOD,RMOD,PITCH,WALL
	DOUBLE PRECISION SYS_DRAND,PPP
	DOUBLE PRECISION RPLATE,TPLATE,XPLATE,YPLATE,DDD,RMP,WPL,QQQ,TTT
        DOUBLE PRECISION BU,BV,BZ
	INTEGER IXF,IYF,IIF
	INTEGER IX,IY,I,KSUR,IQ,NRAD,NWAFFLE,IP
	INTEGER NANUL,NSEC,IAPER
	LOGICAL HIT
C Set spherical radius of the MPO array
	RP=PC(1)
C Use local intersection position to find centre of nearest pore
C The impact point is on the spherical reference surface
C Here we find the rectangular aperture containing the ray
C PC(8) and PC(9) are the local intersection point on spherical surface
	CALL SRT_FINDMPO(PC(8),PC(9),IAPER,XMOD,YMOD,TMOD,
     +	WMOD,HMOD,PL,RMOD,FIBRE,PITCH,WALL,DIQ,BU,BV,BZ,SP,XPLATE,YPLATE)
C Set surface quality for this MPO IQ as integer
	IQ=INT(DIQ)
C Get half width of plate
	WPL=WMOD*0.5
C Set polar coordinates on plate
	RPLATE=SQRT(XPLATE**2+YPLATE**2)
	TPLATE=ATAN2(YPLATE,XPLATE)
	IF(IAPER.GT.0) THEN
C Find centre of pore on spherical surface
		IF(PC(8).GT.0.0) THEN
			X=(INT(PC(8)/PITCH)+0.5)*PITCH
		ELSE
			X=(INT(PC(8)/PITCH)-0.5)*PITCH
		ENDIF
		IF(PC(9).GT.0.0) THEN
			Y=(INT(PC(9)/PITCH)+0.5)*PITCH
		ELSE
			Y=(INT(PC(9)/PITCH)-0.5)*PITCH
		ENDIF
		RXY=SQRT(X**2+Y**2)
C set NRAD for cartesian packing so that dont loose central pores
		NRAD=1
C Now find position within multifibre structure
		XFIB=MOD(ABS(XPLATE),FIBRE)-FIBRE*0.5
		YFIB=MOD(ABS(YPLATE),FIBRE)-FIBRE*0.5
		RFIB=SQRT(XFIB**2+YFIB**2)
		TFIB=ATAN2(YFIB,XFIB)
C Generate unique integer for each multi-fibre
		IXF=INT((XPLATE+WMOD*0.5)/FIBRE)
		IYF=INT((YPLATE+HMOD*0.5)/FIBRE)
		IIF=ABS(IXF+IYF*10000)
	ELSE
C Pore centre outside aperture so move pore so that ray misses
		X=PC(2)*2.0
		Y=PC(2)*2.0
	ENDIF
C X and Y are now the position of the pore axis on spherical surface
C Find Y reference axis at local origin on spherical surface
	CALL SRT_VCRS(CAX,CAR,VY)
	CALL SRT_VNRM(VY,ISTAT)
C Find frame surface normal at pore centre
	Z=SQRT(RP**2-X**2-Y**2)
	ZP(1)=Z*CAX(1)+X*CAR(1)+Y*VY(1)
	ZP(2)=Z*CAX(2)+X*CAR(2)+Y*VY(2)
	ZP(3)=Z*CAX(3)+X*CAR(3)+Y*VY(3)
	CALL SRT_VNRM(ZP,ISTAT)
C Find position of pore centre
	POS(1)=CVR(1)+ZP(1)*RP
	POS(2)=CVR(2)+ZP(2)*RP
	POS(3)=CVR(3)+ZP(3)*RP
C Set effective displacement of centre of sphere
C to allow for difference in RoC of frame and MPO
C DX in direction CAR
C DY in direction NORM X CAR
	DX=(RMOD-RP)*XPLATE/RMOD
	DY=(RMOD-RP)*YPLATE/RMOD
C Add MPO bias angle shift
        DX=DX+BU*RP
        DY=DY+BV*RP
C Set shear error
    	DSHR=0.0
    	DC=0.0
C Set pore rotation TPORE for cartesian plate
        TPORE=0.0
	IF(IDEF(1).GT.0) THEN
C Deformations specified
C Note that the derivatives are not used
C
C DS scaling of intrinsic slump, 1.0 gives model amplitude
		IDEF(2)=1
		CALL SRT_IDFRM(IDEF,IAPER,1,DS,DDX,DDY,ISTAT)
C DDD is the reduction factor applied in the corners of the plate
C QQQ is profile of the thermoelasic errors
                if(RPLATE.LT.WPL) then
		        DDD=1.0
			QQQ=RPLATE/WPL
		else
			RMP=WPL*SQRT(2.0)
			DDD=(RMP-RPLATE)/(RMP-WPL)
			QQQ=DDD
		endif
C Intrinsic slump tilt errors
		DX=DX+DS*COS(TPLATE)*DDD*RPLATE**2/2.0/RMOD
		DY=DY+DS*SIN(TPLATE)*DDD*RPLATE**2/2.0/RMOD
C DE +ve maximum amplitude of thermoelastic axial pointing errors radians
C DE -ve chess-board pattern of tilt errors
		IDEF(2)=2
		CALL SRT_IDFRM(IDEF,IAPER,1,DE,DDX,DDY,ISTAT)
		IF(DE.GE.0.0) THEN
			DX=DX+DE*COS(TPLATE)*QQQ*RMOD
			DY=DY+DE*SIN(TPLATE)*QQQ*RMOD
		ELSE
			TTT=PI/5.0
			PPP=SIN(TTT*XPLATE+TTT*YPLATE/4.0)
     +			*SIN(TTT*YPLATE+TTT*XPLATE/4.0)
			PPP=SIGN(ABS(PPP)**2.0,PPP)
			DX=DX+DE*PPP*RMOD*SIN(RPLATE*TTT)
			DY=DY+DE*PPP*RMOD*COS(RPLATE*TTT)
		ENDIF
C DA is the pore axis rotation error
		IDEF(2)=3
		CALL SRT_IDFRM(IDEF,IAPER,1,DA,DDX,DDY,ISTAT)
		TPORE=TPORE+DA
C DC Shear errors within the multifibre structure
		IDEF(2)=4
		CALL SRT_IDFRM(IDEF,IAPER,1,DC,DDX,DDY,ISTAT)
		SFIB=PITCH*3.0
		DC=DC*EXP(-RFIB**2/2.0/SFIB**2)*SIGN(1.0D0,TAN(TFIB))
C Include shear errors intrinsic to the slump using DS scaling
		DC=DC+DS*RPLATE*ABS(SIN(2.0*TPLATE))/RP*(PITCH-WALL)*
     +		SIGN(1.0D0,TAN(TPLATE))*DDD
		DSHR=DC/(PITCH-WALL)
C DB +ve rms amplitude Gaussian pore tilt/figure errors (2.36*DB=FWHM)
C DB -ve amplitude Cauchy pore tilt/figure errors (2*DB=FWHM)
		IDEF(2)=5
		CALL SRT_IDFRM(IDEF,IAPER,1,DB,DDX,DDY,ISTAT)
		IF(DB.GT.0.0) THEN
			CRFACT=RPLATE/WPL
	                CALL SYS_GAUSS(3,RAN,0.0D0,1.0D0,ISTAT)
			DX=DX+RAN(1)*DB*CRFACT*RMOD
			DX=DX+RAN(2)*DB*CRFACT*RMOD
		ELSE
                        DX=DX+DB*TAN((SYS_DRAND()-0.5)*PI)*RMOD
                        DX=DX+DB*TAN((SYS_DRAND()-0.5)*PI)*RMOD
		ENDIF
	ENDIF
C Apply displacements to change effective centre of sphere
	ZP(1)=Z*CAX(1)+(X-DX)*CAR(1)+(Y-DY)*VY(1)
	ZP(2)=Z*CAX(2)+(X-DX)*CAR(2)+(Y-DY)*VY(2)
	ZP(3)=Z*CAX(3)+(X-DX)*CAR(3)+(Y-DY)*VY(3)
	CALL SRT_VNRM(ZP,ISTAT)
C Now set up axes tangential to surface at pore centre
C check ray direction to see if pore inside of outside spherical surface
	CALL SRT_VDOT(ZP,DIR,COST)
	IF(COST.LT.0.0) THEN
C Reverse direction to align with ray
		ZP(1)=-ZP(1)
		ZP(2)=-ZP(2)
		ZP(3)=-ZP(3)
	ENDIF
C Note that there is no pore centred at the origin so can use the cross product
	CALL SRT_VCRS(CAX,ZP,XP)
	CALL SRT_VNRM(XP,ISTAT)
	CALL SRT_VCRS(ZP,XP,YP)
C For cartesian packing set the pore rotation angle using the reference axes
	CALL SRT_VDOT(XP,CAR,CTH)
	CALL SRT_VDOT(XP,VY,STH)
	TPORE=TPORE+ATAN2(STH,CTH)
C rotate pore to correct angle
	DO I=1,3
		XP(I)=XP(I)*COS(TPORE)+YP(I)*SIN(TPORE)
	ENDDO
	CALL SRT_VCRS(ZP,XP,YP)
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
C Set dimensions of pore aperture
	IF(NRAD.GT.0) THEN
		HX=(PITCH-WALL)*0.5
		HY=(PITCH-WALL)*0.5
	ELSE
C If at centre of radial packing then block off
		HX=0.0
		HY=0.0
	ENDIF
	HL=PL*0.5
        PP(10)=-HX
        PP(11)=-HY
        PP(12)=HX
        PP(13)=HY
	PP(14)=DC
C Find surface index of pore aperture and set that surface element
	KSUR=PC(10)+1
	CALL SRT_SETF(KSUR,27,14,PP,0,0,0,KSUR+1,ISTAT)
C Calculate normal of sides including shear
	XPP(1)=XP(1)-YP(1)*DSHR*0.5
	XPP(2)=XP(2)-YP(2)*DSHR*0.5
	XPP(3)=XP(3)-YP(3)*DSHR*0.5
	CALL SRT_VNRM(XPP,ISTAT)
	YPP(1)=YP(1)-XP(1)*DSHR*0.5
	YPP(2)=YP(2)-XP(2)*DSHR*0.5
	YPP(3)=YP(3)-XP(3)*DSHR*0.5
	CALL SRT_VNRM(YPP,ISTAT)
C Calculate side scaling from shear
	SSHR=SQRT(1.0+(DSHR*0.5)**2)
C Set sides of pore - local x is along pore axis
C Side 1
	PP(1)=XPP(1)
	PP(2)=XPP(2)
	PP(3)=XPP(3)
	PP(4)=ZP(1)
	PP(5)=ZP(2)
	PP(6)=ZP(3)
	PP(7)=POS(1)+ZP(1)*HL-XP(1)*HX
	PP(8)=POS(2)+ZP(2)*HL-XP(2)*HX
	PP(9)=POS(3)+ZP(3)*HL-XP(3)*HX
	PP(10)=-HL
        PP(11)=-HY*SSHR
        PP(12)=HL
        PP(13)=HY*SSHR
	CALL SRT_SETF(KSUR+1,5,13,PP,0,IQ,KSUR+1,-1,ISTAT)
C Side 2
	PP(1)=YPP(1)
	PP(2)=YPP(2)
	PP(3)=YPP(3)
	PP(4)=ZP(1)
	PP(5)=ZP(2)
	PP(6)=ZP(3)
	PP(7)=POS(1)+ZP(1)*HL-YP(1)*HY-XP(1)*DC*0.5
	PP(8)=POS(2)+ZP(2)*HL-YP(2)*HY-XP(2)*DC*0.5
	PP(9)=POS(3)+ZP(3)*HL-YP(3)*HY-XP(3)*DC*0.5
	PP(10)=-HL
        PP(11)=-HX*SSHR
        PP(12)=HL
        PP(13)=HX*SSHR
	CALL SRT_SETF(KSUR+2,5,13,PP,0,IQ,KSUR+1,-1,ISTAT)
C Side 3
	PP(1)=-XPP(1)
	PP(2)=-XPP(2)
	PP(3)=-XPP(3)
	PP(4)=ZP(1)
	PP(5)=ZP(2)
	PP(6)=ZP(3)
	PP(7)=POS(1)+ZP(1)*HL+XP(1)*HX
	PP(8)=POS(2)+ZP(2)*HL+XP(2)*HX
	PP(9)=POS(3)+ZP(3)*HL+XP(3)*HX
	PP(10)=-HL
        PP(11)=-HY*SSHR
        PP(12)=HL
        PP(13)=HY*SSHR
	CALL SRT_SETF(KSUR+3,5,13,PP,0,IQ,KSUR+1,-1,ISTAT)
C Side 4
	PP(1)=-YPP(1)
	PP(2)=-YPP(2)
	PP(3)=-YPP(3)
	PP(4)=ZP(1)
	PP(5)=ZP(2)
	PP(6)=ZP(3)
	PP(7)=POS(1)+ZP(1)*HL+YP(1)*HY+XP(1)*DC*0.5
	PP(8)=POS(2)+ZP(2)*HL+YP(2)*HY+XP(2)*DC*0.5
	PP(9)=POS(3)+ZP(3)*HL+YP(3)*HY+XP(3)*DC*0.5
	PP(10)=-HL
        PP(11)=-HX*SSHR
        PP(12)=HL
        PP(13)=HX*SSHR
	CALL SRT_SETF(KSUR+4,5,13,PP,0,IQ,KSUR+1,-1,ISTAT)
	END
