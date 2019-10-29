*+SRT_PORE     Dynamically set up square pore of slumped MCP
	SUBROUTINE SRT_PORE(DIR,HPOS,CAX,CAR,CVR,PC,IDEF,ISTAT)
	IMPLICIT NONE
	DOUBLE PRECISION DIR(3),HPOS(3),CAX(3),CAR(3),CVR(3),PC(50)
	INTEGER IDEF(2),ISTAT
*DIR	input	direction of ray
*HPOS	input	position of hit on spherical reference surface
*CAX    input   axis of sphere, from centre to local origin on surface
*CAR    input   reference direction
*CVR    input   centre of sphere
*PC     input   radius of sphere and limits etc.
*IDEF   input   deformation indices
*ISTAT  in/out  returned status
* Parameters for pores are:
* PC(1) is the radius of the sphere at front of plate
* PC(2) is the aperture radius/half width
* PC(3) pitchx or pitch r
* PC(4) pitchy or pitch theta
* PC(5) ribx or ribr
* PC(6) riby or rib theta
* PC(7) pore length
* PC(8) local x of ray intersection
* PC(9) local y of ray intersection
* PC(10) this surface index
* PC(11) surface quality
* PC(12) pore packing, 1 for cart, 2 for rad, 3 waff, 4 octag, 5 rand, 6 MIXS
*	               7 NFL
* PC(13) minimum pore length
* PC(14) maximum pore length
* PC(15) fibre bundle size
* If the last PC(13) and PC(14) are different then the length of the
* pores is PC(7)*PC(2)/RADIUS and lengths outside range are clipped
* Added IPACK=6 MIXS packing/aperture RW 2007-Jan-18
* Changed IPACK=6 MIXS aperture RW 2012-Jul-12
* Added IPACK=7 NFL aperture RW 2014-Apr-23
* For IPACK=6 the length of the pores is specified using the aperture
* constellation
* Added slumping errors for IPACK=7 RW 2015-Jun-4
*-Author Dick Willingale 2004-Feb-13
* Note direction of ray is used to determine whether pore should be
* inside of outside the spherical surface
	DOUBLE PRECISION PI,SMALL
	PARAMETER (PI=3.14159265359,SMALL=1.D-4)
	DOUBLE PRECISION VY(3),X,Y,Z,XP(3),YP(3),ZP(3),RXY,CTH,STH
	DOUBLE PRECISION T,POS(3)
	DOUBLE PRECISION RNM(3),RPOS,DA,DC,XPP(3),YPP(3),DB,DB1,DB2,RAN(3)
	DOUBLE PRECISION PP(14),HX,HY,HL,THETA,DTHETA,DX,DY,DT,DDX,DDY
        DOUBLE PRECISION DE,DS
	DOUBLE PRECISION PWAFFLE,X0,XX,SS,Y0,TWAFFLE,RP,RPS,RLO,RHI
	DOUBLE PRECISION GAP,FIBRE,XFIB,YFIB,RFIB,SFIB
	DOUBLE PRECISION ATHETA,T1,T2,T3,T4,PIBY4,ABIT,YY,S45,PL,COST,TGAP
        DOUBLE PRECISION SSHR,DSHR,TFIB,CRFACT,BAND,WFR
	DOUBLE PRECISION RMOD,PMOD,TMOD,WMOD,HMOD,LMOD,CM,GM
	DOUBLE PRECISION SYS_DRAND,PPP,ALPHA
        DOUBLE PRECISION SRT_KINGRAND
        EXTERNAL SRT_KINGRAND
	DOUBLE PRECISION RPLATE,TPLATE,XPLATE,YPLATE,DDD,RMP,WPL,QQQ,TTT
	INTEGER IXF,IYF,IIF
	INTEGER IX,IY,I,KSUR,IQ,IPACK,NRAD,NWAFFLE,IP
	INTEGER NANUL,NSEC,IAPER
	LOGICAL HIT
C Get packing type
	IPACK=INT(PC(12))
C Set spherical radius at front of plate
	RP=PC(1)
C Use local intersection position to find centre of nearest pore
	IF(IPACK.EQ.1) THEN
C Cartesian packing
		IF(PC(8).GT.0.0) THEN
			X=(INT(PC(8)/PC(3))+0.5)*PC(3)
		ELSE
			X=(INT(PC(8)/PC(3))-0.5)*PC(3)
		ENDIF
		IF(PC(9).GT.0.0) THEN
			Y=(INT(PC(9)/PC(4))+0.5)*PC(4)
		ELSE
			Y=(INT(PC(9)/PC(4))-0.5)*PC(4)
		ENDIF
		RXY=SQRT(X**2+Y**2)
		IF(RXY.GT.0.0) THEN
			PL=PC(7)*PC(2)/RXY
		ENDIF
		PL=MAX(MIN(PL,PC(14)),PC(13))
		TWAFFLE=0.0
C set NRAD for cartesian packing so that dont loose central pores
		NRAD=1
	ELSEIF(IPACK.EQ.2) THEN
C Radial packing
		RXY=SQRT(PC(8)**2+PC(9)**2)
		NRAD=INT(RXY/PC(3))
		RXY=(NRAD+0.5)*PC(3)
		DTHETA=2.0*PI/(INT(2.0*PI*RXY/PC(4)))
		THETA=ATAN2(PC(9),PC(8))
		IF(THETA.GT.0.0) THEN
			THETA=(INT(THETA/DTHETA)+0.5)*DTHETA
		ELSE
			THETA=(INT(THETA/DTHETA)-0.5)*DTHETA
		ENDIF
		X=RXY*COS(THETA)
		Y=RXY*SIN(THETA)
		IF(RXY.GT.0.0) THEN
			PL=PC(7)*PC(2)/RXY
		ENDIF
		PL=MAX(MIN(PL,PC(14)),PC(13))
		TWAFFLE=0.0
	ELSEIF(IPACK.EQ.3) THEN
C Waffle packing
		NWAFFLE=400
		PWAFFLE=PC(3)*NWAFFLE*2.0*SQRT(2.0)/PI
C Find start of waffle period which contains intersection point
		IF(PC(8).GT.0.0) THEN
			X0=INT(PC(8)/PWAFFLE)*PWAFFLE
		ELSE
			X0=(INT(PC(8)/PWAFFLE)-1)*PWAFFLE
		ENDIF
C Find length across period and convert this radians
		XX=PC(8)-X0
		IF(XX.GT.PWAFFLE*0.5) THEN
			IP=1
			XX=XX-PWAFFLE*0.5
			X0=X0+PWAFFLE*0.5
		ELSE
			IP=-1
		ENDIF
		XX=XX*4.0/PWAFFLE
		SS=SQRT(-XX**2+2.0*XX+1.0)
		TWAFFLE=ATAN2(XX*0.5-0.5+SS*0.5,-XX*0.5+0.5+SS*0.5)
C Convert radians to distance along arc
		TWAFFLE=TWAFFLE*PWAFFLE/(2.0*SQRT(2.0))
C Convert to centre of pore along the arc
		TWAFFLE=(INT(TWAFFLE/PC(3))+0.5)*PC(3)
C Convert back to radians
		TWAFFLE=TWAFFLE*2.0*SQRT(2.0)/PWAFFLE
C Find corresponding distance in X from start of waffle period
		X=PWAFFLE*(SIN(TWAFFLE)-COS(TWAFFLE)+1.0)/4.0
C Find Y offset at this position in waffle period
		Y0=X*TAN(PI*0.25-TWAFFLE*0.5)
C Include start of waffle period in X displacment
		X=X0+X
C Calculate effective position in waffle y layers
		Y=PC(9)+IP*Y0
C Find centre of pore in Y
		IF(Y.GT.0.0) THEN
			Y=(INT(Y/PC(4))+0.5)*PC(4)
		ELSE
			Y=(INT(Y/PC(4))-0.5)*PC(4)
		ENDIF
C Put back the y layer offset
		Y=Y-IP*Y0
C Change TWAFFLE to the required pore rotation angle
		TWAFFLE=(TWAFFLE-PI*0.25)*IP
C Set NRAD so dont loose central pores
		NRAD=1
C Calculate radius
		RXY=SQRT(X**2+Y**2)
		IF(RXY.GT.0.0) THEN
			PL=PC(7)*PC(2)/RXY
		ENDIF
		PL=MAX(MIN(PL,PC(14)),PC(13))
	ELSEIF(IPACK.EQ.4) THEN
C Ocatagonal sector packing scheme
		T1=PI/8.0
		T2=T1*3.0
		T3=T1*5.0
		T4=T1*7.0
		PIBY4=PI/4.0
		S45=SIN(PIBY4)
		THETA=ATAN2(PC(9),PC(8))
		ATHETA=ABS(THETA)
		ABIT=PC(3)*0.5
		IF((ATHETA.GT.T1.AND.ATHETA.LE.T2).OR.
     +		(ATHETA.GT.T3.AND.ATHETA.LE.T4)) THEN
			IF(PC(8).GT.0.0) THEN
				X=(INT(PC(8)/PC(3))+0.5)*PC(3)
			ELSE
				X=(INT(PC(8)/PC(3))-0.5)*PC(3)
			ENDIF
			IF(PC(9).GT.0.0) THEN
				Y=(INT(PC(9)/PC(4))+0.5)*PC(4)
			ELSE
				Y=(INT(PC(9)/PC(4))-0.5)*PC(4)
			ENDIF
			IF((ABS(Y)-ABS(X)*TAN(T1).LT.ABIT).OR.
     +			(ABS(X)-ABS(Y)*TAN(T1).LT.ABIT)) THEN
C Block pores in region between sectors
				NRAD=0
			ELSE
				NRAD=1
			ENDIF
C set pore rotation angle
			TWAFFLE=ATAN2(Y,X)
		ELSE
C Switch to coodinate system at 45 degrees
			XX=PC(8)*S45+PC(9)*S45
			YY=-PC(8)*S45+PC(9)*S45
			IF(XX.GT.0.0) THEN
				XX=(INT(XX/PC(3))+0.5)*PC(3)
			ELSE
				XX=(INT(XX/PC(3))-0.5)*PC(3)
			ENDIF
			IF(YY.GT.0.0) THEN
				YY=(INT(YY/PC(4))+0.5)*PC(4)
			ELSE
				YY=(INT(YY/PC(4))-0.5)*PC(4)
			ENDIF
			IF((ABS(YY)-ABS(XX)*TAN(T1).LT.ABIT).OR.
     +			(ABS(XX)-ABS(YY)*TAN(T1).LT.ABIT)) THEN
C Block pores in region between sectors
				NRAD=0
			ELSE
				NRAD=1
			ENDIF
C rotate back
			X=XX*S45-YY*S45
			Y=XX*S45+YY*S45
C set pore rotation angle
			TWAFFLE=ATAN2(Y,X)+PIBY4
		ENDIF
C Calculate radius
		RXY=SQRT(X**2+Y**2)
		IF(RXY.GT.0.0) THEN
			PL=PC(7)*PC(2)/RXY
		ENDIF
		PL=MAX(MIN(PL,PC(14)),PC(13))
	ELSEIF(IPACK.EQ.5) THEN
C random packing
		IF(PC(8).GT.0.0) THEN
			X=(INT(PC(8)/PC(3))+0.5)*PC(3)
		ELSE
			X=(INT(PC(8)/PC(3))-0.5)*PC(3)
		ENDIF
		IF(PC(9).GT.0.0) THEN
			Y=(INT(PC(9)/PC(4))+0.5)*PC(4)
		ELSE
			Y=(INT(PC(9)/PC(4))-0.5)*PC(4)
		ENDIF
		RXY=SQRT(X**2+Y**2)
C set NRAD for cartesian packing so that dont loose central pores
		NRAD=1
C generate random rotation
		PIBY4=PI/4.0
		TWAFFLE=SYS_DRAND()*PIBY4
		IF(RXY.GT.0.0) THEN
			PL=PC(7)*PC(2)/RXY
		ENDIF
		PL=MAX(MIN(PL,PC(14)),PC(13))
	ELSEIF(IPACK.EQ.6) THEN
C MIXS packing
C The initial impact point is on the spherical reference surface
C Here we find the annular aperture containing the ray and calculate the
C impact point on the front surface of the MCP
C PC(8) and PC(9) are the local intersection point on spherical surface
C PL is pore length
		CALL SRT_FINDMOD(PC(8),PC(9),IAPER,RMOD,PMOD,TMOD,
     +		WMOD,HMOD,PL,CM,GM,XX,YY)
C RP is the radius of curvature of the spherical surface at this point
		RP=PC(1)
		IF(IAPER.GT.0.AND.PL.NE.PC(14)) THEN
			RPS=PC(1)
			PC(1)=PC(1)-(PC(14)-PL)
			CALL SRT_SPHR(DIR,HPOS,CAX,CAR,CVR,
     +			PC,0,5,HIT,POS,RNM,ISTAT)
			CALL SRT_FINDMOD(PC(8),PC(9),IAPER,RMOD,PMOD,TMOD,
     +			WMOD,HMOD,PL,CM,GM,XX,YY)
			RP=PC(1)
			PC(1)=RPS
		ENDIF
C RPOS is the local radius of this point from the axis
		RPOS=SQRT(PC(8)**2+PC(9)**2)
C Set cartesian and polar coordinates on plate
C XX is mm in radial direction, YY is azimuthal angle in plate sector
		XPLATE=XX
		YPLATE=YY*RPOS
		RPLATE=SQRT(XPLATE**2+YPLATE**2)
		TPLATE=ATAN2(YPLATE,XPLATE)
C Set half width of plate
		WPL=WMOD*0.5
		IF(IAPER.GT.0) THEN
C Now we use the impact position to find rotation angle of the fibre ATHETA
			FIBRE=PC(15)
			NRAD=INT(RPOS/FIBRE)
			RXY=(NRAD+0.5)*FIBRE
			DTHETA=2.0*PI/(INT(2.0*PI*RXY/FIBRE))
			ATHETA=ATAN2(PC(9),PC(8))
			IF(ATHETA.GT.0.0) THEN
				ATHETA=(INT(ATHETA/DTHETA)+0.5)*DTHETA
			ELSE
				ATHETA=(INT(ATHETA/DTHETA)-0.5)*DTHETA
			ENDIF
C Use impact position to find centre of pore and rotation angle
			NRAD=INT(RPOS/PC(3))
			RXY=(NRAD+0.5)*PC(3)
			DTHETA=2.0*PI/(INT(2.0*PI*RXY/PC(4)))
			THETA=ATAN2(PC(9),PC(8))
			IF(THETA.GT.0.0) THEN
				THETA=(INT(THETA/DTHETA)+0.5)*DTHETA
			ELSE
				THETA=(INT(THETA/DTHETA)-0.5)*DTHETA
			ENDIF
			X=RXY*COS(THETA)
			Y=RXY*SIN(THETA)
C Set TWAFFLE as angular difference between pore angle and fibre angle
			TWAFFLE=ATHETA-THETA
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
	ELSEIF(IPACK.EQ.7) THEN
C NFL packing
C The impact point is on the spherical reference surface
C Here we find the rectangular aperture containing the ray
C PC(8) and PC(9) are the local intersection point on spherical surface
C PL is pore length
		CALL SRT_FINDMOD(PC(8),PC(9),IAPER,RMOD,PMOD,TMOD,
     +		WMOD,HMOD,PL,CM,GM,XX,YY)
C Get half width of plate
		if(HMOD.EQ.0.0) THEN
			WPL=WMOD*0.5
		else
			WPL=MIN(WMOD,HMOD)*0.5
		endif
C Set cartesian and polar coordinates on plate
		XPLATE=XX
		YPLATE=YY
		RPLATE=SQRT(XX**2+YY**2)
		TPLATE=ATAN2(YY,XX)
		IF(IAPER.GT.0) THEN
C Cartesian packing
			IF(PC(8).GT.0.0) THEN
				X=(INT(PC(8)/PC(3))+0.5)*PC(3)
			ELSE
				X=(INT(PC(8)/PC(3))-0.5)*PC(3)
			ENDIF
			IF(PC(9).GT.0.0) THEN
				Y=(INT(PC(9)/PC(4))+0.5)*PC(4)
			ELSE
				Y=(INT(PC(9)/PC(4))-0.5)*PC(4)
			ENDIF
			RXY=SQRT(X**2+Y**2)
			PL=MAX(MIN(PL,PC(14)),PC(13))
C Set TWAFFLE for cartesian plate
		        TWAFFLE=0.0
C set NRAD for cartesian packing so that dont loose central pores
			NRAD=1
C Now find position within multifibre structure
			FIBRE=PC(15)
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
	ENDIF
C X and Y are now the position of the pore axis on spherical surface
C Find Y reference axis at local origin on spherical surface
	CALL SRT_VCRS(CAX,CAR,VY)
	CALL SRT_VNRM(VY,ISTAT)
C Find surface normal at pore centre
	Z=SQRT(RP**2-X**2-Y**2)
	ZP(1)=Z*CAX(1)+X*CAR(1)+Y*VY(1)
	ZP(2)=Z*CAX(2)+X*CAR(2)+Y*VY(2)
	ZP(3)=Z*CAX(3)+X*CAR(3)+Y*VY(3)
	CALL SRT_VNRM(ZP,ISTAT)
C Find position of pore centre
	POS(1)=CVR(1)+ZP(1)*RP
	POS(2)=CVR(2)+ZP(2)*RP
	POS(3)=CVR(3)+ZP(3)*RP
	IF(IDEF(1).GT.0) THEN
C Deformations specified
C DX  is effective displacement of centre of sphere in direction CAR
C DY is effective displacement of centre of sphere in direction NORM X CAR
C DA is small change in pore rotation
C DC is small change in pore shape
C DB mean value for Gaussian in-plane figure errors radians
C If IPACK=7 then DX and DY are read as
C     DS scaling of intrinsic slump, 1.0 gives model amplitude
C     DE maximum amplitude of thermoelastic axial pointing errors radians
C     DS and DE then used to calculate DX and DY
C Note that the derivatives are not used
	    IF(IPACK.NE.6.AND.IPACK.NE.7) THEN
C For most packing use the aperture position X,Y
		IDEF(2)=1
		CALL SRT_DFRM(IDEF,X,Y,DX,DDX,DDY,ISTAT)
		IDEF(2)=2
		CALL SRT_DFRM(IDEF,X,Y,DY,DDX,DDY,ISTAT)
		IDEF(2)=3
		CALL SRT_DFRM(IDEF,X,Y,DA,DDX,DDY,ISTAT)
		IDEF(2)=4
		CALL SRT_DFRM(IDEF,X,Y,DC,DDX,DDY,ISTAT)
		IDEF(2)=5
		CALL SRT_DFRM(IDEF,X,Y,DB,DDX,DDY,ISTAT)
	        CALL SYS_GAUSS(2,RAN,0.0D0,1.0D0,ISTAT)
C in-plane figure errors of individual pores
		DB1=RAN(1)*DB
		DB2=RAN(2)*DB
	    ELSEIF(IPACK.EQ.6) THEN
C For IPACK=6 (MIXS-T) use the aperture index
C Combine the pore axis deformation with systematic error
C predicted from the slumping compression/expansion
C DX and DY from deformation are the rms tilt error per multi-fibre
		IDEF(2)=1
		CALL SRT_IDFRM(IDEF,IAPER,1,DX,DDX,DDY,ISTAT)
		IDEF(2)=2
		CALL SRT_IDFRM(IDEF,IAPER,1,DY,DDX,DDY,ISTAT)
C DDD is the reduction factor applied in the corners of the plate
		DDD=1.0
		if(RPLATE.GT.WPL) then
			RMP=WPL*SQRT(2.0)
			DDD=DDD*(RMP-RPLATE)/(RMP-WPL)
		endif
C Pick random numbers to set multi-fibre tilts
	        CALL SYS_GAUSS(2,RAN,0.0D0,1.0D0,ISTAT)
		DX=DX*RAN(1)+COS(TPLATE)*DDD*RPLATE**2/2.0/RP
		DY=DY*RAN(2)+SIN(TPLATE)*DDD*RPLATE**2/2.0/RP
C DA is the pore axis rotation error
		IDEF(2)=3
		CALL SRT_IDFRM(IDEF,IAPER,1,DA,DDX,DDY,ISTAT)
C Combine the shear deformation with the shear expected from the slump
C Also include the likelihood of a shear error within the multifibre
C structure
		IDEF(2)=4
		CALL SRT_IDFRM(IDEF,IAPER,1,DC,DDX,DDY,ISTAT)
		SFIB=PC(3)*4.0
		DC=DC*EXP(-RFIB**2/2.0/SFIB**2)*SIGN(1.0D0,TAN(TFIB))
C *(SYS_DRAND()*0.5+0.5)
		DC=DC+RPLATE*ABS(SIN(2.0*TPLATE))/RP*(PC(4)-PC(6))*
     +		SIGN(1.0D0,TAN(TPLATE))*DDD
		DSHR=DC/(PC(4)-PC(6))
C DB is the pore rms figure error radians
		IDEF(2)=5
		CALL SRT_IDFRM(IDEF,IAPER,1,DB,DDX,DDY,ISTAT)
	        CALL SYS_GAUSS(2,RAN,0.0D0,1.0D0,ISTAT)
		DB1=RAN(1)*DB
		DB2=RAN(2)*DB
	    ELSEIF(IPACK.EQ.7) THEN
C For IPACK=7 (NFL) use the aperture index
C Combine the pore axis deformation with systematic error
C predicted from the slumping compression/expansion
C DS scaling of intrinsic slump, 1.0 gives model amplitude
C DE maximum amplitude of thermoelastic axial pointing errors radians
		IDEF(2)=1
		CALL SRT_IDFRM(IDEF,IAPER,1,DS,DDX,DDY,ISTAT)
		IDEF(2)=2
		CALL SRT_IDFRM(IDEF,IAPER,1,DE,DDX,DDY,ISTAT)
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
C Intrinsic slump tilt error 
		DX=DS*COS(TPLATE)*DDD*RPLATE**2/2.0/RP
		DY=DS*SIN(TPLATE)*DDD*RPLATE**2/2.0/RP
C Set correlation factor
C		WFR=0.5
C		IF(SYS_DRAND().LT.WFR) THEN
C			CRFACT=4.0
C		ELSE
C			CRFACT=1.0
C		ENDIF
		CRFACT=RPLATE/WPL
C Thermoelastic axial pointing errors if DE +ve
		IF(DE.GE.0.0) THEN
			DX=DX+DE*COS(TPLATE)*QQQ*RP
			DY=DY+DE*SIN(TPLATE)*QQQ*RP
		ELSE
C If negative then use as rms pointing error in radians on just
C one axis chosen at random
C	        	CALL SYS_GAUSS(2,RAN,0.0D0,1.0D0,ISTAT)
C			IF(RAN(1).GT.0.0) THEN
C				DX=DX+DE*RAN(2)*Z*CRFACT
C			ELSE
C				DY=DY+DE*RAN(2)*Z*CRFACT
C			ENDIF
C Tilt errors with minima on diagonals, maxima at middle of edges
C			TTT=-DE*Z*(RPLATE/WPL)**2*(COS(TPLATE*2.0)**2-0.5)
C			DX=DX+TTT*COS(TPLATE)
C			DY=DY+TTT*SIN(TPLATE)
C Chess board pattern of tilt errors
			TTT=2.0*PI/15.0
			PPP=SIN(TTT*XPLATE+TTT*YPLATE/8.0)
     +			*SIN(TTT*YPLATE+TTT*XPLATE/8.0)
			PPP=SIGN(ABS(PPP)**2.0,PPP)
			DX=DX+DE*PPP*Z*SIN(RPLATE*TTT)
			DY=DY+DE*PPP*Z*COS(RPLATE*TTT)
		ENDIF
C DA is the pore axis rotation error
		IDEF(2)=3
		CALL SRT_IDFRM(IDEF,IAPER,1,DA,DDX,DDY,ISTAT)
C Fix pattern twist 
C		DA=DA*RPLATE/WPL
C Combine the shear deformation with the shear expected from the slump
C Also include the likelihood of a shear error within the multifibre
C structure
		IDEF(2)=4
		CALL SRT_IDFRM(IDEF,IAPER,1,DC,DDX,DDY,ISTAT)
		SFIB=PC(3)*3.0
C Shear errors within the multifibre structure
		DC=DC*EXP(-RFIB**2/2.0/SFIB**2)*SIGN(1.0D0,TAN(TFIB))
C Shear errors intrinsic to the slump
		DC=DC+DS*RPLATE*ABS(SIN(2.0*TPLATE))/RP*(PC(4)-PC(6))*
     +		SIGN(1.0D0,TAN(TPLATE))*DDD
		DSHR=DC/(PC(4)-PC(6))
C DB is the pore rms figure error radians
C Increase rms figure errors for a fraction of the pores
C		BAND=5.0
C		IF(ABS(XFIB).LT.PC(3)*BAND.OR.ABS(YFIB).LT.PC(3)*BAND) THEN
C			CRFACT=4.0
C		ELSE
C			CRFACT=1.0
C		ENDIF
		IDEF(2)=5
		CALL SRT_IDFRM(IDEF,IAPER,1,DB,DDX,DDY,ISTAT)
		IF(DB.GT.0.0) THEN
	                CALL SYS_GAUSS(3,RAN,0.0D0,1.0D0,ISTAT)
			DB1=RAN(1)*DB*CRFACT
			DB2=RAN(2)*DB*CRFACT
		ELSE
C Use Cauchy PDF to set tilt errors - DB is the scaling (2*DB=FWHM)
C                        DB1=DB*TAN((SYS_DRAND()-0.5)*PI)
C                        DB2=DB*TAN((SYS_DRAND()-0.5)*PI)
C Use modified Lorentzian (King) profile. When alpha=1 then Lorentzian
C and FWHM=2*DB
   		        ALPHA=1.4
                        DB1=DB*SRT_KINGRAND(ALPHA)
                        DB2=DB*SRT_KINGRAND(ALPHA)
C			DB1=0.0
C			DB2=0.0
C			DX=DX*(1.0+DB*RAN(1))
C			DY=DY*(1.0+DB*RAN(2))
C			DSHR=DSHR*(1.0+DB*RAN(3))
		ENDIF
	    ENDIF
C Apply deformations to change effective centre of sphere
	    ZP(1)=Z*CAX(1)+(X-DX+DB1*Z)*CAR(1)+(Y-DY+DB2*Z)*VY(1)
	    ZP(2)=Z*CAX(2)+(X-DX+DB1*Z)*CAR(2)+(Y-DY+DB2*Z)*VY(2)
	    ZP(3)=Z*CAX(3)+(X-DX+DB1*Z)*CAR(3)+(Y-DY+DB2*Z)*VY(3)
	    CALL SRT_VNRM(ZP,ISTAT)
C Change rotation of pore
	    TWAFFLE=TWAFFLE+DA
C DC is pore shear (used in setting up pore surfaces)
	ELSE
	    DSHR=0.0
	    DC=0.0
	ENDIF
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
	IF(IPACK.EQ.1.OR.IPACK.EQ.7) THEN
C For cartesian packing set the pore rotation angle using the reference axes
		CALL SRT_VDOT(XP,CAR,CTH)
		CALL SRT_VDOT(XP,VY,STH)
		TWAFFLE=TWAFFLE+ATAN2(STH,CTH)
	ENDIF
C rotate pore to correct angle
	DO I=1,3
		XP(I)=XP(I)*COS(TWAFFLE)+YP(I)*SIN(TWAFFLE)
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
		HX=(PC(3)-PC(5))*0.5
		HY=(PC(4)-PC(6))*0.5
		IF(IPACK.EQ.4.OR.IPACK.EQ.5) THEN
			HX=HX/SQRT(2.0)
			HY=HY/SQRT(2.0)
		ENDIF
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
C Get surface quality of pore walls
	IQ=INT(PC(11))
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
