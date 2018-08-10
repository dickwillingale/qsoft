*+SRT_REFS	Reflection, refraction, scattering and diffraction
	SUBROUTINE SRT_REFS(POS,SP,DIR,DNM,ISU,DRF,QRY,ISTAT)
	IMPLICIT NONE
	DOUBLE PRECISION POS(3),SP(7),DIR(3),DNM(3),DRF(3),QRY(3)
	INTEGER ISU,ISTAT
*POS	input	position on surface
*SP	input	surface parameters
*DIR	input	direction of incident ray
*DNM	input	normal to surface
*ISU	input	index to surface quality table
*DRF	output	reflected, refracted, scattered direction
*QRY	in/out	quality factors for ray
*ISTAT	in/out	returned status
* QRY(1)	wavelength A
* QRY(2)	unpolarized area
*-Author Dick Willingale 1996-Nov-21
* Modified for ITYPE=4 grating 2005-Jul-13
	INCLUDE 'SRT_COM'
C Surface quality parameters are found using common ISQP(3,MAXST)
C type of parameter list 	ITYPE=ISQP(1,ISU)
C parameter index 		IP=ISQP(2,ISU)
C number of parameters		NP=ISQP(3,ISU)
C If ITYPE=1, X-ray reflection and scattering using refractive index
C 	P(IP) 	wavelength A
C 	P(IP+1)	specific roughness A**2 mm
C		if -ve then rms figure gradient error radius
C	P(IP+2)	minimum surface frequency for roughness mm-1
C	P(IP+3)	surface roughness power law index
C	P(IP+4)	unused
C	P(IP+5)	unused
C	P(IP+6)	unused
C 	P(IP+7)	alpha, real part of dielectric constant
C 	P(IP+8)	gamma, imaginary part of dielectric constant
C If ITYPE=2, X-ray reflection and scattering using a look-up table
C	P(IP)	wavelength A
C 	P(IP+1)	specific roughness A**2 mm
C		if -ve then rms figure gradient error radius
C	P(IP+2)	minimum surface frequency for roughness mm-1
C	P(IP+3)	surface roughness power law index
C	P(IP+4)	unused
C	P(IP+5)	unused
C	P(IP+6)	unused
C	P(IP+7) unused
C	P(IP+8) unused
C 	P(IP+9)	start of incidence angle and reflectivity sample pairs
C	P(IP+NP-1)	end of sample pairs, angle samples increasing 
C 	If all parameters 0 then perfect reflecting surface
C If ITYPE=3, refraction using refractive index
C 	P(IP) 	wavelength A
C 	P(IP+1)	specific roughness A**2 mm
C		if -ve then rms figure gradient error radius
C	P(IP+2)	minimum surface frequency for roughness mm-1
C	P(IP+3)	surface roughness power law index
C	P(IP+4)	unused
C	P(IP+5)	unused
C	P(IP+6)	unused
C 	P(IP+7)	refractive index n2 on side of -ve normal or n2/n1
C	P(IP+8) unused
C If ITYPE=4, diffraction grating
C	P(IP)	wavelength A
C 	P(IP+1)	specific roughness A**2 mm
C		if -ve then rms figure gradient error radius
C	P(IP+2)	minimum surface frequency for roughness mm-1
C	P(IP+3)	surface roughness power law index
C	P(IP+4) d spacing of grating mm
C	P(IP+5)	spacing gradient along ruling (hub distance off-plane)
C               if <1.0 then spacing gradient across ruling in-plane
C		ruling direction specified by surface axis
C	P(IP+6)	diffraction order
C	P(IP+7) unused
C	P(IP+8) unused
C 	P(IP+9)	start of incidence angle and reflectivity sample pairs
C	P(IP+NP-1)	end of sample pairs, angle samples increasing 
C 	If all parameters 0 then perfect reflecting surface
	DOUBLE PRECISION SYS_DRAND
	DOUBLE PRECISION COSA,ANGL,GRAZ,WL,SPECROUGH,RINDROUGH,BREAKROUGH
	DOUBLE PRECISION GAMMA0,GAMMA1,RND,OMEGA,RR,ADEL,AB,COSB,COSAB
	DOUBLE PRECISION RSI,RPI,REF,SINA,SINB,SR
	DOUBLE PRECISION RULE(3),YAX(3),COSG,SING,GAMMA,BETA,ALPHA
	DOUBLE PRECISION VI(3),VD(3)
	DOUBLE PRECISION DSPACE,DHUB,DBYR,DBYA,VP(3),ORDER
	DOUBLE PRECISION XX,YY
	DOUBLE PRECISION HUB(3),TIS
	INTEGER IP,NP,J,IABS
C
	IF(ISTAT.NE.0) RETURN
C
	IF(ISU.LT.1.OR.ISU.GT.MAXST) THEN
		WRITE(*,*) 'SRT_REFS error - quality index out of range'
		ISTAT=1
		RETURN
	ENDIF
	IP=ISQP(2,ISU)
	IF(IP.EQ.0.OR.IP.GT.MAXPAR) THEN
		WRITE(*,*) 'SRT_REFS error - parameter index out of range'
		ISTAT=1
		RETURN
	ENDIF
C
	IABS=0
C
	WL=PAR(IP)
C save wavelength in ray quality (1)
	QRY(1)=WL
	SPECROUGH=PAR(IP+1)
	BREAKROUGH=PAR(IP+2)
	RINDROUGH=PAR(IP+3)
C Find incidence angle 
	CALL SRT_VDOT(DIR,DNM,COSA)
	ANGL=ACOS(COSA)
C Calculate grazing angle
	IF(ANGL.GT.PIBY2) THEN
		GRAZ=ANGL-PIBY2
	ELSE
		GRAZ=PIBY2-ANGL
	ENDIF
	IF(ISQP(1,ISU).EQ.1.OR.ISQP(1,ISU).EQ.2) THEN
C Perform reflection
		IF(ANGL.LT.0.0.OR.ANGL.GT.PIBY2) THEN
			COSB=COSA
		ELSE
			COSB=COS(ANGL)
		ENDIF
		COSAB=COSA+COSB
		DRF(1)=DIR(1)-COSAB*DNM(1)
		DRF(2)=DIR(2)-COSAB*DNM(2)
		DRF(3)=DIR(3)-COSAB*DNM(3)
		CALL SRT_VNRM(DRF,ISTAT)
C Calculate scattering angle
	        CALL SRT_SCAT(WL,SPECROUGH,BREAKROUGH,RINDROUGH,
     +		GRAZ,1,ADEL,ISTAT)
		IF(ADEL.NE.0.0) THEN
C Limit scattering range - absorb if into surface
			IF(ADEL.LT.-GRAZ*0.99) THEN
				ADEL= -GRAZ*0.99
				IABS=1
			ENDIF
			IF(ADEL.GT.(PI-GRAZ)*0.99) THEN
				ADEL= (PI-GRAZ)*0.99
				IABS=1
			ENDIF
C find vector perpendicular to reflection plane
			CALL SRT_VCRS(DRF,DNM,YAX)
C find vector in reflection plane perpendicular to reflected ray
			CALL SRT_VCRS(YAX,DRF,VP)
			CALL SRT_VNRM(VP,ISTAT)
C Scatter ray
			ADEL=TAN(ADEL)
			DRF(1)=DRF(1)+VP(1)*ADEL
			DRF(2)=DRF(2)+VP(2)*ADEL
			DRF(3)=DRF(3)+VP(3)*ADEL
		ENDIF
		CALL SRT_VNRM(DRF,ISTAT)
	ELSEIF(ISQP(1,ISU).EQ.3) THEN
C Perform refraction, scattering by roughness not included
		SINA=SIN(ANGL)
		IF(COSA.LT.0.0) THEN
			SINB=SINA/PAR(IP+6)
		ELSE
			SINB=SINA*PAR(IP+6)
		ENDIF
		IF(ABS(SINB).LT.1.0) THEN
C refraction
			COSB=SQRT(1.D0-SINB**2)
			IF(SINA.NE.0.0) THEN
				SR=SINB/SINA
			ELSE
				SR=1.0
			ENDIF
			IF(COSA.LT.0.0) THEN
				DRF(1)=(DIR(1)+DNM(1))*SR-DNM(1)*COSB
				DRF(2)=(DIR(2)+DNM(2))*SR-DNM(2)*COSB
				DRF(3)=(DIR(3)+DNM(3))*SR-DNM(3)*COSB
			ELSE
				DRF(1)=(DIR(1)-DNM(1))*SR+DNM(1)*COSB
				DRF(2)=(DIR(2)-DNM(2))*SR+DNM(2)*COSB
				DRF(3)=(DIR(3)-DNM(3))*SR+DNM(3)*COSB
			ENDIF
		ELSE
C Total internal reflection
			COSAB=COSA+COSA
			DRF(1)=DIR(1)-DNM(1)*COSAB
			DRF(2)=DIR(2)-DNM(2)*COSAB
			DRF(3)=DIR(3)-DNM(3)*COSAB
		ENDIF
		CALL SRT_VNRM(DRF,ISTAT)
	ELSEIF(ISQP(1,ISU).EQ.4) THEN
C Diffraction grating
		DSPACE=PAR(IP+4)
		DHUB=PAR(IP+5)
		ORDER=PAR(IP+6)
C Set ruling direction and d-spacing gradients
		IF(DHUB.GT.1.0) THEN
C find hub position. Note ruling direction from hub to grating
			DO J=1,3
				HUB(J)=SP(6+J)-SP(3+J)*DHUB
			ENDDO
                	CALL SRT_DIDI(HUB,POS,RULE,RR)
			CALL SRT_VNRM(RULE,ISTAT)
			DBYR=DSPACE/DHUB
			DBYA=0.0
		ELSE
			DO J=1,3
				RULE(J)=SP(3+J)
			ENDDO
			DBYR=0.0
			DBYA=DHUB
		ENDIF
C Find other reference axis on surface
		CALL SRT_VCRS(DNM,RULE,YAX)
		CALL SRT_VNRM(YAX,ISTAT)
C Find vector from origin on plane to intersection
                CALL SRT_DIDI(SP(7),POS,VP,RR)
C Calculate local x and y on plane (x is ruling direction)
                CALL SRT_VDOT(VP,SP(4),XX)
                XX=XX*RR
                CALL SRT_VDOT(VP,YAX,YY)
		YY=YY*RR
		DSPACE=DSPACE+DBYR*XX+DBYA*YY
C find angle between ruling direction and ray
		CALL SRT_VDOT(DIR,RULE,COSG)
		GAMMA=ACOS(COSG)
		SING=SIN(GAMMA)
C find incident cone vector
		DO J=1,3
			VI(J)=RULE(J)*COSG-DIR(J)
		ENDDO
		CALL SRT_VNRM(VI,ISTAT)
C find angle between normal and incident cone vector
		CALL SRT_VDOT(DNM,VI,COSA)
		ALPHA=ACOS(COSA)
C use grating equation to find angle between normal and diffracted cone vector
		SINB=ORDER*WL*1.0D-7/(DSPACE*SING)-SIN(ALPHA)
		IF(SINB.GT.1.0.OR.SINB.LT.-1.0) THEN
			IABS=1
		ELSE
			IABS=0
		ENDIF
		SINB=MIN(MAX(SINB,-1.0),1.0)
		BETA=ASIN(SINB)
		COSB=COS(BETA)
C	write(*,*) 'grating',xx,yy,dspace,dbyr*xx,dbya*yy,order
C	write(*,*) 'nominal ruling',sp(4),sp(5),sp(6)
C	write(*,*) 'hub ruling',rule
C	write(*,*) 'gamma alpha beta',gamma*180./pi,alpha*180./pi,beta*180./pi
C	write(*,*) 'grazing angle',asin(sin(gamma)*cos(alpha))*180.0/pi
C find diffracted cone vector
		DO J=1,3
			VD(J)=-YAX(J)*SINB+DNM(J)*COSB
		ENDDO
		CALL SRT_VNRM(VD,ISTAT)
C find diffracted direction
		DO J=1,3
			DRF(J)=RULE(J)*COSG+VD(J)*SING
		ENDDO
		CALL SRT_VNRM(DRF,ISTAT)
	ENDIF
	IF(ISQP(1,ISU).EQ.1.AND.PAR(IP+7).NE.0.0) THEN
C Use Fresnels equations to calculate reflectivities
		IF(ANGL.GT.PIBY2) THEN
			GRAZ=ANGL-PIBY2
		ELSE
			GRAZ=PIBY2-ANGL
		ENDIF
		CALL SRT_FRNL(GRAZ,PAR(IP+7),PAR(IP+8),RSI,RPI)
C		write(*,*) graz,par(ip+7),par(ip+8),rsi,rpi
		REF=(RSI+RPI)*0.5
		QRY(2)=QRY(2)*REF
	ELSEIF(ISQP(1,ISU).EQ.2.OR.ISQP(1,ISU).EQ.4) THEN
C Use look-up table to estimate reflectivity
		NP=(ISQP(3,ISU)-9)/2
		IF(ANGL.GT.PIBY2) THEN
			ANGL=PIBY2*2-ANGL
		ENDIF
C New version of table has angles in degrees
		ANGL=ANGL*180.0/PI
		CALL SRT_RELK(ANGL,NP,PAR(IP+9),REF)
		QRY(2)=QRY(2)*REF
	ENDIF
	IF(IABS.EQ.1) THEN
		QRY(2)=0.0
	ENDIF
	END
