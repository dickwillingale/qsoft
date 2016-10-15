*+SRT_SFINDMOD    Find module aperture on spherical surface
	SUBROUTINE SRT_SFINDMOD(SRAD,XR,YR,MHIT,RM,PM,TM,WM,HM,LM,CM,GM,XX,YY)
	IMPLICIT NONE
      	INTEGER MHIT
	DOUBLE PRECISION SRAD,XR,YR,RM,PM,TM,WM,HM,LM,CM,GM,XX,YY
*SRAD	input	radius of curvature of spherical surface
*XR,YR	input	ray position on aperture plane
*MHIT	output	index of module (0 if missed)
*RM	output	radius of module
*PM	output	azimuth of module
*TM	output	rotation of module wrt radius vector
*WM	output	width of module x (radial)
*HM	output	height of module y (azimuthal)
*LM	output	length of module (axial)
*CM	output	curvature signature of module
*GM	output	grazing angle ratio of module
*XX,YY	output	local coordinates of ray position
*-Author Dick Willingale 2012-Jul-6
	INCLUDE 'SRT_COM'
	DOUBLE PRECISION XC,YC,PHI,WID,HIG,PHIS,THETA,RP
	DOUBLE PRECISION X,Y,ELW,ELH,RPHI,RRAD
	INTEGER ITRY
C Find module
	ITRY=1
	MHIT=0
	DO WHILE(MHIT.EQ.0.AND.ITRY.LE.nrofelts)
C Position of module centre
		PHI=philist(ITRY)
		RP=rlist(ITRY)
C Dimensions of module aperture
		ELW=wlist(ITRY)*0.5
		ELH=hlist(ITRY)*0.5
		IF(iptype.EQ.1) THEN
C Rectangular modules (width x radial, height y azimuthal)
			XC=RP*COS(PHI)
			YC=RP*SIN(PHI)
C Module rotation wrt to full aperture axes
			THETA=PHI+thlist(ITRY)
C Position of ray wrt centre of module
			X=XR-XC
			Y=YR-YC
C Transform position to module coordinates
			XX=X*cos(THETA)+Y*sin(THETA)
			YY=-X*sin(THETA)+Y*cos(THETA)
		ELSE
C Sector modules (width x radius, height y azimuth)
			RRAD=SQRT(XR**2+YR**2)
			RPHI=ATAN2(YR,XR)
			XX=RRAD-RP
			YY=RPHI-PHI
		ENDIF
C test to see if ray inside module aperture
		IF(ABS(XX).LE.ELW) THEN
			IF(ABS(YY).LE.ELH) THEN
				MHIT=ITRY
				RM=rlist(ITRY)
				PM=philist(ITRY)
				TM=thlist(ITRY)
				WM=wlist(ITRY)
				HM=hlist(ITRY)
				LM=llist(ITRY)
				CM=clist(ITRY)
				GM=glist(ITRY)
			ENDIF
		ENDIF
		ITRY=ITRY+1
	ENDDO
	END
