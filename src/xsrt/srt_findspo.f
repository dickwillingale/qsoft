*+SRT_FINDSPO    Find SPO module aperture in array
	SUBROUTINE SRT_FINDSPO(XR,YR,MHIT,RM,PM,TM,WM,HM,LM,CU,GM,
     +  PP,WA,RR,WR,FW,SQ,SP,XX,YY)
	IMPLICIT NONE
      	INTEGER MHIT
	DOUBLE PRECISION XR,YR,RM,PM,TM,WM,HM,LM,XX,YY
	DOUBLE PRECISION CU,GM,PP,WA,RR,WR,FW,SQ,SP
*XR,YR	input	ray position on aperture plane
*MHIT	output	index of module (0 if missed)
*RM	output	radius of module
*PM	output	azimuth of module
*TM	output	rotation of module wrt radius vector
*WM	output	width of module (radial)
*HM	output	height of module (azimuthal)
*LM	output	axial length of module
*CU	output	axial curvature index of module
*GM	output	grazing angle ratio of module
*PP	output	radial pitch of pores (pore size + wall)
*WA	output	membrane thickness between rows of pores
*RR	output	azimuthal pitch of pores (rib spacing)
*WR	output	width of ribs along rows of pores
*FW	output	module frame width
*SQ	output	module reflecting surface quality index
*SP	outout	spare parameter
*XX,YY	output	local coordinates of ray position
*-Author Dick Willingale 2018-Mar-20
	INCLUDE 'SRT_COM'
	DOUBLE PRECISION XC,YC,PHI,WID,HIG,PHIS,THETA,RP
	DOUBLE PRECISION X,Y,ELW,ELH,RPHI,RRAD,RXY
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
C Assume modules have rectangular aperture
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
C test to see if ray inside module aperture
		IF(ABS(XX).LE.ELW) THEN
			IF(ABS(YY).LE.ELH) THEN
C If in module aperture then pull out module parameters
				MHIT=ITRY
				RM=rlist(ITRY)
				PM=philist(ITRY)
				TM=thlist(ITRY)
				WM=wlist(ITRY)
				HM=hlist(ITRY)
				LM=llist(ITRY)
			    	CU=clist(ITRY)
				GM=glist(ITRY)
				PP=olist(ITRY)
				WA=plist(ITRY)
				RR=qlist(ITRY)
				WR=ulist(ITRY)
				FW=vlist(ITRY)
				SQ=zlist(ITRY)
				SP=slist(ITRY)
			ENDIF
		ENDIF
		ITRY=ITRY+1
	ENDDO
	END
