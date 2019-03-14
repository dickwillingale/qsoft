*+SRT_FINDMPO    Find MPO aperture in array
	SUBROUTINE SRT_FINDMPO(XR,YR,MHIT,XM,YM,TM,WM,HM,LM,CU,MF,
     +  PP,WA,SQ,BU,BV,BZ,SP,XX,YY)
	IMPLICIT NONE
      	INTEGER MHIT
	DOUBLE PRECISION XR,YR,XM,YM,TM,WM,HM,LM,XX,YY
	DOUBLE PRECISION CU,MF,PP,WA,SQ,BU,BV,BZ,SP
*XR,YR	input	ray position on aperture plane
*MHIT	output	index of MPO (0 if missed)
*XM	output	x position of MPO
*YM	output	y position of MPO
*TM	output	rotation of MPO wrt radius vector
*WM	output	width of MPO delx
*HM	output	height of MPO dely
*LM	output	axial length of MPO (thickness)
*CU	output	radius of curvature of MPO
*MF	output	MPO multifibre size 
*PP	output	MPO pitch of pores (pore size + wall)
*WA	output	MPO wall thickness between pores
*SQ	output	MPO reflecting surface quality index
*BU	output	MPO bias angle x radians
*BV	output	MPO bias angle y radians
*BZ	output	MPO efficiency wrt to theory for 1 reflection in pore
*SP	outout	spare MPO parameter
*XX,YY	output	local coordinates of ray position
*-Author Dick Willingale 2017-Oct-23
	INCLUDE 'SRT_COM'
	DOUBLE PRECISION XC,YC,PHI,WID,HIG,PHIS,THETA,RP
	DOUBLE PRECISION X,Y,ELW,ELH,RPHI,RRAD,RXY
	INTEGER ITRY
C zero all parameters
	XM=0.0
	YM=0.0
	TM=0.0
	WM=0.0
	HM=0.0
	LM=0.0
	CU=0.0
	MF=0.0
	PP=0.0
	WA=0.0
	SQ=0.0
	BU=0.0
	BV=0.0
	BZ=0.0
	SP=0.0
C Find module
	ITRY=1
	MHIT=0
	DO WHILE(MHIT.EQ.0.AND.ITRY.LE.nrofelts)
C Position of module centre
		XC=rlist(ITRY)
		YC=philist(ITRY)
C Dimensions of module aperture
		ELW=wlist(ITRY)*0.5
		ELH=hlist(ITRY)*0.5
C Module rotation wrt to full aperture axes
		THETA=thlist(ITRY)
C Position of ray wrt centre of module
		X=XR-XC
		Y=YR-YC
C Transform position to module coordinates
		XX=X*cos(THETA)+Y*sin(THETA)
		YY=-X*sin(THETA)+Y*cos(THETA)
C test to see if ray inside module aperture
		IF(ABS(XX).LE.ELW) THEN
			IF(ABS(YY).LE.ELH) THEN
C If in MPO aperture then pull out MPO parameters
				MHIT=ITRY
				XM=rlist(ITRY)
				YM=philist(ITRY)
				TM=thlist(ITRY)
				WM=wlist(ITRY)
				HM=hlist(ITRY)
				LM=llist(ITRY)
			    	CU=clist(ITRY)
				MF=glist(ITRY)
				PP=olist(ITRY)
				WA=plist(ITRY)
				SQ=qlist(ITRY)
				BU=ulist(ITRY)
				BV=vlist(ITRY)
				BZ=zlist(ITRY)
				SP=slist(ITRY)
			ENDIF
		ENDIF
		ITRY=ITRY+1
	ENDDO
	END
