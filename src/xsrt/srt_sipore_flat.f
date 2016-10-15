*+SRT_SIPORE     Dynamically set up a pair of square Si pores for HPO
	SUBROUTINE SRT_SIPORE(DIR,HPOS,CVR,CAX,CAR,PC,IDEF,ISTAT)
	IMPLICIT NONE
	DOUBLE PRECISION DIR(3),HPOS(3),CAX(3),CAR(3),CVR(3),PC(50)
	INTEGER IDEF(2),ISTAT
*DIR	input	direction of ray
*HPOS	input	position of hit on plane reference surface
*CVR    input   centre of aperture (in front of principal plane)
*CAX    input   axis of optic
*CAR    input   reference direction in aperture
*PC     input   remaining parameters
*IDEF   input   deformation indices
*ISTAT  in/out  returned status
* Parameters for pores are:
* PC(1) inner radius of aperture
* PC(2) outer radius of aperture
* PC(3) focal length
* PC(4) pitch of pores
* PC(5) wall thickness of pores
* PC(6) minimum pore length
* PC(7) maximum pore length
* PC(8) local x of ray intersection
* PC(9) local y of ray intersection
* PC(10) this surface index
* PC(11) surface quality
*-Author Dick Willingale 2008-Oct-2
	DOUBLE PRECISION PI,SMALL
	PARAMETER (PI=3.14159265359,SMALL=1.D-4)
	DOUBLE PRECISION VY(3),X,Y,Z,XP(3),YP(3),ZP(3),RXY,CTH,STH,T,POS(3)
	DOUBLE PRECISION RNM(3),RPOS,CVRN(3)
	DOUBLE PRECISION PP(13),HX,HY,HL,THETA,DTHETA,DX,DY,DT,DDX,DDY,THETAG
	DOUBLE PRECISION PWAFFLE,X0,XX,SS,Y0,TWAFFLE,FL,RPS,RLO,RHI,GAP,FIBRE
	DOUBLE PRECISION ATHETA,T1,T2,T3,T4,PIBY4,ABIT,YY,S45,PL,TGAP,RP
	INTEGER IX,IY,I,KSUR,IQ,IPACK,NRAD,NWAFFLE,IP
	INTEGER NANUL,NSEC
	LOGICAL HIT
C
C Get focal length
	FL=PC(3)
C Find radial position
	RXY=SQRT(PC(8)**2+PC(9)**2)
C Calculate initial guess at grazing angle and length of pore
	THETAG=ATAN2(RXY,FL)/4.0
	PL=MIN(MAX(PC(6),(PC(4)-PC(5))/THETAG),PC(7))
C Move along ray to aperture plane near top of pore
	DO I=1,3
		CVRN(I)=CVR(I)-CAX(I)*(PC(7)-PL-PC(4))
	ENDDO
	CALL SRT_PLNA(DIR,HPOS,CVRN,CAX,CAR,PC,0,9,HIT,POS,RNM,ISTAT)
C Find new radius
	RXY=SQRT(PC(8)**2+PC(9)**2)
C Radial packing of pores - find centre of pore for ray
	NRAD=INT(RXY/PC(4))
	RXY=(NRAD+0.5)*PC(4)
	DTHETA=2.0*PI/(INT(2.0*PI*RXY/PC(5)))
	THETA=ATAN2(PC(9),PC(8))
	IF(THETA.GT.0.0) THEN
		THETA=(INT(THETA/DTHETA)+0.5)*DTHETA
	ELSE
		THETA=(INT(THETA/DTHETA)-0.5)*DTHETA
	ENDIF
	X=RXY*COS(THETA)
	Y=RXY*SIN(THETA)
C check pore position is inside radial aperture
	IF(RXY.LT.PC(1).OR.RXY.GT.PC(2)) THEN
C Pore centre outside aperture so move pore so that ray misses
		X=X+100.0
	ENDIF
C The end of the pore on principal plane is at a slightly smaller
C radius that the entrance to the pore
	X=X*(1.0-ABS(PL*TAN(THETAG)*COS(THETA)/X))
	Y=Y*(1.0-ABS(PL*TAN(THETAG)*SIN(THETA)/Y))
	RXY=SQRT(X**2+Y**2)
C Re-calculate grazing angle and axis length of pore 
C Note that effective radius is a little larger because the reflecting wall
C of the pore is offset from the centre line.
	THETAG=ATAN2(RXY+PC(4)-PC(5),FL)/4.0
	PL=MIN(MAX(PC(6),(PC(4)-PC(5))/TAN(THETAG)),PC(7))
C Find distance to vertex of 1st surface from principal plane
	Z=RXY/TAN(THETAG)
C Find vertex of 1st surface
	DO I=1,3
		CVRN(I)=CVR(I)-CAX(I)*(Z+PC(7))
	ENDDO
C Find Y reference axis at local origin on aperture
	CALL SRT_VCRS(CAX,CAR,VY)
	CALL SRT_VNRM(VY,ISTAT)
C Find direction of pore axis
	ZP(1)=Z*CAX(1)+X*CAR(1)+Y*VY(1)
	ZP(2)=Z*CAX(2)+X*CAR(2)+Y*VY(2)
	ZP(3)=Z*CAX(3)+X*CAR(3)+Y*VY(3)
	CALL SRT_VNRM(ZP,ISTAT)
C Find position of pore centre on principal plane
	RP=SQRT(Z**2+RXY**2)
	POS(1)=CVRN(1)+ZP(1)*RP
	POS(2)=CVRN(2)+ZP(2)*RP
	POS(3)=CVRN(3)+ZP(3)*RP
	IF(IDEF(1).GT.0) THEN
C Get deformations (note deformation gradients are not used)
C DX is effective displacement of centre of sphere in direction CAR
C DY is effective displacement of centre of sphere in direction NORM X CAR
		IDEF(2)=1
		CALL SRT_DFRM(IDEF,X,Y,DX,DDX,DDY,ISTAT)
		IDEF(2)=2
		CALL SRT_DFRM(IDEF,X,Y,DY,DDX,DDY,ISTAT)
C Apply deformation to change effective centre of sphere
		ZP(1)=Z*CAX(1)+(X-DX)*CAR(1)+(Y-DY)*VY(1)
		ZP(2)=Z*CAX(2)+(X-DX)*CAR(2)+(Y-DY)*VY(2)
		ZP(3)=Z*CAX(3)+(X-DX)*CAR(3)+(Y-DY)*VY(3)
		CALL SRT_VNRM(ZP,ISTAT)
	ENDIF
C Now set up axes tangential to surface at pore centre
C Note that there is no pore centred at the origin so can use the cross product
	CALL SRT_VCRS(CAX,ZP,XP)
	CALL SRT_VNRM(XP,ISTAT)
	CALL SRT_VCRS(ZP,XP,YP)
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
C Set dimensions of pore aperture
	IF(NRAD.GT.0) THEN
		HX=(PC(4)-PC(5))*0.5
		HY=(PC(4)-PC(5))*0.5
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
C Find surface index of pore aperture and set that surface element
	KSUR=PC(10)+1
	CALL SRT_SETF(KSUR,1,13,PP,0,0,0,KSUR+1,ISTAT)
C Get surface quality of pore walls
	IQ=INT(PC(11))
C Set sides of pore - local x is along pore axis
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
	PP(1)=YP(1)
	PP(2)=YP(2)
	PP(3)=YP(3)
	PP(4)=ZP(1)
	PP(5)=ZP(2)
	PP(6)=ZP(3)
	PP(7)=POS(1)+ZP(1)*HL-YP(1)*HY
	PP(8)=POS(2)+ZP(2)*HL-YP(2)*HY
	PP(9)=POS(3)+ZP(3)*HL-YP(3)*HY
	PP(10)=-HL
        PP(11)=-HX
        PP(12)=HL
        PP(13)=HX
	CALL SRT_SETF(KSUR+2,5,13,PP,0,IQ,KSUR+1,-1,ISTAT)
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
	CALL SRT_SETF(KSUR+3,5,13,PP,0,0,KSUR+1,-1,ISTAT)
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
	CALL SRT_SETF(KSUR+4,5,13,PP,0,0,KSUR+1,-1,ISTAT)
C Find distance for vertex of 2nd surface from Principal plane
	Z=RXY/TAN(THETAG*3.0)
C Find position of vertex of 2nd surface
	DO I=1,3
		CVRN(I)=CVR(I)-CAX(I)*(Z+PC(7))
	ENDDO
C Find direction of 2nd pore axis
	ZP(1)=Z*CAX(1)+X*CAR(1)+Y*VY(1)
	ZP(2)=Z*CAX(2)+X*CAR(2)+Y*VY(2)
	ZP(3)=Z*CAX(3)+X*CAR(3)+Y*VY(3)
	CALL SRT_VNRM(ZP,ISTAT)
C Find position of 2nd pore centre at principal plane
C Shift to give a small gap between 1st and 2nd pore
	RP=SQRT(Z**2+RXY**2)
	POS(1)=CVRN(1)+ZP(1)*(RP-PC(5))
	POS(2)=CVRN(2)+ZP(2)*(RP-PC(5))
	POS(3)=CVRN(3)+ZP(3)*(RP-PC(5))
	IF(IDEF(1).GT.0) THEN
C Get deformations (note deformation gradients are not used)
C DX is effective displacement of centre of sphere in direction CAR
C DY is effective displacement of centre of sphere in direction NORM X CAR
		IDEF(2)=1
		CALL SRT_DFRM(IDEF,X,Y,DX,DDX,DDY,ISTAT)
		IDEF(2)=2
		CALL SRT_DFRM(IDEF,X,Y,DY,DDX,DDY,ISTAT)
C Apply deformation to change effective centre of sphere
		ZP(1)=Z*CAX(1)+(X-DX)*CAR(1)+(Y-DY)*VY(1)
		ZP(2)=Z*CAX(2)+(X-DX)*CAR(2)+(Y-DY)*VY(2)
		ZP(3)=Z*CAX(3)+(X-DX)*CAR(3)+(Y-DY)*VY(3)
		CALL SRT_VNRM(ZP,ISTAT)
	ENDIF
C Now set up axes tangential to surface at pore centre
C Note that there is no pore centred at the origin so can use the cross product
	CALL SRT_VCRS(CAX,ZP,XP)
	CALL SRT_VNRM(XP,ISTAT)
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
        PP(10)=-HX
        PP(11)=-HY
        PP(12)=HX
        PP(13)=HY
C Find surface index of 2nd pore aperture and set that surface element
	KSUR=PC(10)+6
	CALL SRT_SETF(KSUR,1,13,PP,0,0,0,KSUR+1,ISTAT)
C Get surface quality of pore walls
	IQ=INT(PC(11))
C Set sides of pore - local x is along pore axis
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
	CALL SRT_SETF(KSUR+2,5,13,PP,0,IQ,KSUR+1,-1,ISTAT)
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
