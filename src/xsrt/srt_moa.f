*+SRT_MOA     Dynamically set up slot of a cylindrical MOA
	SUBROUTINE SRT_MOA(DIR,HPOS,CAX,CAR,CVR,PC,IDEF,ISTAT)
	IMPLICIT NONE
	DOUBLE PRECISION DIR(3),HPOS(3),CAX(3),CAR(3),CVR(3),PC(50)
	INTEGER IDEF(2),ISTAT
*DIR	input	direction of ray
*HPOS	input	position of hit on cylindrical reference surface
*CAX    input   axis of cylinder (x axis on MOA surface)
*CAR    input   reference direction (y axis on MOA surface)
*CVR    input   centre of cylinder below aperture
*PC     input   coefficient and limits etc.
*IDEF   input   deformation indices
*ISTAT  in/out  returned status
* Parameters for cylindrical MOA are:
* PC(1) t**2 coeff=0.0
* PC(2) t coeff=0.0
* PC(3) constant coeff = radius**2
* PC(4) xmin
* PC(5) ymin
* PC(6) xmax
* PC(7) ymax
* PC(8) slot pitch in y
* PC(9) rib width between slots in y
* PC(10) depth of slot (thickness of MOA)
* PC(11) this surface index (i.e. index of cylindrical surface)
* PC(12) surface quality for reflecting sides of slots
*-Author Dick Willingale 2007-Feb-01
	DOUBLE PRECISION RNM(3),CEN(3),RAD,XH,YH,POS(3),X,Y,H
	DOUBLE PRECISION Z,XP(3),YP(3),ZP(3),COST
	DOUBLE PRECISION PP(13),HX,HY,HL,DX,DY,DT,DDX,DDY,THETA
	DOUBLE PRECISION WALL, PITCH
	INTEGER I,KSUR,IQ
C Find radius of cylinder
	RAD=SQRT(PC(3))
C Find surface normal at centre of MOA
	CALL SRT_VCRS(CAX,CAR,RNM)
	CALL SRT_VNRM(RNM,ISTAT)
C Find position of centre of MOA on cylindrical surface
	DO I=1,3
		CEN(I)=CVR(I)+RNM(I)*RAD
	ENDDO
C Calculate local coordinates of hit on MOA surface
	DO I=1,3
		POS(I)=HPOS(I)-CEN(I)
	ENDDO
	CALL SRT_VDOT(POS,CAX,XH)
	CALL SRT_VDOT(POS,CAR,YH)
C check ray direction to see if pore inside of outside cylindrical surface
        CALL SRT_VDOT(RNM,DIR,COST)
	THETA=ACOS(ABS(COST))
	IF(COST.LT.0.0) THEN
		THETA=-THETA
	ENDIF
C Get half length (half depth) of slots
	HL=PC(10)*0.5
C Get wall thickness and pitch of slots
	WALL=PC(9)
	IF(COST.LT.0.0) THEN
C Convex surface
		PITCH=PC(8)*(1.0+HL/RAD)
	ELSE
C Concave surface
		PITCH=PC(8)*(1.0-HL/RAD)
	ENDIF
C Use local intersection position to find centre of nearest slot on surface
	X=0.0
	Y=(INT(ABS(YH)/PITCH)+0.5)*PITCH
	IF(YH.LT.0.0) THEN
		Y=-Y
	ENDIF
C X,Y is now the local position of the slot centre on cylindrial surface
C Find surface normal at slot centre
	Z=SQRT(RAD**2-Y**2)
	ZP(1)=Z*RNM(1)+X*CAX(1)+Y*CAR(1)
	ZP(2)=Z*RNM(2)+X*CAX(2)+Y*CAR(2)
	ZP(3)=Z*RNM(3)+X*CAX(3)+Y*CAR(3)
	CALL SRT_VNRM(ZP,ISTAT)
C Find position of slot centre
	POS(1)=CVR(1)+ZP(1)*RAD
	POS(2)=CVR(2)+ZP(2)*RAD
	POS(3)=CVR(3)+ZP(3)*RAD
	IF(IDEF(1).GT.0) THEN
C Get deformations (note deformation gradients are not used)
C DX is effective displacement of centre CVR  in direction CAX
C DY is effective displacement of centre CVR in direction NORM X CAR
		IDEF(2)=1
		CALL SRT_DFRM(IDEF,X,Y,DX,DDX,DDY,ISTAT)
		IDEF(2)=2
		CALL SRT_DFRM(IDEF,X,Y,DY,DDX,DDY,ISTAT)
C Apply deformation to change effective centre of cylinder
		ZP(1)=Z*RNM(1)+(X-DX)*CAX(1)+(Y-DY)*CAR(1)
		ZP(2)=Z*RNM(2)+(X-DX)*CAX(2)+(Y-DY)*CAR(2)
		ZP(3)=Z*RNM(3)+(X-DX)*CAX(3)+(Y-DY)*CAR(3)
		CALL SRT_VNRM(ZP,ISTAT)
	ENDIF
C The slot axis must be flipped if we are at a concave surface
	if(cost.gt.0.0) then
		ZP(1)=-ZP(1)
		ZP(2)=-ZP(2)
		ZP(3)=-ZP(3)
	endif
C Now set up axes tangential to surface at slot centre
	XP(1)=CAX(1)
	XP(2)=CAX(2)
	XP(3)=CAX(3)
	CALL SRT_VCRS(XP,ZP,YP)
        CALL SRT_VNRM(YP,ISTAT)
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
C Set dimensions of slot aperture
	HX=(PC(6)-PC(4))*0.5
	HY=(PITCH-WALL)*0.5
        PP(10)=-HX
        PP(11)=-HY
        PP(12)=HX
        PP(13)=HY
C Find surface index of slot aperture and set that surface element
	KSUR=PC(11)+1
	CALL SRT_SETF(KSUR,1,13,PP,0,0,0,KSUR+1,ISTAT)
C Get surface quality of slot walls
	IQ=INT(PC(12))
C Set sides of slot - local x of sides is along slot axis
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
	CALL SRT_SETF(KSUR+1,5,13,PP,0,IQ,KSUR+1,-1,ISTAT)
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
	CALL SRT_SETF(KSUR+3,5,13,PP,0,IQ,KSUR+1,-1,ISTAT)
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
	CALL SRT_SETF(KSUR+4,5,13,PP,0,IQ,KSUR+1,-1,ISTAT)
	END
