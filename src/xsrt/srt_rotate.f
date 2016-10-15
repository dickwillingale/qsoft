*+SRT_ROTATE	Rotate surface element
	SUBROUTINE SRT_ROTATE(IS,PL,AX,ANGLE,ISTAT)
	IMPLICIT NONE
	INTEGER IS,ISTAT
	DOUBLE PRECISION PL(3),AX(3),ANGLE
*IS	input	surface element number
*PL	input	position of rotation centre
*AX	input	axis of rotation
*ANGLE	input	rotation angle about AX (degrees)
*ISTAT	in/out	returned status
*-Author Dick Willingale 2003-Sep-11
	INCLUDE 'SRT_COM'
	INTEGER I
	DOUBLE PRECISION POS(3),ANM(3),AXL(3),PP(3)
	DOUBLE PRECISION X(3),Y(3),RADS,RN,CA,SA
	DOUBLE PRECISION CP,R,PZ
C
	IF(ISTAT.NE.0) RETURN
C
	IF(IS.LT.1.OR.IS.GT.NSUR) THEN
		WRITE(*,*) 'SRT_ROTATE error - element ',IS,' not set'
		ISTAT=1
		RETURN
	ENDIF
C get position, normal and reference axis of element
	I=IPAR(IS)
	ANM(1)=PAR(I)
	ANM(2)=PAR(I+1)
	ANM(3)=PAR(I+2)
	AXL(1)=PAR(I+3)
	AXL(2)=PAR(I+4)
	AXL(3)=PAR(I+5)
	POS(1)=PAR(I+6)
	POS(2)=PAR(I+7)
	POS(3)=PAR(I+8)
C calculate sin and cosine of rotation angle
	RADS=ANGLE*ASIN(1.D0)/90.D0
	CA=COS(RADS)
	SA=SIN(RADS)
C change position
	PP(1)=POS(1)-PL(1)
	PP(2)=POS(2)-PL(2)
	PP(3)=POS(3)-PL(3)
	RN=PP(1)**2+PP(2)**2+PP(3)**2
	CALL SRT_VDOT(AX,PP,PZ)
	R=SQRT(RN-PZ**2)
	IF(R.GT.0.0) THEN
		CALL SRT_VNRM(PP,ISTAT)
		CALL SRT_VCRS(AX,PP,X)
		RN=X(1)**2+X(2)**2+X(3)**2
		IF(RN.GT.0.0) THEN
			CALL SRT_VCRS(AX,X,Y)
			POS(1)=POS(1)+X(1)*R*SA+Y(1)*R*(1.0D0-CA)
			POS(2)=POS(2)+X(2)*R*SA+Y(2)*R*(1.0D0-CA)
			POS(3)=POS(3)+X(3)*R*SA+Y(3)*R*(1.0D0-CA)
		ENDIF
	ENDIF
C change normal
	CALL SRT_VDOT(AX,ANM,CP)
	R=SQRT(1.0D0-CP**2)
	CALL SRT_VCRS(AX,ANM,X)
	RN=X(1)**2+X(2)**2+X(3)**2
	IF(RN.GT.0) THEN
		CALL SRT_VCRS(AX,X,Y)
		ANM(1)=ANM(1)+X(1)*SA*R+Y(1)*R*(1.0D0-CA)
		ANM(2)=ANM(2)+X(2)*SA*R+Y(2)*R*(1.0D0-CA)
		ANM(3)=ANM(3)+X(3)*SA*R+Y(3)*R*(1.0D0-CA)
	ENDIF
C change reference axis
	CALL SRT_VDOT(AX,AXL,CP)
	R=SQRT(1.0D0-CP**2)
	CALL SRT_VCRS(AX,AXL,X)
	RN=X(1)**2+X(2)**2+X(3)**2
	IF(RN.GT.0) THEN
		CALL SRT_VCRS(AX,X,Y)
		AXL(1)=AXL(1)+X(1)*SA*R+Y(1)*R*(1.0D0-CA)
		AXL(2)=AXL(2)+X(2)*SA*R+Y(2)*R*(1.0D0-CA)
		AXL(3)=AXL(3)+X(3)*SA*R+Y(3)*R*(1.0D0-CA)
	ENDIF
C update position, normal and reference axis of element
	I=IPAR(IS)
	PAR(I)=ANM(1)
	PAR(I+1)=ANM(2)
	PAR(I+2)=ANM(3)
	PAR(I+3)=AXL(1)
	PAR(I+4)=AXL(2)
	PAR(I+5)=AXL(3)
	PAR(I+6)=POS(1)
	PAR(I+7)=POS(2)
	PAR(I+8)=POS(3)
	END
