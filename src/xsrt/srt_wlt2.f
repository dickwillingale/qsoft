*+SRT_WLT2      Calculate parameters for Wolter Type II
	SUBROUTINE SRT_WLT2(RP,GP,RH,GH,RM,FOVR,AX,AR,FF,PP,PY)
	IMPLICIT NONE
	DOUBLE PRECISION RP,GP,RH,GH,RM,FOVR,AX(3),AR(3),FF(3)
	DOUBLE PRECISION PP(16),PY(16)
*RP	input	maximum radius of paraboloid
*GP	input	grazing angle at RP on paraboloid (radians)
*RH	input	maximum radius of hyperboloid
*GH	input	grazing angle at RH on hyperboloid (radians)
*RM	input	minimum radius on paraboloid
*FOVR	input	radius of FOV radians
*AX	input	axis direction
*AR	input	reference axis (perpendicular to AX)
*FF	input	focus position
*PP     output  parameters for the paraboloid
*PY     output  parameters for the hyperboloid
* In both cases (PP(16) and PY(16) the parameters are:
* r**2=P(10).x**2+P(11).x+P(12) where x is axial distance from focus
* xmin=P(13), rmin=P(14), xmax=P(15), rmax=P(16)
*-Author Dick Willingale 1997-Jul-22
	INTEGER J,IROOT
	DOUBLE PRECISION TANP,TAN2P,TAN2PH,TAN2P2H
	DOUBLE PRECISION D,A2,XH,XP,XM,DX,A,B,C,G,XX
	DOUBLE PRECISION R1,R2,SRT_QUFU
	EXTERNAL SRT_QUFU
C Transfer vectors
	DO J=1,3
		PP(J)=AX(J)
		PY(J)=AX(J)
	ENDDO
	DO J=1,3
		PP(J+3)=AR(J)
		PY(J+3)=AR(J)
	ENDDO
	DO J=1,3
		PP(J+6)=FF(J)
		PY(J+6)=FF(J)
	ENDDO
C 
	TANP=TAN(GP)
	TAN2P=TAN(GP*2.D0)
	TAN2PH=TAN(2.D0*GP-GH)
	TAN2P2H=TAN(2.D0*(GP-GH))
C Calculate quadratic coefficients r**2=A.x**2+B.x+C for paraboloid
	PP(10)=0.0
	PP(11)=RP*TANP*2.D0
	PP(12)=RP**2-PP(11)*(RH/TAN2P2H+(RP-RH)/TAN2P)
C Set limits of parabola
	XM=(RM**2-PP(12))/PP(11)
	PP(13)=XM
	PP(14)=RM
	XP=(RP**2-PP(12))/PP(11)
	PP(15)=XP
	PP(16)=RP
C Calculate quadratic coefficients r**2=A.x**2+B.x+C for hyperboloid
	D=XP-RP/TAN2P
	XH=D+RH/TAN2P
	A2=0.25D0*D**2/(2.D0*RH*TAN2PH/(2.D0*XH-D)+1.D0)
	PY(10)=(D**2-4.D0*A2)/(4.D0*A2)
	PY(11)=-D*PY(10)
	PY(12)=RH**2-PY(10)*XH**2-PY(11)*XH
C Set limits of hyperbola allowing for FOV
	DX=(XP-XH)*FOVR/(2.D0*GP-GH)
	G=TAN(2.0D0*ATAN(PP(11)/(2.D0*RM)))
	XX=RM-G*XM
	A=PY(10)-G**2
	B=PY(11)-2.D0*G*XX
	C=PY(12)-XX**2
	CALL SRT_QURT(A,B,C,IROOT,R1,R2)
	IF(IROOT.EQ.1) THEN
		PY(13)=R1
	ELSEIF(IROOT.EQ.2) THEN
		PY(13)=MAX(R1,R2)
	ENDIF
	PY(13)=PY(13)-DX
C Check that lower limit in x is not beyond the intersection of the
C hyperbola with the axis
	A=PY(10)
	B=PY(11)
	C=PY(12)
	CALL SRT_QURT(A,B,C,IROOT,R1,R2)
	R1=MAX(R1,R2)
	PY(13)=MAX(R1,PY(13))
	PY(14)=SQRT(SRT_QUFU(PY(10),PY(13)))
C Now upper limit
	G=TAN(2.0D0*ATAN(PP(11)/(2.D0*RP)))
	XX=RP-G*XP
	A=PY(10)-G**2
	B=PY(11)-2.D0*G*XX
	C=PY(12)-XX**2
	CALL SRT_QURT(A,B,C,IROOT,R1,R2)
	IF(IROOT.EQ.1) THEN
		PY(15)=R1
	ELSEIF(IROOT.EQ.2) THEN
		PY(15)=MAX(R1,R2)
	ENDIF
	PY(15)=PY(15)+DX
	PY(16)=SQRT(SRT_QUFU(PY(10),PY(15)))
	END
