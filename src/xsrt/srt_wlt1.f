*+SRT_WLT1      Calculate parameters for Wolter Type I
	SUBROUTINE SRT_WLT1(XJ,RJ,RA,PL,PH,HL,HH,AX,AR,FF,PP,PY)
	IMPLICIT NONE
	DOUBLE PRECISION XJ,RJ,RA,PL,PH,HL,HH,AX(3),AR(3),FF(3)
	DOUBLE PRECISION PP(16),PY(16)
*XJ     input   axial position of join plane (wrt focal plane)
*RJ     input   radius at join plane
*RA     input   ratio of grazing angles (usually 1.0)
*PL     input   axial position of start of paraboloid (near join)
*PH     input   axial position of end of paraboloid (entrance)
*HL     input   axial position of start of hyperboloid (exit)
*HH     input   axial position of end of hyperboloid (near join)
*AX	input	axis direction
*AR	input	reference axis (perpendicular to AX)
*FF	input	focus position
*PP     output  parameters for the paraboloid
*PY     output  parameters for the hyperboloid
* In both cases (PP(16) and PY(16) the parameters are:
* r**2=P(10).x**2+P(11).x+P(12) where x is axial distance from focus
* xmin=P(13), rmin=P(14), xmax=P(15), rmax=P(16)
*-Author Dick Willingale 1996-Nov-20
	DOUBLE PRECISION SRT_QUFU,TANALP4,COSALP4,ALP4,TANGRP,TANGRH
	DOUBLE PRECISION PWI,EWI2,DWI
	INTEGER J
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
C Calculate nominal grazing angle at join plane
	TANALP4=RJ/XJ
	ALP4=ATAN(TANALP4)
	COSALP4=COS(ALP4)
C Calculate tangents of grazing angles using ratio
	TANGRP=TAN(ALP4*0.5*RA/(1.0+RA))
	TANGRH=TAN(ALP4*0.5*(1.0+2.0*RA)/(1.0+RA))
C Calculate standard parameters
	PWI=RJ*TANGRP
	EWI2=(COSALP4*(1.0+TANALP4*TANGRH))**2
	DWI=RJ*(TANALP4-TANGRH)/(1.0+TANALP4*TANGRH)
C Calculate quadratic coefficients for r**2=A.x**2+B.x+C
	PP(10)=0.0
	PP(11)=2.0*PWI
	PP(12)=PWI**2+4.0*PWI*DWI*EWI2/(EWI2-1.0)
	PP(13)=PL
	PP(14)=SQRT(SRT_QUFU(PP(10),PL))
	PP(15)=PH
	PP(16)=SQRT(SRT_QUFU(PP(10),PH))
	PY(10)=EWI2-1.0
	PY(11)=2.0*DWI*EWI2
	PY(12)=EWI2*DWI**2
	PY(13)=HL
	PY(14)=SQRT(SRT_QUFU(PY(10),HL))
	PY(15)=HH
	PY(16)=SQRT(SRT_QUFU(PY(10),HH))
	END
