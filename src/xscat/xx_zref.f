
*+XX_ZREF	Calculate X-ray reflectivity from optical constants (Zombeck)
	SUBROUTINE XX_ZREF(GRAZ,DELTA,BETA,RSIG,RPI)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*GRAZ	input	grazing angle
*DELTA	input	real part of refractive index
*BETA	input	imaginary part of refractive index
*RSIG	output	reflectivity sigma polarization
*RPI	output	reflectivity pi polarization
*-Author Dick Willingale 1986-Jan-14
*Modified to run with CROMER etc. 1986-Feb-7 RW
	REAL PI/3.1415927/
C Use Fresnel Equations given by Zombeck 13-36
	THETA=PI*0.5-GRAZ
	ST=SIN(THETA)
	CT=COS(THETA)
	TT=TAN(THETA)
	RN=1.0-DELTA
	RK=BETA
	X=RN**2-RK**2-ST**2
	Y=2.0*RN*RK
	A2=0.5*(SQRT(X**2+Y**2)+X)
	A=SQRT(A2)
	B2=0.5*(SQRT(X**2+Y**2)-X)
	RSIG=(A2+B2-2.0*A*CT+CT**2)/(A2+B2+2.0*A*CT+CT**2)
	RPI=RSIG*(A2+B2-2.0*A*ST*TT+ST**2*TT**2)
	RPI=RPI/(A2+B2+2.0*A*ST*TT+ST**2*TT**2)
	RETURN
	END
