*+XX_REF	Calculate X-ray reflectivity from optical constants
	SUBROUTINE XX_REF(GRAZ,DELTA,BETA,RSIG,RPI)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*GRAZ	input	grazing angle
*DELTA	input	real part of refractive index
*BETA	input	imaginary part of refractive index
*RSIG	output	reflectivity sigma polarization
*RPI	output	reflectivity pi polarization
*-Author Dick Willingale 1986-Jan-10
C Modified to use DP 1986-Feb-7
	REAL PI/3.1415927/
C Uses exact form of Fresnel equations
	ALF=2.*DELTA-DELTA**2+BETA**2
	GAM=2.*(1.-DELTA)*BETA
	GAM2=GAM**2
	SING=SIN(GRAZ)
	SING2=SING**2
	COSG=COS(GRAZ)
	COSG2=COSG**2
	A2=.5*(SING2-ALF+SQRT((SING2-ALF)**2+GAM2))
	A=SQRT(A2)
	RSIG=(4.*A2*(SING-A)**2+GAM2)/(4.*A2*(SING+A)**2+GAM2)
	RPI=RSIG*(4.*A2*(A-COSG2/SING)**2+GAM2)
	RPI=RPI/(4.*A2*(A+COSG2/SING)**2+GAM2)
	RETURN
	END
