*+XX_HREF	Calculate X-ray reflectivity from optical constants (Henke)
	SUBROUTINE XX_HREF(GRAZ,ALPHA,GAMMA,RSIG,RPI)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*GRAZ	input	grazing angle radians
*ALPHA	input	real part of dielectric constant
*GAMMA	input	imaginary part of dielectric constant
*RSIG	output	reflectivity sigma polarization
*RPI	output	reflectivity pi polarization
*-Author Dick Willingale 1986-Jun-4
*Uses exact form of Fresnel equations given by Henke et al.
	DOUBLE PRECISION PI
	PARAMETER(PI=3.1415926535898)
C
	THETA=GRAZ
	IF(THETA.LE.0.0D0) THEN
		RSIG=1.0D0
		RPI=1.0D0
		RETURN
	ENDIF
	ST=SIN(THETA)
	CT=COS(THETA)
	TT=TAN(THETA)
C
	SALP=ST**2-ALPHA
	G2=GAMMA**2
	ROW=SQRT((SALP+SQRT(SALP**2+G2))*0.5D0)
	SPROW=ST+ROW
	SMROW=ST-ROW
	ROW4=4.D0*ROW**2
	RSIG=(ROW4*SMROW**2+G2)/(ROW4*SPROW**2+G2)
	CC=(ROW-CT/TT)**2
	CPT=(ROW+CT/TT)**2
C
	RPI=RSIG*(ROW4*CC+G2)/(ROW4*CPT+G2)
	RETURN
	END
