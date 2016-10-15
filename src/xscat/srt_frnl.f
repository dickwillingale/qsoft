*+SRT_FRNL	Calculate X-ray reflectivity from optical constants
	SUBROUTINE SRT_FRNL(GRAZ,ALPHA,GAMMA,RSIG,RPI)
	IMPLICIT NONE
	DOUBLE PRECISION GRAZ,ALPHA,GAMMA,RSIG,RPI
*GRAZ	input	grazing angle
*ALPHA	input	real part of dielectric constant
*GAMMA	input	imaginary part of dielectric constant
*RSIG	output	reflectivity sigma polarization
*RPI	output	reflectivity pi polarization
*-Author Dick Willingale 1986-Jun-4
*Uses exact form of Fresnel equations given by Henke et al.
	DOUBLE PRECISION THETA,ST,CT,TT,SALP,G2
	DOUBLE PRECISION ROW,SPROW,SMROW,ROW4,CC,CPT
C
	THETA=GRAZ
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
	RPI=RSIG*(ROW4*CC+G2)/(ROW4*CPT+G2)
	END
