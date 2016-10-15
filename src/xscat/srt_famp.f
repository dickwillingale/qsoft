*+SRT_FAMP	Fresnel reflection and transmission amplitudes
	SUBROUTINE SRT_FAMP(S0,N1,N2,RS,TS,RP,TP,C2)
	IMPLICIT NONE
	DOUBLE COMPLEX S0,N1,N2,RS,TS,RP,TP,C2
*S0	input	sine of angle of incidence at free space boundary
*N1	input	refractive index of medium 1 
*N2	input	refractive index of medium 2
*RS	output	sigma reflection amplitude
*TS	output	sigma transmission amplitude
*RP	output	pi reflection amplitude
*TP	output	pi transmission amplitude
*C2	output	cosine of refraction angle
* Use general form of the Fresnel Equations with complex refractive
* indices and complex incidence and refraction angles
*-Author Dick Willingale 1997-Feb-24
	DOUBLE COMPLEX C1,N1N2,S1,S2,AA,BB,CC
C Use Snell's law to calculate angle of incidence in medium 1
	S1=S0/N1
	C1=SQRT(CMPLX(1.0,0.0)-S1**2)
C Find ratio of refractive indices n1/n2
	N1N2=N1/N2
C Use Snell's law to get sin and cos angle in medium 2
	S2=S1*N1N2
	C2=SQRT(CMPLX(1.0,0.0)-S2**2)
C Calculate sigma amplitude coefficients
	CC=N1N2*C1
	BB=CC+C2
	IF(ABS(BB).NE.0.0) THEN
		RS=(CC-C2)/BB
		TS=CMPLX(2.0,0.0)*CC/BB
	ELSE
		RS=CMPLX(1.0,0.0)
		TS=CMPLX(0.0,0.0)
	ENDIF
C Calculate pi amplitude coefficients
	AA=N1N2*C2
	BB=AA+C1
	IF(ABS(BB).NE.0.0) THEN
		RP=(AA-C1)/BB
		TP=CMPLX(2.0,0.0)*CC/BB
	ELSE
		RP=CMPLX(1.0,0.0)
		TP=CMPLX(0.0,0.0)
	ENDIF
	END
