*+SRT_FPOW	Fresnel reflection and transmission powers
	SUBROUTINE SRT_FPOW(S0,N1,N2,RSP,TSP,RPP,TPP)
	IMPLICIT NONE
	COMPLEX S0,N1,N2
	DOUBLE PRECISION RSP,TSP,RPP,TPP
*S0	input	sine of angle of incidence at free space boundary
*N1	input	refractive index of medium 1 
*N2	input	refractive index of medium 2
*RSP	output	sigma reflection 
*TSP	output	sigma transmission
*RPP	output	pi reflection 
*TPP	output	pi transmission
*-Author Dick Willingale 1997-Feb-24
	COMPLEX RS,TS,RP,TP,C2
	REAL CC,DD,PR
C
	CALL SRT_FAMP(S0,N1,N2,RS,TS,RP,TP,C2)
	RSP=ABS(RS)**2
	RPP=ABS(RP)**2
C Calculate ratio of Poynting vectors required for transmission
	CC=REAL(SQRT(N1**2-S0**2))
	DD=REAL(SQRT(N2**2-S0**2))
	IF(CC.NE.0.0) THEN
		PR=DD/CC
	ELSE
		PR=0.0
	ENDIF
	TSP=PR*ABS(TS)**2
	TPP=PR*ABS(TP)**2
	END
