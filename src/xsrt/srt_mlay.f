*+SRT_MLAY	Generate characteristic matrix for an interface	
	SUBROUTINE SRT_MLAY(S0,WV,D,F1,F2,CMS,CMP,TSF,TPF)
	IMPLICIT NONE
	DOUBLE PRECISION WV,D
	DOUBLE COMPLEX S0,F1,F2,CMS(2,2),CMP(2,2)
	LOGICAL TSF,TPF
*S0	input	sine of the incidence angle in free space
*WV	input	free space wavevector
*D	input	thickness of layer 2 (0.0 for infinity, no reflection)
*F1	input	refractive index of layer 1 
*F2	input	refractive index of layer 2
*CMS	output 	characteristic matrix sigma polarization
*CMP	output 	characteristic matrix pi polarization
*TSF	in/out	sigma polarization flag (true for transmission)
*TPF	in/out	pi polarization flag (true for transmission)
*-Author Dick Willingale 1997-Feb-17
	DOUBLE COMPLEX C2,DR,PP,PM,RS,RP,PH,TS,TP
C Calculate Fresnel amplitudes etc. at interface
	CALL SRT_FAMP(S0,F1,F2,RS,TS,RP,TP,C2)
C Calculate phase factor across layer 2
	IF(D.NE.0.0) THEN
		DR=CMPLX(0.0,WV*D)*F2*C2
		PP=EXP(DR)
		PM=EXP(-DR)
	ELSE
		PP=1.0
		PM=1.0
	ENDIF
C Set matrix coefficients for sigma polarization at interface
	PH=CMPLX(1.0,0.0)
	IF(ABS(TS).GT.0.0) THEN
		PH=PH/TS
		CMS(1,1)=PP*PH
		CMS(1,2)=RS*PM*PH
		CMS(2,1)=RS*PP*PH
		CMS(2,2)=PM*PH
	ELSE
		TSF=.FALSE.
	ENDIF
C Set matrix coefficients for pi polarization at interface
	PH=CMPLX(1.0,0.0)
	IF(ABS(TP).GT.0.0) THEN
		PH=PH/TP
		CMP(1,1)=PP*PH
		CMP(1,2)=RP*PM*PH
		CMP(2,1)=RP*PP*PH
		CMP(2,2)=PM*PH
	ELSE
		TPF=.FALSE.
	ENDIF
	END
