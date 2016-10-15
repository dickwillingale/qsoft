*+SPF_BPLAW	Broken power law continuum spectrum
	FUNCTION SPF_BPLAW(E,ALPHA,BETA,EB,A)
	REAL SPF_BAND,E,ALPHA,BETA,EB,A
*E		Energy keV
*ALPHA		first photon index
*BETA		second photon index
*EB		break energy keV
*A		Normalisation
*SPF_BPLAW	continuum at E (photons/keV)
*-Author Dick Willingale 2007-Mar-15
	REAL DIFF,R
	DIFF=ALPHA-BETA
	IF(E.LT.EB) THEN
		R=E**(-ALPHA)
	ELSE
		R=(E**(-BETA))*(EB**(-DIFF))
	ENDIF
	SPF_BPLAW=A*R
	END
