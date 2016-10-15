*+SPF_BAND	D.Band GRB continuum spectrum
	FUNCTION SPF_BAND(E,ALPHA,BETA,EC,A)
	REAL SPF_BAND,E,ALPHA,BETA,EC,A
*E		Energy keV
*ALPHA		first photon index
*BETA		second photon index
*EC		cut energy keV
*A		Normalisation
*SPF_BAND	continuum at E (photons/keV)
*-Author Dick Willingale 2007-Mar-15
	REAL DIFF,EP,R
	DIFF=ALPHA-BETA
	EP=-DIFF*EC
	IF(E.LT.EP) THEN
		R=(E**(-ALPHA))*EXP(-E/EC)
	ELSE
		R=(E**(-BETA))*EXP(DIFF)*(EP**(-DIFF))
	ENDIF
	SPF_BAND=A*R
	END
