*+SPF_PLAW	Simple power law
	FUNCTION SPF_PLAW(E,ALPHA,A)
	REAL SPF_PLAW,E,T,A
*E		Energy keV
*ALPHA		Index (photon)
*A		Normalisation
*SPF_PLAW	power law continuum at E (photons/keV)
*-Author Dick Willingale 1986-Sep-4
C Trap very large index
	IF(ALPHA.LT.10.) THEN
		R=1./(E**ALPHA)
	ELSE
		R=0.
	ENDIF
	SPF_PLAW=A*R
	END

