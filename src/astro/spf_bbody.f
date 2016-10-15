*+SPF_BBODY	Black body spectrum
	FUNCTION SPF_BBODY(E,T,A)
	REAL SPF_BBODY,E,T,A
*E		Energy keV
*T		Temperature keV
*A		Normalisation
*SPF_BBODY	Black body continuum at energy E (photons/keV)
*-Author Dick Willingale 1986-Sep-4
CTrap out very low values in tail
	FAC=E/T
	IF(FAC.LT.30.) THEN
		SPF_BBODY=A*E*E/(EXP(FAC)-1.0)
	ELSE
		SPF_BBODY=0.
	ENDIF
	END

