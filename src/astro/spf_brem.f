*+SPF_BREM	Bremsstrahlung with Gaunt factor
	FUNCTION SPF_BREM(E,T,A)
	REAL SPF_BREM,E,T,A
*E		Energy keV
*T		Temperature keV
*A		Normalisation
*SPF_BREM	Bremsstrahlung continuum at energy E (photons/keV)
*-Author Dick Willingale 1986-Sep-4
	IF(T.GT.1.E-15) THEN
		SPF_BREM=A*EXP(-E/T)*SPF_GAUNT(E,T)/E
	ELSE
		SPF_BREM=0.
	ENDIF
	END
