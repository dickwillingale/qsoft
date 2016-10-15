	SUBROUTINE SPX_HYABS(CD,EN,NE,FACT)
	INTEGER NE
	REAL CD,EN(NE),FACT(NE)
*CD		input	hydgrogen column 10**21 cm-2
*EN(NE)		input	energies keV
*FACT(NE)	output	absorption factor
*-Author Dick Willingale 1991-Apr-22
	DO J=1,NE
		FACT(J)=SPF_HYABS(EN(J),CD)
	ENDDO
	END
