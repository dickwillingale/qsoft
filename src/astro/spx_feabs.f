	SUBROUTINE SPX_FEABS(ITYPE,CD,EN,NE,FACT)
	INTEGER ITYPE,NE
	REAL CD,EN(NE),FACT(NE)
*ITYPE		input	cold 1, He like 2, H like 3
*CD		input	hydgrogen column 10**16.52 cm-2
*EN(NE)		input	energies keV
*FACT(NE)	output	absorption factor
*-Author Dick Willingale 1991-Apr-22
	DO J=1,NE
		FACT(J)=SPF_FEABS(EN(J),ITYPE,CD)
	ENDDO
	END
