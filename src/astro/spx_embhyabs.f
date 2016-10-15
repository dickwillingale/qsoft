	SUBROUTINE SPX_EMBHYABS(CD,EN,NE,FACT)
	INTEGER NE
	REAL CD,EN(NE),FACT(NE)
*CD		input	hydgrogen column 10**21 cm-2
*EN(NE)		input	energies keV
*FACT(NE)	output	absorption factor
*-Author Dick Willingale 2002-May-09
	REAL TAU
	DO J=1,NE
		TAU=SPF_HYTAU(EN(J),CD)
		IF(TAU.NE.0.0) THEN
			FACT(J)=(1.0-EXP(-TAU))/TAU
		ELSE
			FACT(J)=0.0
		ENDIF
	ENDDO
	END
