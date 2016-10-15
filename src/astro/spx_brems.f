	SUBROUTINE SPX_BREMS(AN,TK,BS,EN,NE,FLUX)
	INTEGER NE
	REAL AN,PI,BS,EN(NE),FLUX(NE)
*AN		input	normalisation
*TK		input	temperature keV
*BS		input	blueshift
*EN(NE)		input	energy keV
*FLUX(NE)	output	spectrum
*-Author Dick Willingale 1991-Apr-26
	DO J=1,NE
		FLUX(J)=SPF_BREM(EN(J)*BS,TK,AN)
	ENDDO
	END
