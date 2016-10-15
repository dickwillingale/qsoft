	SUBROUTINE SPX_BAND(AN,ALPHA,BETA,EC,BS,EN,NE,FLUX)
	INTEGER NE
	REAL AN,ALPHA,BETA,EC,BS,EN(NE),FLUX(NE)
*AN		input	normalization
*ALPHA		input	1st photon index
*BETA		input	2nd photon index
*EC		input	cut-off energy keV
*BS		input	shift
*EN(NE)		input	array of energies keV
*FLUX(NE)	output	spectrum
*-Author Dick Willingale 2007-Mar-15
	DO J=1,NE
		FLUX(J)=SPF_BAND(EN(J)*BS,ALPHA,BETA,EC,AN)
	ENDDO
	END
