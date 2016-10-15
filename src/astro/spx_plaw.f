	SUBROUTINE SPX_PLAW(AN,PI,BS,EN,NE,FLUX)
	INTEGER NE
	REAL AN,PI,BS,EN(NE),FLUX(NE)
*AN		input	normalization
*PI		input	photon index
*BS		input	blueshift
*EN(NE)		input	array of energies keV
*FLUX(NE)	output	spectrum
*-Author Dick Willingale 1991-Apr-22
	DO J=1,NE
		FLUX(J)=SPF_PLAW(EN(J)*BS,PI,AN)
	ENDDO
	END
