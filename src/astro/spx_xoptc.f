*+SPX_XOPTC	Calculate X-ray properties of a material
	SUBROUTINE SPX_XOPTC(ITYPE,IN,NAT,ICOMP,RHO,GRAZ,NE,EKEV,WL,ALPHA,
     +	GAMMA,ABSL,F1,F2,RSIG,RPI,ISTAT)
	INTEGER ITYPE,IN,NAT,ICOMP(NAT),NE,ISTAT
	DOUBLE PRECISION RHO,GRAZ,EKEV(NE),WL(NE),ALPHA(NE)
	DOUBLE PRECISION GAMMA(NE),ABSL(NE),RSIG(NE),RPI(NE),F1(NE),F2(NE)
*ITYPE		input	source of data (0=Cromer and Liberman)
*IN		input	data channel (from xx_get_cromer)
*NAT		input	number of atomic types on IO
*ICOMP		input	composition
*RHO		input	density (gm/cm**3)
*GRAZ		input	grazing angle (degrees)
*NE		input	number of energies
*EKEV		input	array of energies keV
*WL		output	array wavelengths Angstoms
*ALPHA		output	real part of dielectric constant
*GAMMA		output	imaginary part of dielectric constant
*ABSL		output	absorption length cm-1
*F1		output	real part of scattering factor
*F2		output	imaginary part of scattering factor
*RSIG		output	reflectivity sigma polarization
*RPI		output	reflectivity pi polarization
*ISTAT		in/out	returned status, 0 if OK
*-Author Dick Willingale 1992-Dec-18
	DOUBLE PRECISION AKEV,DTOR,GRAN
	PARAMETER (AKEV=12.397639)
	PARAMETER (DTOR=1.745329252E-2)
C
	IF(ISTAT.NE.0) RETURN
C
	GRAN=GRAZ*DTOR
C Loop for all energies
	DO J=1,NE
		IF(EKEV(J).GT.0.0) THEN
C Convert keV to Angstoms
			WL(J)=AKEV/EKEV(J)
C Calculate optical constants using atomic data
			IO=0
			REWIND IN
			CALL XX_OPTL(ITYPE,IN,IO,ICOMP,RHO,NAT,WL(J),
     +			ALPHA(J),GAMMA(J),ABSL(J),F1(J),F2(J))
C Calculate reflectivites using optical constants
			CALL XX_HREF(GRAN,ALPHA(J),GAMMA(J),RSIG(J),RPI(J))
		ELSE
			WL(J)=0.0
			ALPHA(J)=0.0
			GAMMA(J)=0.0
			ABSL(J)=0.0
			RSIG(J)=0.0
			RPI(J)=0.0
			F1(J)=0.0
			F2(J)=0.0
		ENDIF
	ENDDO
	END
