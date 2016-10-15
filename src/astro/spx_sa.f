*+SPX_SA	Stellar atmospheres
	SUBROUTINE SPX_SA(SAT,SV,SR,SD,ST,SG,SP,SH,NE,EN,FLUX,
     +	FT,FG,FP,FH,ISTAT)
	CHARACTER SAT*(*)
	INTEGER NE,ISTAT
	REAL SV,SR,SD,ST,SG,SP,SH,EN(NE),FLUX(NE),FT,FG,FP,FH
*SAT	input	type of stellar atmosphere (file name excluding .DAT and dir)
*SV	input	V magnitude used for normalisation if SR=0.0
*SR	input	stellar radius (Solar radii) used for normalisation if >0.0
*SD	input	stellar distance parsecs used for normalisation if SR>0.0
*ST	input	stellar temperature Kelvin (nearest model chosen)
*SG	input	log stellar surface gravity (nearest model chosen)
*SP	input	log stella pressure
*SH	input	He/H ratio
*NE	input	number of energy samples
*EN	input	energy samples, keV
*FLUX	output	stellar flux at Earth unabsorbed photons cm-2 s-1 keV-1
*FT	output	stellar temperature Kelvin found
*FG	output	log stellar temperature found
*FP	output	log stellar pressure found
*FH	output	log stellar He/H ratio found
*ISTAT	in/out	returned status
*-Author Dick Willingale 1992-Nov-11
	PARAMETER (PCONST=6.6262E-27)	!Plancks constant erg sec
	PARAMETER (SRAD=6.9599E10)	!Solar radius cm
	PARAMETER (PCTOCM=3.0856E18)	!parsec to cm
	PARAMETER (AKVERG=1.602E-9)	!keV to erg
	PARAMETER (VKEV=2.2563E-3)	!V band keV
	CHARACTER SPECFILE*80
	PARAMETER (MAXSPEC=4000)	!Maximum number of energy samples
	REAL SPEC(MAXSPEC),ECENTRE(MAXSPEC)
	LOGICAL SEARCH
	CHARACTER SATU*40
C
	IF(ISTAT.NE.0) RETURN
C
	SATU=SAT
	CALL SYS_UPCASE(SATU)
	SPECFILE =SAT
	CALL SYS_GETLUN(INS,ISTAT)
	OPEN (UNIT=INS,FILE=SPECFILE,STATUS='OLD',IOSTAT=ISTAT)
	IF(ISTAT.NE.0) THEN
		WRITE(*,*) 'SPX_SA failed to open '//SPECFILE
	ENDIF
	IF(ISTAT.NE.0) RETURN
C Read in emergent flux and frequency from file
C T=effective temperature of model
C GLOG=log surface gravity of model
C PLOG=log pressure model
C HEH=He/H ratio model
C NPSEC=number of spectral points
	HEH = 0
	PLOG = 0
	SEARCH=.TRUE.
	DO WHILE(SEARCH)
		IF (SATU.EQ.'HEH') THEN
			READ (INS,99003,IOSTAT=ISTAT) T,GLOG,NSPEC,HEH
		ELSE IF (SATU.EQ.'SHEH') THEN
			READ (INS,99003,IOSTAT=ISTAT) T,GLOG,NSPEC,PLOG
		ELSE
			READ (INS,99001,IOSTAT=ISTAT) T,GLOG,NSPEC
		END IF
		READ (INS,99002,IOSTAT=ISTAT) (ECENTRE(I),SPEC(I),I=1,NSPEC)
99001    	FORMAT (F11.2,F5.2,I5)
99002    	FORMAT (1X,E11.4,E11.4)
99003    	FORMAT (F11.2,F5.2,I5,E11.4)
		IF(ISTAT.NE.0) THEN
			WRITE(*,*) 'SPX_SA failed to find',SAT,ST,SG
		ELSE
			T_DIF = ABS((T-ST)/ST)
			G_DIF = ABS((GLOG-SG)/SG)
		ENDIF
		SEARCH=(ISTAT.EQ.0) .AND.
     +		((T.LT.ST .AND. T_DIF.GT.0.0001) .OR.
     +        	(GLOG.LT.SG .AND. G_DIF.GT.0.05) .OR.
     +		(HEH.LT.SH) .OR.
     +		(PLOG.LT.SP))
		IF(.NOT.SEARCH.AND.ISTAT.EQ.0) THEN
			FT=T
			FG=GLOG
			FP=PLOG
			FH=HEH
C Convert erg cm-2 s-1 Hz-1 to 1E+30 keV cm-2 s-1 keV-1 and frequency
C to keV. Also derive 5500A (2.2563E-3keV) emergent flux for normalisation
C to V magnitude.
			IVFLAG = 0
			DO I = 1,NSPEC
				SPEC(I) = SPEC(I)/(PCONST*1E30)
				ECENTRE(I) = ECENTRE(I)*PCONST/AKVERG
				IF(VKEV.LT.ECENTRE(I).AND.IVFLAG.EQ.0) IVFLAG=I
			END DO
 			IF (IVFLAG.NE.0) THEN
				DL = (ECENTRE(IVFLAG)-VKEV)
     &				/(ECENTRE(IVFLAG)-ECENTRE(IVFLAG-1))
				VEM=SPEC(IVFLAG)*(SPEC(IVFLAG-1)
     +				/SPEC(IVFLAG))**DL
			ELSE
				VEM=0.0
			END IF
C Normalise to flux at earth (keV cm-2 s-1 keV-1)
			IF (SR.EQ.0) THEN
C Normalise using Vmag
C SV is the V magnitude of the star
C PCONST is Plancks constant in ergs sec
C VOBS is the observed flux at 5550 A (V band)
C VEM is the model flux at 5550 A (V band)
        			VOBS = (10.0**(-22.419-0.4*SV))*1.0E3/PCONST
				IF(VEM.NE.0.0) THEN
					CONST = VOBS/VEM
				ELSE
					WRITE(*,*) 'SPX_SA V band not in range'
					ISTAT=1
				ENDIF
			ELSE
C Normalise using stellar radius and distance
C SRAD is Solar radius in cm
C SR is the radius of star in Solar radii
C SD is distance to star in parsec
C DIST is distance to star in cm
				CONST = 1.0E30*((SR*SRAD)/(SD*PCTOCM))**2
			ENDIF
			DO I = 1,NSPEC
				SPEC(I) = SPEC(I)*CONST
			ENDDO
		ENDIF
	ENDDO
	CLOSE(UNIT=INS)
	IF(ISTAT.EQ.0) THEN
C Interpolate and convert to photons cm-2 s-1 keV-1
C Samples outside energy range of file =0.0
		J=1
		I=2
		DO WHILE(J.LE.NE.AND.EN(J).LT.ECENTRE(1))
			FLUX(J)=0.0
			J=J+1
		ENDDO
		DO WHILE(J.LE.NE.AND.EN(J).LE.ECENTRE(NSPEC))
			DO WHILE(I.LE.NSPEC.AND.EN(J).GT.ECENTRE(I))
				I=I+1
			ENDDO
			R=((ECENTRE(I)-EN(J)))/(ECENTRE(I)-ECENTRE(I-1))
			FLUX(J)=SPEC(I)-(SPEC(I)-SPEC(I-1))*R
			FLUX(J)=FLUX(J)/EN(J)
			J=J+1
		ENDDO
		DO WHILE(J.LE.NE)
			FLUX(J)=0.0
			J=J+1
		ENDDO
	ENDIF
	END
