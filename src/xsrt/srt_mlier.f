*+SRT_MLIER	Calculate integrated energy response of a multilayer
	SUBROUTINE SRT_MLIER(DMIN,DMAX,HRAT,NPER,NE,EKEV,NA,ANGS,
     +	HNR,HNI,LNR,LNI,NREF,ARRAY,ISTAT)
	IMPLICIT NONE
	INTEGER NPER,NE,NA,NREF,ISTAT
	DOUBLE PRECISION DMIN,DMAX,HRAT,EKEV(NE),ANGS(NA),HNR(NE),HNI(NE)
	DOUBLE PRECISION LNR(NE),LNI(NE),ARRAY(NE,NA)
*DMIN	input	minimum d-spacing A
*DMAX	input	maximum d-spacing A
*HRAT	input	heavy thickness fraction
*NPER	input	number of periods (0 for a single interface)
*NE	input	number of energies
*EKEV	input	array of energies keV
*NA	input	number of incidence angles
*ANGS	input	incidence angles (degrees)
*HNR	input	real part of refractive index of heavy layers
*HNI	input	imaginary part of refractive index of heavy layers
*LNR	input	real part of refractive index of light layers
*LNI	input	imaginary part of refractive index of light layers
*NREF	input	number of reflections 
*ARRAY	output	array of integrated energy responses
*ISTAT	in/out	returned status
*Note that the light element values are ignored if NPER=0
*-Author Dick Willingale 1997-Apr-29
	INCLUDE 'SRT_COM'
	DOUBLE PRECISION ANG,WL,DWL,D(5),SA,RSIG,RPI,TSIG,TPI
	DOUBLE PRECISION ELO,EHI,ESAM,CONV,EN
	INTEGER I,J,K,NLAY,NR
	DOUBLE COMPLEX FR(5)
C 
	IF(ISTAT.NE.0.0) RETURN
C Set number of layers
	IF(NPER.EQ.0) THEN
		NLAY=2
	ELSE
		NLAY=5
	ENDIF
C Number of samples to be used for integration at each energy
	NR=200
C Conversion from keV to A^(-1)
	CONV=12.39854
C Set complex refractive index of first medium
	FR(1)=CMPLX(1.0,0.0)
C Set thicknesses of first medium and substrate
	D(1)=0.0
	D(5)=0.0
C Loop for incidence angles
	DO J=1,NA
C Convert angle to radians
		ANG=ANGS(J)*PI/180.0
		SA=SIN(PI*0.5-ANG)
C Loop for energies
		DO K=1,NE
C Convert kev to wavelength A
			WL=CONV/EKEV(K)	
C Calculate d spacing for wavelength angle combination
			IF(SA.NE.0.0) THEN
				DWL=WL*0.5/SA
			ELSE
				DWL=0.0
			ENDIF
			DWL=MIN(DWL,DMAX)
C Initialize integrated response
			ARRAY(K,J)=0.0
			IF(NLAY.EQ.5.AND.DWL.GT.DMIN) THEN
C Calculate integrated reflectivity over range of energies
				ELO=EKEV(K)*0.5
				EHI=EKEV(K)*1.5
				ESAM=(EHI-ELO)/(NR-1)
C Set thicknesses of periodic layer
				D(2)=DWL*HRAT
				D(3)=DWL*(1.0-HRAT)
				D(4)=D(2)
				DO I=1,NR
					EN=DBLE(I-1)*ESAM+ELO
					WL=CONV/EN	
C Set refractive indices of layers
					CALL SRT_LIRI(NE,EKEV,HNR,HNI,LNR,LNI,
     +					EN,FR(2),FR(3),ISTAT)
					FR(4)=FR(2)
					FR(5)=FR(3)
					CALL SRT_MLTI(ANG,WL,NLAY,FR,D,
     +					NPER,RSIG,RPI,TSIG,TPI,ISTAT)
					ARRAY(K,J)=ARRAY(K,J)+
     +					((RSIG+RPI)*0.5)**NREF
				ENDDO
C Weight integral by sin of grazing angle
				ARRAY(K,J)=ARRAY(K,J)*SA*ESAM/EKEV(K)
			ELSEIF(NLAY.EQ.2) THEN
C Fresnel reflectivity from  a single interface
				D(2)=0.0
				FR(2)=CMPLX(HNR(K),HNI(K))
				CALL SRT_MLTI(ANG,WL,NLAY,FR,D,
     +				NPER,RSIG,RPI,TSIG,TPI,ISTAT)
				ARRAY(K,J)=SA*((RSIG+RPI)*0.5)**NREF
			ENDIF
		ENDDO
		write(*,*) 'done angle',ang
	ENDDO	
	END
