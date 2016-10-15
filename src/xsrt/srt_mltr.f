*+SRT_MLTR	Calculate X-ray reflectivities of a multilayer
	SUBROUTINE SRT_MLTR(NG,ANGS,WL,NLAY,NR,NI,D,NPER,FR,
     +	RSIG,RPI,TSIG,TPI,ISTAT)
	IMPLICIT NONE
	INTEGER NG,NLAY,ISTAT,NPER
	DOUBLE PRECISION ANGS(NG),WL,NR(NLAY),NI(NLAY),D(NLAY)
	DOUBLE PRECISION RSIG(NG),RPI(NG),TSIG(NG),TPI(NG)
	DOUBLE COMPLEX FR(NLAY)
*ANGS	input	incidence angles (degrees)
*WL	input	wavelength (same units as D)
*NLAY	input	number of layers (including substrate if there is one)
*NR	input	real part of refractive index of layers
*NI	input	imaginary part of refractive index of layers
*D	input	thickness of layers (same units as WL)
*NPER	input	number of periods
*FR	output	complex refractive index of layers
*RSIG	output	reflectivity sigma polarization
*RPI	output	reflectivity pi polarization
*TSIG	output	transmission sigma polarization
*TPI	output	transmission pi polarization
*ISTAT	in/out	returned status
*-Author Dick Willingale 1997-Feb-17
	INCLUDE 'SRT_COM'
	DOUBLE PRECISION ANG
	INTEGER J
C 
	IF(ISTAT.NE.0.0) RETURN
C Set complex refractive index
	DO J=1,NLAY
		FR(J)=CMPLX(NR(J),NI(J))
	ENDDO
C Loop for incidence angles
	DO J=1,NG
C Convert angle to radians
		ANG=ANGS(J)*PI/180.0
		CALL SRT_MLTI(ANG,WL,NLAY,FR,D,NPER,RSIG(J),RPI(J),TSIG(J),
     +		TPI(J),ISTAT)
	ENDDO	
	END
