*+SRT_FPOR	Calculate reflectivities and transmissions of interface
	SUBROUTINE SRT_FPOR(NG,ANGS,NR1,NI1,NR2,NI2,RSIG,RPI,
     +	TSIG,TPI,ISTAT)
	IMPLICIT NONE
	INTEGER NG,ISTAT
	DOUBLE PRECISION ANGS(NG),NR1,NI1,NR2,NI2
	DOUBLE PRECISION RSIG(NG),RPI(NG),TSIG(NG),TPI(NG)
*NG	input	number of incidence angles
*ANGS	input	incidence angles (degrees)
*NR1	input	real part of refractive index of medium 1
*NI1	input	imaginary part of refractive index medium 1
*NR2	input	real part of refractive index of medium 2
*NI2	input	imaginary part of refractive index medium 2
*RSIG	output	reflectivity sigma polarization
*RPI	output	reflectivity pi polarization
*TSIG	output	transmission sigma polarization
*TPI	output	transmission pi polarization
*ISTAT	in/out	returned status
*-Author Dick Willingale 1997-Feb-17
	INCLUDE 'SRT_COM'
	REAL ANG
	INTEGER J
	COMPLEX FR1,FR2,S0
C 
	IF(ISTAT.NE.0.0) RETURN
C Set complex refractive index
	FR1=CMPLX(REAL(NR1),REAL(NI1))
	FR2=CMPLX(REAL(NR2),REAL(NI2))
C Loop for incidence angles
	DO J=1,NG
C Convert angle to radians
		ANG=ANGS(J)*PI/180.0
		S0=SIN(ANG)*FR1
		CALL SRT_FPOW(S0,FR1,FR2,RSIG(J),TSIG(J),RPI(J),TPI(J),ISTAT)
	ENDDO	
	END
