*SRT_LIRI	Linear interpolation of refractive index values
	SUBROUTINE SRT_LIRI(NE,EKEV,HNR,HNI,LNR,LNI,EN,HFR,LFR,ISTAT)
	IMPLICIT NONE
	INTEGER NE,ISTAT
	DOUBLE PRECISION EKEV(NE),HNR(NE),HNI(NE),LNR(NE),LNI(NE),EN
	DOUBLE COMPLEX HFR,LFR
*NE	input	number of energies
*EKEV	input	array of energies (strictly ascending or descending)
*HNR	input	real part of heavy element refractive index
*HNI	input	imaginary part of heavy element refractive index
*LNR	input	real part of light element refractive index
*LNI	input	imaginary part of light element refractive index
*EN	input	energy at which refractive index required
*HFR	output	refractive index of heavy element at energy EN
*LFR	output	refractive index of light element at energy EN
*ISTAT	in/out	returned status
*-Author Dick Willingale 1997-Apr-29
	INTEGER J,JJ
	REAL F
	LOGICAL ASCEND
C
	IF(ISTAT.NE.0) RETURN
C	
	J=1
	JJ=0
	ASCEND=(EKEV(2).GT.EKEV(1))
	DO WHILE(J.LT.NE.AND.JJ.EQ.0)	
		J=J+1
		IF(ASCEND) THEN
			IF(EN.LT.EKEV(J).OR.J.EQ.NE) JJ=J
		ELSE
			IF(EN.GE.EKEV(J).OR.J.EQ.NE) JJ=J
		ENDIF
	ENDDO
C Interpolate or extrapolate refractive index values
	F=(EN-EKEV(JJ-1))/(EKEV(JJ)-EKEV(JJ-1))
	HFR=HNR(JJ-1)+F*(HNR(JJ)-HNR(JJ-1))
	LFR=LNR(JJ-1)+F*(LNR(JJ)-LNR(JJ-1))
	END
