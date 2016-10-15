*+SRT_FLEG	Fourier-Legendre function for cylindirical coordinates
	FUNCTION SRT_FLEG(X,TH,C,NN,MM,IFAIL)
	IMPLICIT NONE
	INTEGER NN,MM,IFAIL
	DOUBLE PRECISION SRT_FLEG,X,TH,C(NN,MM)
*X	axial position (normalised to range -1 to +1)
*TH	azimuth radians
*C	coefficients
*NN	maximum index for Legendre Polynomial
*MM	maximum index for Fourier coefficients (COS + SIN terms)
*IFAIL	status returned, 0 if OK
*-Author Dick Willingale 1996-May-28
	DOUBLE PRECISION SQ2,F1,F2,F3,T1,T2,T3,T4,T5,T6
C
	IF(NN.NE.3.OR.MM.NE.7) THEN
		IFAIL=1
		write(*,*) 'srt_fleg failed - wrong dimensions for coeffs'
		RETURN
	ENDIF
C
	SQ2=SQRT(2.0D0)
	F1=SQ2
	F2=SQRT(6.0D0)*X
	F3=SQRT(10.0D0)*(3*X*X-1.0D0)*0.5D0
C
	T1=COS(TH)
	T2=COS(TH*2.D0)
	T3=COS(TH*3.D0)
	T4=SIN(TH)
	T5=SIN(TH*2.D0)
	T6=SIN(TH*3.D0)
C
	SRT_FLEG=(C(1,1)*F1+C(2,1)*F2+C(3,1)*F3)/SQ2
	SRT_FLEG=SRT_FLEG+(C(1,2)*F1+C(2,2)*F2+C(3,2)*F3)*T1
	SRT_FLEG=SRT_FLEG+(C(1,3)*F1+C(2,3)*F2+C(3,3)*F3)*T2
	SRT_FLEG=SRT_FLEG+(C(1,4)*F1+C(2,4)*F2+C(3,4)*F3)*T3
	SRT_FLEG=SRT_FLEG+(C(1,5)*F1+C(2,5)*F2+C(3,5)*F3)*T4
	SRT_FLEG=SRT_FLEG+(C(1,6)*F1+C(2,6)*F2+C(3,6)*F3)*T5
	SRT_FLEG=SRT_FLEG+(C(1,7)*F1+C(2,7)*F2+C(3,7)*F3)*T6
	END	
