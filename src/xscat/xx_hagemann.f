*+XX_HAGEMANN	Interpolate optical constants from Hagemann et al
	SUBROUTINE XX_HAGEMANN(IN,IO,XRAY,ATWT,FP,FPP,ABS)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	INTEGER IN,IO
*IN		input	channel for atomic data
*IO		input	channel for messages (0 for no messages)
*XRAY		input	wavelength (Angstroms)
*ATWT		output	atomic weight
*FP		output	real part of atomic scatterimg factor
*FPP		output	imaginary part of atomic scattering factor
*ABS		output	absorption 10E5 cm-1
*Note FP  is not corrected for scattering angle term sin(theta)/lambda
*-Author Dick Willingale 1986-Jun-2
	DOUBLE PRECISION KWRK(200),W(2000),CWRK(200)
	DOUBLE PRECISION ENG(200),AB(200),DD(200),BB(200)
	CHARACTER NAME*2
10	FORMAT(A10,F10.5,F10.5,I10)
11	FORMAT(1X,'   Hagemann ',A2,' ATWT ',F10.3,' RHO ',F10.3,
     +  '  Wavelength ',G13.6,/,'  Absorption ',G13.6,' ALPHA ',G13.6,
     +  ' GAMMA ',G13.6,/)
C Find energy in eV from wavelength in Angstoms
	XR=12.397639E3/XRAY
	READ(IN,10) NAME,ATWT,RHO,NSAM
	DO J=1,NSAM
		READ(IN,*) ENG(J),AB(J),DD(J),BB(J)
C Use logs for absorption and imaginary part of dielectric constant to
C get better match to edges.
		AB(J)=LOG(AB(J))
		BB(J)=LOG(BB(J))
	ENDDO
C Interpolate to find value of absorption length, ALPHA and GAMMA at requested
C wavelength.
	IF(XR.LE.ENG(1)) THEN
		ABS=0.0
		ALPHA=0.0
		GAMMA=0.0
	ELSEIF(XR.GE.ENG(NSAM)) THEN
		ABS=0.0
		ALPHA=0.0
		GAMMA=0.0
	ELSE
C Use linear interpolation
C First absorption
		DO J=2,NSAM
			IF(XR.LT.ENG(J)) THEN
				J1=J-1
				J2=J
				GOTO 103
			ENDIF
		ENDDO
103		ABS=AB(J1)+(AB(J2)-AB(J1))*(XR-ENG(J1))/(ENG(J2)-ENG(J1))
		ABS=EXP(ABS)
C Then alpha and gamma
		ALPHA=DD(J1)+(DD(J2)-DD(J1))*(XR-ENG(J1))/(ENG(J2)-ENG(J1))
		GAMMA=BB(J1)+(BB(J2)-BB(J1))*(XR-ENG(J1))/(ENG(J2)-ENG(J1))
		GAMMA=EXP(GAMMA)
C		LWRK=2000
C		LCK=200
C		KI=NSAM
C		KPRIME=KI+4
C 		IFAIL=0
C		CALL E01BAF(KI,ENG,AB,KWRK,CWRK,LCK,W,LWRK,IFAIL)
C		CALL E02BBF(KPRIME,KWRK,CWRK,XR,ABS,IFAIL)
C		ABS=EXP(ABS)
C 		IFAIL=0
C		CALL E01BAF(KI,ENG,DD,KWRK,CWRK,LCK,W,LWRK,IFAIL)
C		CALL E02BBF(KPRIME,KWRK,CWRK,XR,ALPHA,IFAIL)
C 		IFAIL=0
C		CALL E01BAF(KI,ENG,BB,KWRK,CWRK,LCK,W,LWRK,IFAIL)
C		CALL E02BBF(KPRIME,KWRK,CWRK,XR,GAMMA,IFAIL)
	ENDIF
	IF(IO.GT.0) THEN
		WRITE(6,11) NAME,ATWT,RHO,XRAY,ABS,ALPHA,GAMMA
	ENDIF
C
	A=(5.4018D-6)*RHO*XRAY**2
	FP=ALPHA*ATWT/A
	FPP=GAMMA*ATWT/A
	ABS=ABS/RHO
	IF(IO.GT.0) THEN
		WRITE(6,*) 'FP ',FP,'     FPP ',FPP
	ENDIF
	RETURN
	END
