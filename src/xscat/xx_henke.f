*+XX_HENKE	Interpolate X-ray scattering parameters from Henke data
	SUBROUTINE XX_HENKE(IN,IO,XRAY,ATWT,LEMENT,FP,FPP,ABS)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	INTEGER IN,IO
*IN		input	channel for atomic data
*IO		input	channel for messages (0 for no messages)
*XRAY		input	wavelength (Angstroms)
*ATWT		output	atomic weight
*FP		output	real part of atomic scattering factor
*FPP		output	imaginary part of atomic scattering factor
*ABS		output	absorption (cm**2/gm)
*LEMENT         input   tells the routine which type to work with
*NAT            input   number of atomic types
*Note FP is not corrected for scattering angle term sin(theta)/lambda
*-Author Dick Willingale 1986-Mar-17
	PARAMETER (NMAX=1000)
	PARAMETER (IMAX=10)
	DOUBLE PRECISION ENG(NMAX),F1(NMAX),F2(NMAX),DATA
	DOUBLE PRECISION C,C1,C2,RX,PI,AO
	CHARACTER*2 AA(100)
	COMMON/HENKE/NATOMA(IMAX),ATWTA(IMAX),NSAMA(IMAX),DATA(IMAX,NMAX,3)
	COMMON/HENKEC/AA
C C is velocity of light in atomic units
C C1 constant to convert wavelength to atomic units
C C2 constant to convert barns to cm**2/gm
	DATA C/137.0367/
	DATA C1/0.02721/
	DATA C2/0.602472/	
	DATA AO/2.80022E+7/
	DATA PI/3.14159265/
	CHARACTER NAME*2,TEMP*40
C10	FORMAT(A2,1X,I4,F8.2,I5)
11	FORMAT(1X,'   HENKE ',A2,
     +  '  Wavelength ',G13.6,/,'  Absorption ',G13.6,' F1 ',G13.6,
     +  ' F2 ',G13.6,/)
C Find energy in eV from wavelength in Angstoms
	XR=12.397639E3/XRAY
	NAME=AA(LEMENT)
	NATOM=NATOMA(LEMENT)
	ATWT=ATWTA(LEMENT)
	NSAM=NSAMA(LEMENT)
	DO J=1,NSAM
		ENG(J)=DATA(LEMENT,J,1)
		F1(J)=DATA(LEMENT,J,2)
		F2(J)=DATA(LEMENT,J,3)
C Use logs for absorption and F2 to get better fit to edges.
		F2(J)=LOG(F2(J))
	ENDDO
C Interpolate to find value of absorption length, F1 and F2 at requested
C wavelength.
	IF(XR.LE.ENG(1)) THEN
		ABS=0.0
		FP=0.0
		FPP=0.0
	ELSEIF(XR.GE.ENG(NSAM)) THEN
		ABS=0.0
		FP=0.0
		FPP=0.0
	ELSE
C Use linear interpolation
		DO J=2,NSAM
			IF(XR.LT.ENG(J)) THEN
				J1=J-1
				J2=J
				GOTO 103
			ENDIF
		ENDDO
C Then FP and FPP
103		FP=F1(J1)+(F1(J2)-F1(J1))*(XR-ENG(J1))/(ENG(J2)-ENG(J1))
		FPP=F2(J1)+(F2(J2)-F2(J1))*(XR-ENG(J1))/(ENG(J2)-ENG(J1))
		FPP=EXP(FPP)
C
		ABS=FPP*AO*C2*4.D0*PI/(C*ATWT*XR*1.E-3/C1)
	ENDIF
	IF(IO.GT.0) THEN
		WRITE(6,11) NAME,XRAY,ABS,FP,FPP
	ENDIF
	RETURN
	END
