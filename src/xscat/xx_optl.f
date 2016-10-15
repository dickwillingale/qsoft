*+XX_OPTL	 Calculate X-ray optical properties of a medium
	SUBROUTINE XX_OPTL(ITYPE,IN,IO,ICOMP,RHO,NAT,XRAY,ALPHA,GAMMA,ABSL,
     +	F1,F2)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	INTEGER ITYPE,IN,IO,NAT
	INTEGER ICOMP(NAT)
*ITYPE		input	source of atomic data (0 Cromer, 1 Henke, 2 Hagemann)
*IN		input	atomic data channel number (see XX_GET_ routines)
*IO		input	channel for messages (0 for no messages)
*ICOMP(NAT)	input	composition array
*RHO		input	density (gms/cm**3)
*NAT		input	number of atomic types on IN
*XRAY		input	wavelength Angstroms
*ALPHA		output	real part of dielectric constant
*GAMMA		output	imaginary part of dielectric constant
*ABSL		output	absorption length (cm-1)
*F1		output	real part of scattering factor
*F2		output	imaginary part of scattering factor
*-Author Dick Willingale 1986-Mar-17
*Modified to accept Hagemann data 1986-Jun-2
C Common block to pick up atomic weight from Cromer
	COMMON/ATOMIC/ATWT,ATOM,ANU
C Zero values before looping for atomic types
	F1=0.D0
	F2=0.D0
	ABSL=0.D0
	RMOLWT=0.D0
	LEMENT=1
	DO J=1,NAT
		IF(ITYPE.EQ.0) THEN
			CALL XX_CROMER(IN,IO,XRAY,0.0D0,F0,FP,FPP,ABS)
		ELSEIF(ITYPE.EQ.1) THEN
			   CALL XX_HENKE(IN,IO,XRAY,ATWT,LEMENT,FP,FPP,ABS)
			   F0=0.
		ELSEIF(ITYPE.EQ.2) THEN
			CALL XX_HAGEMANN(IN,IO,XRAY,ATWT,FP,FPP,ABS)
			F0=0.
		ENDIF
		COMP=ICOMP(J)
		RMOLWT=RMOLWT+ATWT*COMP
		F1=F1+(F0+FP)*COMP
		F2=F2+FPP*COMP
		ABSL=ABSL+ABS*COMP*ATWT
		LEMENT=LEMENT+1
	ENDDO
	A=(5.4018D-6)*RHO*XRAY**2
	ABSL=ABSL*RHO/RMOLWT
	ALPHA=F1*A/RMOLWT
	GAMMA=F2*A/RMOLWT
	IF(IO.GT.0) THEN
		WRITE(IO,*) 'Alpha ',ALPHA,'  Gamma ',GAMMA
	ENDIF
	RETURN
	END
