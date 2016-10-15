*+XX_DIELEC	Interpolate optical constants from file XRAY_DATA
	SUBROUTINE XX_DIELEC(XRAY,ALPHA,GAMMA,ABS)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*XRAY		input	wavelength (Angstroms)
*ALPHA		output	real part of dielectric constant
*GAMMA		output	imaginary part of dielectric constant
*ABS		output	absorption 10E5 cm-1
*-Author Dick Willingale 1986-Nov-14
	DOUBLE PRECISION KWRK(200),W(2000),CWRK(200)
	DOUBLE PRECISION ENG(200),AB(200),DD(200),BB(200)
	CHARACTER NAME*2
	DATA IN/0/
C Formats
10	FORMAT(A10,F10.5,F10.5,I10)
11	FORMAT(1X,'   Hagemann ',A2,' ATWT ',F10.3,' RHO ',F10.3,
     +  '  Wavelength ',G13.6,/,'  Absorption ',G13.6,' ALPHA ',G13.6,
     +  ' GAMMA ',G13.6,/)
C Buffer in data if first pass
	IF(IN.EQ.0) THEN
		ISTAT=0
		CALL SYS_GETLUN(IN,ISTAT)
		OPEN(UNIT=IN,FILE='XRAY_DATA',STATUS='OLD',IOSTAT=IERR)
		IF(IERR.NE.0) THEN
			WRITE(*,*) 'XX_DIELEC failed to open XRAY_DATA'
			STOP
		ENDIF
		READ(IN,10) NAME,ATWT,RHO,NSAM
		DO J=1,NSAM
			READ(IN,*) ENG(J),AB(J),DD(J),BB(J)
C Use logs for absorption and imaginary part of dielectric constant to
C get better match to edges.
			AB(J)=LOG(AB(J))
			BB(J)=LOG(BB(J))
		ENDDO
		CLOSE(UNIT=IN)
	ENDIF
C Find energy in eV from wavelength in Angstoms
	XR=12.397639E3/XRAY
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
	ENDIF
	ABS=ABS/RHO
	END
