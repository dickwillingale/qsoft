*+SRT_SCAT	Scattering from surface roughness
	SUBROUTINE SRT_SCAT(WL,SPECROUGH,BREAKROUGH,RINDROUGH,GRAZ,NA,
     {	ADEL,ISTAT)
	IMPLICIT NONE
	INTEGER NA,ISTAT
	DOUBLE PRECISION WL,SPECROUGH,BREAKROUGH,RINDROUGH,GRAZ,ADEL(NA)
* WL		input	wavelength A
* SPECROUGH	input	specific roughness A^2 mm
*			if -ve then rms figure gradient error radians
* BREAKROUGH	input	mm^-1
* RINDROUGH	input	power law index
* GRAZ		input	grazing angle radians
* NA		input	number of scattering angles
* ADEL		output	array of scattering angles radians
* ISTAT		in/out	status
*-Author Dick Willingale 2007-Jan-31
	INCLUDE 'SRT_COM'
	INTEGER I
	DOUBLE PRECISION TIS,OMEGA,GAMMA1,RND,SIGMAS,RR,AB
        DOUBLE PRECISION SYS_DRAND
	IF(ISTAT.NE.0) RETURN
C
	IF(SPECROUGH.GT.0.0.AND.WL.GT.0.0) THEN
		GAMMA1=RINDROUGH-1.0D0
C Calculate rms roughness and TIS
		SIGMAS=SPECROUGH*(1.0+1.0/GAMMA1)*(BREAKROUGH**(-GAMMA1))
		TIS=1.0-EXP(-SIGMAS*(4.0*PI*SIN(GRAZ)/WL)**2)
		AB=(BREAKROUGH*(WL*1.D-7)/SIN(GRAZ))*180.0*3600/PI
C		write(*,*) 'wl graz',WL,GRAZ
C		write(*,*) 'rms roughness A',sqrt(sigmas)
C		write(*,*) 'TIS',TIS
C		write(*,*) 'scattering angle arc secs at break',AB
		DO I=1,NA
C Calculate scattering angle
			IF(SYS_DRAND().LE.TIS) THEN
C Scatter
				RND=SYS_DRAND()
				IF(RND.GT.(GAMMA1/RINDROUGH)) THEN
C Above the break frequency in the power spectrum
				    RND=MIN(RND,0.99)
				    RR=RINDROUGH*(1.0-RND)
				    OMEGA=BREAKROUGH*RR**(-1.0/GAMMA1)
				ELSE
				    OMEGA=RND*BREAKROUGH*(RINDROUGH/GAMMA1)
				ENDIF
C OMEGA is spatial frequency - calculate scattering angle using grating eq.
				ADEL(I)=OMEGA*(WL*1.D-7)/SIN(GRAZ)
C set either +ve or negative
				IF(SYS_DRAND().LT.0.5) THEN
					ADEL(I)=-ADEL(I)
				ENDIF
			ELSE
				ADEL(I)=0.0
			ENDIF
		ENDDO
	ELSEIF(SPECROUGH.LT.0.0) THEN
		CALL SYS_GAUSS(NA,ADEL,0.0D0,ABS(SPECROUGH),ISTAT)
	ELSE
		DO I=1,NA
			ADEL(I)=0.0
		ENDDO
	ENDIF
	END
