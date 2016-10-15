*+SRT_SETT	Set surface quality parameters
	SUBROUTINE SRT_SETT(IS,IT,WL,SR,FM,FI,GD,DHUB,ORDER,ALPHA,GAMMA,
     +	NREFS,ANGS,REFS,ISTAT)
	IMPLICIT NONE
	INTEGER IS,IT,NREFS,ISTAT
	DOUBLE PRECISION WL,SR,FM,FI,GD,DHUB,ORDER,ALPHA,GAMMA
	DOUBLE PRECISION ANGS(*),REFS(*)
*IS	input	surface index
*IT	input	surface type
*WL	input	wavelength
*SR	input	specific roughness A**2 mm
*FM	input	minimum surface frequency mm-1
*FI	input	roughness power law index
*GD	input	grating spacing mm
*DHUB	input	if >1.0 mm then hub distance for off-plane grating
*               if <1.0 mm grating spacing gradient across ruling
*ORDER	input	diffraction order
*ALPHA	input	real part of dielectric constant (or refractive index)
*GAMMA	input	imaginary part of dielectric constant
*NREFS	input	number of angle-reflectivity pairs
*ANGS	input	array of angles
*REFS	input	array of reflectivities
*ISTAT	in/out	returned status
*-Author Dick Willingale 2012-Apr-30
	INCLUDE 'SRT_COM'
	INTEGER NPN,NPI,J,J1,J2
C
	IF(ISTAT.NE.0) RETURN
C Check surface index
	IF(IS.LT.0.OR.IS.GT.MAXST) THEN
		WRITE(*,*) 'SRT_SETT error - surface index out of range'
		ISTAT=1
		RETURN
	ENDIF
C Check parameter space
	NPN=9+NREFS*2
	IF(NPAR+NPN.GT.MAXPAR) THEN
		WRITE(*,*) 'SRT_SETT error - parameter space full;'
		ISTAT=1
		RETURN
	ENDIF
	NPI=NPAR+1
	NPAR=NPAR+NPN
C Set parameters
	ISQP(1,IS)=IT
	ISQP(2,IS)=NPI
	ISQP(3,IS)=NPN
	PAR(NPI)=WL
	PAR(NPI+1)=SR
	PAR(NPI+2)=FM
	PAR(NPI+3)=FI
	PAR(NPI+4)=GD
	PAR(NPI+5)=DHUB
	PAR(NPI+6)=ORDER
	PAR(NPI+7)=ALPHA
	PAR(NPI+8)=GAMMA
	DO J=1,NREFS
      		J2=NPI+8+J*2
      		J1=J2-1
		PAR(J1)=ANGS(J)
		PAR(J2)=REFS(J)
	ENDDO
	END
