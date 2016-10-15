*+SRT_SETS	Set source parameters
	SUBROUTINE SRT_SETS(IT,IDEF,SD,SP,SA,AP,AN,AR,PL,ISTAT)
	IMPLICIT NONE
	INTEGER IT,IDEF,ISTAT
	DOUBLE PRECISION SD(3),SP(3),SA,AP(3),AN(3),AR(3),PL(6)
*IT	input	type of source
*IDEF	input	deformation index (used to convert point source into pixels)
*SD	input	direction of source
*SP	input	position of source
*SA	input	aperture area per ray (normal)
*AP	input	reference postion in aperture
*AN	input	normal to aperture
*AR	input	reference axis in aperture
*PL	input	limits and other source parameters
*ISTAT	in/out	returned status
*-Author Dick Willingale 1996-Dec-6
	INCLUDE 'SRT_COM'
	INTEGER NPI,J
C	
	IF(ISTAT.NE.0) RETURN
C Check source type
	IF(IT.LT.0.OR.IT.GT.8) THEN
		WRITE(*,*) 'SRT_SETD error - unknown source type'
		ISTAT=1
		RETURN
	ENDIF
C Check parameter space
	IF(NPAR+22.GT.MAXPAR) THEN
		WRITE(*,*) 'SRT_SETS error - parameter space full'
		ISTAT=1
		RETURN
	ENDIF
	NPI=NPAR+1
	NPAR=NPAR+22
C Set parameters
	ISRC(1)=IT
	ISRC(2)=NPI
	ISRC(3)=IDEF
	DO J=1,3
		PAR(NPI+J-1)=SD(J)
	ENDDO
	DO J=1,3
		PAR(NPI+J+2)=SP(J)
	ENDDO
	PAR(NPI+6)=SA
	DO J=1,3
		PAR(NPI+J+6)=AN(J)
	ENDDO
	DO J=1,3
		PAR(NPI+J+9)=AR(J)
	ENDDO
	DO J=1,3
		PAR(NPI+J+12)=AP(J)
	ENDDO
	DO J=1,6
		PAR(NPI+J+15)=PL(J)
	ENDDO
	END
