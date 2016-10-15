*+SRT_MLTI	Calculate reflectivity of a multilayer
	SUBROUTINE SRT_MLTI(ANG,WL,NLAY,FR,D,NPER,RSIG,RPIE,
     +	TSIG,TPIE,ISTAT)
	IMPLICIT NONE
	INTEGER NLAY,NPER,ISTAT
	DOUBLE COMPLEX FR(NLAY)
	DOUBLE PRECISION ANG,WL,D(NLAY),RSIG,RPIE,TSIG,TPIE
*ANG	input	incidence angle (radians from normal in layer 1)
*WL	input	free space wavelength (same units as D)
*NLAY	input	number of layers 
*FR	input	complex refractive index of layers
*D	input	thickness of layers (same units as WL)
*NPER	input	number of periods of internal interfaces
*RSIG	output	reflectivity sigma polarization
*RPIE	output	reflectivity pi polarization
*TSIG	output	transmission sigma polarization
*TPIE	output	transmission pi polarization
*ISTAT	in/out	returned status
*
* Uses characteristic matrix formulation for the layers.
* Principles of Optics, Born and Wolf, 5th Ed. Pergamon Press, 1975
* Crook A.W., J. Opt. Soc. Am., Vol 38, Number 11, 1948
* Lee P., X-ray diffraction from multilayers,
* Optics Communications, Vol 37, Number 3, pages 159-164, 1981
*
* The interfaces 1 to 2 and NLAY-1 to NLAY represent the entrance
* and exit to the multilayer. Therefore the thicknesses of layer 1 and
* layer NLAY are ignored.
*
* Interfaces 2 to 3 through to NLAY-2 to NLAY-1 constitute the internal
* structure of the multilayer. If NPER>1 then layer 2 must be the
* same as layer NLAY-1 and the total number of layers must be odd.
* Then the internal interfaces will be repeated NPER times to construct
* a periodic multilayer.
*
*-Author Dick Willingale 1997-Feb-17
	INCLUDE 'SRT_COM'
	DOUBLE COMPLEX AS(2,2),AP(2,2),AI(2,2),BI(2,2)
	DOUBLE COMPLEX S0,TOS(2,2),TOP(2,2),TRAS,TRAP,SRT_CPII,DR,CR
	EXTERNAL SRT_CPII
	DOUBLE PRECISION WV,CC,DD,PR
	INTEGER J,K,JJ
	LOGICAL TSF,TPF,TSFI,TPFI,TSFR,TPFR
C Check status on input
	IF(ISTAT.NE.0) RETURN
C Check wavelength
	IF(WL.LE.0.0) THEN
		WRITE(*,*) 'SRT_MLTI error - wavelength <= 0'
		ISTAT=1
		RETURN
	ENDIF
C Check incidence angle
	IF(ABS(ANG).GT.PIBY2) THEN
		WRITE(*,*) 'SRT_MLTI error - incidence angle > pi/2'
		ISTAT=1
		RETURN
	ENDIF
	IF(SIN(ANG).EQ.1.0) THEN
		WRITE(*,*) 'SRT_MLTI error - incidence angle = pi/2'
		ISTAT=1
		RETURN
	ENDIF
C Check number of layers
	IF(NLAY.LT.2) THEN
		WRITE(*,*) 'SRT_MLTI error - NLAY must be >1'
		ISTAT=1
		RETURN
	ENDIF
	IF(ABS(NPER).GT.1) THEN
C Check number of layers for periodic structure
	    IF(MOD(NLAY,2).EQ.0) THEN
		WRITE(*,*) 'SRT_MLTI error - NLAY must be odd for NPER>1'
		ISTAT=1
		RETURN
	    ENDIF
	    IF(NLAY.LT.5) THEN
		WRITE(*,*) 'SRT_MLTI error - NLAY too small for NPER>1'
		ISTAT=1
		RETURN
	    ENDIF
C Check refractive indices for periodic structures
	    IF(FR(2).NE.FR(NLAY-1)) THEN
		WRITE(*,*) 'SRT_MLTI error - n(2) <> n(nlay-1) for NPER>1'
		ISTAT=1
		RETURN
	    ENDIF
	ENDIF
C Take sine of angle of incidence in free space
	S0=SIN(ANG)*FR(1)
C Calculate incidence free space wavevector amplitude
	WV=2.0*PI/WL
C Calculate characteristic matrix for entrance interface
	IF(NLAY.GT.2) THEN
		DD=D(2)
	ELSE
		DD=0.0
	ENDIF
	TSF=.TRUE.
	TPF=.TRUE.
	CALL SRT_MLAY(S0,WV,DD,FR(1),FR(2),TOS,TOP,TSF,TPF)
C Initialize characteristic matrices for internal interfaces
	DO J=1,2
		DO K=1,2
			AS(K,J)=CMPLX(0.0,0.0)
			AP(K,J)=CMPLX(0.0,0.0)
		ENDDO
	ENDDO
	AS(1,1)=CMPLX(1.0,0.0)
	AS(2,2)=CMPLX(1.0,0.0)
	AP(1,1)=CMPLX(1.0,0.0)
	AP(2,2)=CMPLX(1.0,0.0)
C Loop for all internal interfaces
	TSFR=TSF
	TPFR=TPF
	DO J=3,NLAY-1
C Calculate characteristic matrix for interface
		DD=D(J)
		CALL SRT_MLAY(S0,WV,DD,FR(J-1),FR(J),AI,BI,TSFR,TPFR)
C Multiply with current product of previous matrices
		CALL SRT_MCM(AS,AI,TSFR)
		CALL SRT_MCM(AP,BI,TPFR)
	ENDDO
	IF(NPER.GT.1) THEN
C Loop for all periods of internal interfaces
		DO J=1,2
			DO K=1,2
				AI(K,J)=AS(K,J)
				BI(K,J)=AP(K,J)
			ENDDO
		ENDDO
		DO J=2,NPER
			CALL SRT_MCM(AS,AI,TSFR)
			CALL SRT_MCM(AP,BI,TPFR)
		ENDDO
	ELSEIF(NPER.LT.-1) THEN
C Use Abeles result to take nth power of internal matrices if periodic
		IF(TSFR) THEN
			TRAS=(AS(1,1)+AS(2,2))*0.5
			DR=SRT_CPII(ABS(NPER)-1,TRAS)
			CR=SRT_CPII(ABS(NPER)-2,TRAS)
			AS(1,1)=AS(1,1)*DR-CR
			AS(1,2)=AS(1,2)*DR
			AS(2,1)=AS(2,1)*DR
			AS(2,2)=AS(2,2)*DR-CR
		ENDIF
		IF(TPFR) THEN
			TRAP=(AP(1,1)+AP(2,2))*0.5
			DR=SRT_CPII(ABS(NPER)-1,TRAP)
			CR=SRT_CPII(ABS(NPER)-2,TRAP)
			AP(1,1)=AP(1,1)*DR-CR
			AP(1,2)=AP(1,2)*DR
			AP(2,1)=AP(2,1)*DR
			AP(2,2)=AP(2,2)*DR-CR
		ENDIF
	ENDIF
C Multiply by top layer
	CALL SRT_MCM(TOS,AS,TSF)
	CALL SRT_MCM(TOP,AP,TPF)
	IF(NLAY.GT.2) THEN
C Calculate characteristic matrix for exit interface
		CALL SRT_MLAY(S0,WV,0.0,FR(NLAY-1),FR(NLAY),AI,BI,TSF,TPF)
C Multiply with current product of previous matrices
		CALL SRT_MCM(TOS,AI,TSF)
		CALL SRT_MCM(TOP,BI,TPF)
	ENDIF
C Calculate ratio of Poynting vectors required for transmission
	CC=REAL(SQRT(FR(1)**2-S0**2))
	DD=REAL(SQRT(FR(NLAY)**2-S0**2))
	IF(CC.NE.0.0) THEN
		PR=DD/CC
	ELSE
		PR=0.0
	ENDIF
C Calculate reflectivities and transmissions
	CC=ABS(TOS(2,1))
	DD=ABS(TOS(1,1))
	IF(DD.NE.0.0) THEN
		RSIG=(CC/DD)**2
		IF(TSF.AND.DD.LT.1.E15.AND.DD.GT.1.E-15) THEN
			TSIG=PR*(1.0/DD)**2
		ELSE
			TSIG=0.0
		ENDIF
	ELSE
		RSIG=1.0
		TSIG=0.0
	ENDIF
	CC=ABS(TOP(2,1))
	DD=ABS(TOP(1,1))
	IF(DD.NE.0.0) THEN
		RPIE=(CC/DD)**2
		IF(TPF.AND.DD.LT.1.E15.AND.DD.GT.1.E-15) THEN
			TPIE=PR*(1.0/DD)**2
		ELSE
			TPIE=0.0
		ENDIF
	ELSE
		RPIE=1.0
		TPIE=0.0
	ENDIF
	END
