*+XX_GET_HENKE	Read Henke et al atomic x-section data from file
	SUBROUTINE XX_GET_HENKE(ID,IO,ATOM,ICOMP,NAT)
	INTEGER ID,IO,NAT,ICOMP(NAT)
	DOUBLE PRECISION DATA,ATWT,ATWTA
	CHARACTER ATOM*(*)
	DIMENSION ATOM(NAT)
*ID	input channel
*IO	output channel : Dummy variable - Ignored.
*ATOM	array of atomic types
*ICOMP	composition
*NAT	number of atomic types
*-Author Dick Willingale 1986-MAR-17
	CHARACTER A*2,XTEMP*80
	CHARACTER*2 AA(100)
	PARAMETER(NMAX=1000)
	PARAMETER(IMAX=10)
	COMMON/HENKE/NATOMA(IMAX),ATWTA(IMAX),NSAMA(IMAX),DATA(IMAX,NMAX,3)
	COMMON/HENKEC/AA
	L=1
C
C Selects atomic data from Henke data file for specified elements.
C
C Formats
1	FORMAT(A2,1X,I4,F8.2,I5)
2	FORMAT(A)
	DO K=1,NAT
		REWIND ID
150		READ(ID,1,IOSTAT=IERR) A,NATOM,ATWT,NSAM
		IF(IERR.NE.0) THEN
			WRITE(6,*) 'Atomic data not available for ',ATOM(K)
			STOP
		ENDIF
		NL=NSAM
		IF(ATOM(K).NE.A) THEN
			DO J=1,NL
				READ(ID,2) XTEMP
			ENDDO
			GOTO 150
		ELSE
			AA(L)=A
			NATOMA(L)=NATOM
			ATWTA(L)=ATWT
			NSAMA(L)=NSAM
			DO J=1,NL
				READ(ID,*) DATA(L,J,1), DATA(L,J,2), DATA(L,J,3)				
			ENDDO
			L=L+1
		ENDIF
	ENDDO
	END
