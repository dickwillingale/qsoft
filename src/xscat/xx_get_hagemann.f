*+XX_GET_HAGEMANN	Read Hagemann et al atomic x-section data from file
	SUBROUTINE XX_GET_HAGEMANN(ID,IO,ATOM,ICOMP,NAT)
	INTEGER ID,IO,NAT,ICOMP(NAT)
	CHARACTER ATOM*(*)
	DIMENSION ATOM(NAT)
*ID	input channel
*IO	output channel
*ATOM	array of atomic types
*ICOMP	composition
*NAT	number of atomic types
*-Author Dick Willingale 1986-June-2
	CHARACTER A*2,XTEMP*80
C
C Selects atomic data from Hagemann data file for specified elements.
C
C Formats
1	FORMAT(A2,8X,F10.3,F10.3,I10)
2	FORMAT(A)
	DO K=1,NAT
		REWIND ID
150		READ(ID,1,IOSTAT=IERR) A,ATWT,RHO,NSAM
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
			WRITE(IO,1) A,ATWT,RHO,NSAM
			DO J=1,NL
				READ(ID,2) XTEMP
				WRITE(IO,2) XTEMP
			ENDDO
		ENDIF
	ENDDO
	RETURN
	END
