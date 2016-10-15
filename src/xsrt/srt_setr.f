*+SRT_SETR	Transfer angle-reflectivity samples to a single array
	SUBROUTINE SRT_SETR(NA,ANGS,NR,REFS,NP,P,ISTAT)
	IMPLICIT NONE
	INTEGER NA,NR,NP,ISTAT
	DOUBLE PRECISION ANGS(NA),REFS(NR),P(NP)
*NA	input	number of angles
*ANGS	input	incidence angles (degrees)
*NR	input	number of reflectivities
*REFS	input	reflectvities
*NP	input	number of pairs (angles in radians)
*P	output	single array of pairs
*ISTAT	in/out	returned status
*-Author Dick Willingale 1997-Mar-4
	INCLUDE 'SRT_COM'
	INTEGER J,JJ
C
	IF(ISTAT.NE.0) RETURN
C
	IF(NA.NE.NR.OR.NP.NE.NA+NR) THEN
		WRITE(*,*) 'SRT_SETR error - inconsistent dimensions'
		ISTAT=1
		RETURN
	ENDIF
	DO J=1,NA
		JJ=J*2
		P(JJ-1)=ANGS(J)*PI/180.
		P(JJ)=REFS(J)
	ENDDO
	END
