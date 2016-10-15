*+SRT_RELK	Reflectivity from look up table
	SUBROUTINE SRT_RELK(ANG,NP,P,REF)
	IMPLICIT NONE
	INTEGER NP
	DOUBLE PRECISION ANG,P(2,NP),REF
*ANG	input	incidence angle radians, range 0 to pi/2
*NP	input	number of angle-reflectivity pairs in table
*P	input	look-up table, angles in increasing order
*REF	output	reflectivity interpolated from table
*-Author Dick Willingale 1996-Nov-28
	DOUBLE PRECISION FR
	INTEGER J
C Search for incidence angle in table and interpolate
C If angle out of range of table then returns REF=0.0
	REF=0.0
	DO J=1,NP-1
		IF(P(1,J).LE.ANG.AND.P(1,J+1).GT.ANG) THEN
			FR=(ANG-P(1,J))/(P(1,J+1)-P(1,J))
			REF=FR*(P(2,J+1)-P(2,J))+P(2,J)
			RETURN
		ENDIF
	ENDDO
	END
