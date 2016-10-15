*+SRT_SU11      Intersection of ray with surface type 11
        SUBROUTINE SRT_SU11(IP,IDEF,ORG,DIR,HIT,POS,RNM,ISTAT)
        INTEGER IP,IDEF(2),ISTAT
        DOUBLE PRECISION ORG(3),DIR(3),POS(3),RNM(3)
        LOGICAL HIT
*IP   	input   index of surface parameters
*IDEF	input	deformation indices
*ORG    input   origin of ray
*DIR    input   direction of ray
*HIT    output  true if hit surface
*POS    output  position on surface
*RNM    output  normal to surface
*ISTAT  in/out  returned status, 0 is OK
*Type 11 is a conic at normal incidence with axial deformation
* and cartesian limits
* Parameters in common array PAR(IP)
* (3) axis of surface
* (3) radial reference axis
* (3) vertex of conic 
* (3) quadratic coefficients ax**2+bx*c
* (1) minimum x
* ()  minimum y
* (1) maximum x
* (1) maximum y
*-Author Dick Willingale 1996-Nov-10
	INCLUDE 'SRT_COM'
C Check status
        IF(ISTAT.NE.0) RETURN
C Use configuration 2
	CALL SRT_CNIN(DIR,ORG,PAR(IP),PAR(IP+3),PAR(IP+6),PAR(IP+9),
     +	IDEF,2,HIT,POS,RNM,ISTAT)
        END
