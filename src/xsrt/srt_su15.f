*+SRT_SU15      Intersection of ray with surface type 15
        SUBROUTINE SRT_SU15(IP,IDEF,ORG,DIR,HIT,POS,RNM,ISTAT)
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
*Type 15 is a sphere at normal incidence with radial deformation
* and radial limits
* Parameters in common array PAR(IP)
* (3) axis of surface
* (3) tangential reference axis
* (3) centre of sphere
* (1) radius
* (1) minimum r
* (1) maximum r
*-Author Dick Willingale 1996-Nov-10
	INCLUDE 'SRT_COM'
C Check status
        IF(ISTAT.NE.0) RETURN
C Use configuration 4
	CALL SRT_SPHR(DIR,ORG,PAR(IP),PAR(IP+3),PAR(IP+6),PAR(IP+9),
     +	IDEF,4,HIT,POS,RNM,ISTAT)
        END
