*+SRT_SU1      Intersection of ray with surface type 1
        SUBROUTINE SRT_SU1(IP,IDEF,ORG,DIR,HIT,POS,RNM,ISTAT)
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
*Type 1 is a plane aperture with normal deformation and cartesian limits
* Parameters in common array PAR(IP)
* (3) normal to surface
* (3) reference axis on surface
* (3) reference position on surface
* (1) minimum x
* (1) minimum y
* (1) maximum x
* (1) maximum y
*-Author Dick Willingale 1996-Nov-10
	INCLUDE 'SRT_COM'
	INTEGER KHIT
C Check status
        IF(ISTAT.NE.0) RETURN
C Use plane aperture with configuration 1
	CALL SRT_PLNA(DIR,ORG,PAR(IP+6),PAR(IP),PAR(IP+3),PAR(IP+9),
     +	IDEF,1,KHIT,POS,RNM,ISTAT)
	HIT=(KHIT.EQ.0)
        END
