*+SRT_SU26      Intersection of ray with surface type 26
        SUBROUTINE SRT_SU26(IP,IDEF,ORG,DIR,HIT,POS,RNM,ISTAT)
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
*Type 26 is a cylindrical surface profile
* Parameters in common array PAR(IP)
* (3) axis of cylindrical surface (also x axis on surface)
* (3) tangential reference axis (also y axis on surface)
* (3) centre of cylinder above aperture
* (1) t**2 coefficient=0.0
* (1) t coefficient=0.0
* (1) const coefficient cylinder radius squared 
* (1) minimum x
* (1) minimum y
* (1) maximum x
* (1) maximum y
*-Author Dick Willingale 2010-Apr-27
	INCLUDE 'SRT_COM'
	integer KHIT
C Check status
        IF(ISTAT.NE.0) RETURN
C Use configuration 3 for conic surface radial deformation axial/circum limits
	CALL SRT_CNIC(DIR,ORG,PAR(IP),PAR(IP+3),PAR(IP+6),PAR(IP+9),
     +	0,3,KHIT,POS,RNM,ISTAT)
     	HIT=(KHIT.GT.0.0)
        END
