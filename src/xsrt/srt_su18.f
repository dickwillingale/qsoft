*+SRT_SU18      Intersection of ray with surface type 18
        SUBROUTINE SRT_SU18(IP,IDEF,ORG,DIR,HIT,POS,RNM,ISTAT)
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
*Type 18 is a plane aperture with normal deformation and cartesian grid limits
* Parameters in common array PAR(IP)
* (3) normal to surface
* (3) reference axis on surface
* (3) reference position on surface
* (1) pitch x
* (1) pitch y
* (1) rib x (< pitch x)
* (1) rib y (< pitch y)
*-Author Dick Willingale 2001-Nov-9
	INCLUDE 'SRT_COM'
	INTEGER KHIT
C Check status
        IF(ISTAT.NE.0) RETURN
C Use plane aperture with configuration 6
	CALL SRT_PLNA(DIR,ORG,PAR(IP+6),PAR(IP),PAR(IP+3),PAR(IP+9),
     +	IDEF,6,KHIT,POS,RNM,ISTAT)
	HIT=(KHIT.EQ.0)
        END
