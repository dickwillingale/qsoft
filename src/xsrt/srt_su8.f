*+SRT_SU8      Intersection of ray with surface type 8
        SUBROUTINE SRT_SU8(IP,IDEF,ORG,DIR,HIT,POS,RNM,ISTAT)
        INTEGER IP,IDEF(2),ISTAT
	LOGICAL HIT
        DOUBLE PRECISION ORG(3),DIR(3),POS(3),RNM(3)
*IP   	input   index of surface parameters
*IDEF	input	deformation indices, IDEF(2)=number of apertures
*ORG    input   origin of ray
*DIR    input   direction of ray
*HIT    output  true if hit
*POS    output  position on surface
*RNM    output  normal to surface
*ISTAT  in/out  returned status, 0 is OK
*Type 8 is a plane with radial deformation and radial limits, nested
* Parameters in common array PAR(IP)
* (3) normal to surface
* (3) reference axis on surface
* (3) reference position on surface
* (1) minimum radius
* (1) maximum radius
* (1) minimum radius
* (1) maximum radius
* ...
*-Author Dick Willingale 1996-Nov-10
	INCLUDE 'SRT_COM'
	INTEGER IKON,KHIT
C Check status
        IF(ISTAT.NE.0) RETURN
C Use plane with configuration -ve
C The number of annuli is given in IDEF(2)
	IKON=-IDEF(2)
	CALL SRT_PLNE(DIR,ORG,PAR(IP+6),PAR(IP),PAR(IP+3),PAR(IP+9),
     +	IDEF,IKON,KHIT,POS,RNM,ISTAT)
	HIT=(KHIT.NE.0)
        END
