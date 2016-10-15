*+SRT_SU7      Intersection of ray with surface type 7
        SUBROUTINE SRT_SU7(IP,IDEF,ORG,DIR,HIT,POS,RNM,ISTAT)
        INTEGER IP,IDEF(2),ISTAT
	LOGICAL HIT
        DOUBLE PRECISION ORG(3),DIR(3),POS(3),RNM(3)
*IP   	input   index of surface parameters
*IDEF	input	deformation indices, IDEF(2)=number of apertures
*ORG    input   origin of ray
*DIR    input   direction of ray
*KHIT   output  number of annulus hit (0 if not hit)
*POS    output  position on surface
*RNM    output  normal to surface
*ISTAT  in/out  returned status, 0 is OK
*Type 7 is a plane with radial deformation and radial limits, single 
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
	INTEGER KHIT
C Check status
        IF(ISTAT.NE.0) RETURN
C Use configuration 3
	CALL SRT_PLNE(DIR,ORG,PAR(IP+6),PAR(IP),PAR(IP+3),PAR(IP+9),
     +	IDEF,3,KHIT,POS,RNM,ISTAT)
	HIT=(KHIT.NE.0)
        END
