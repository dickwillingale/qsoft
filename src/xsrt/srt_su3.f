*+SRT_SU3      Intersection of ray with surface type 3
        SUBROUTINE SRT_SU3(IP,IDEF,ORG,DIR,HIT,POS,RNM,ISTAT)
        INTEGER IP,IDEF(2),ISTAT
	LOGICAL HIT
        DOUBLE PRECISION ORG(3),DIR(3),POS(3),RNM(3)
*IP   	input   index of surface parameters
*IDEF	input	deformation indices, IDEF(2)=number of apertures
*ORG    input   origin of ray
*DIR    input   direction of ray
*KHIT   output  aperture number hit (0 if missed aperture)
*POS    output  position on surface
*RNM    output  normal to surface
*ISTAT  in/out  returned status, 0 is OK
*Type 3 is a plane aperture with radial deformation and radial limits, single
* Parameters in common array PAR(IP)
* (3) normal to surface
* (3) reference axis on surface
* (3) reference position on surface
* (1) axial reference position of limits
* (1) minimum radius
* (1) maximum radius
*-Author Dick Willingale 1996-Nov-10
	INCLUDE 'SRT_COM'
C Check status
        IF(ISTAT.NE.0) RETURN
C Use configuration 3
	CALL SRT_PLNA(DIR,ORG,PAR(IP+6),PAR(IP),PAR(IP+3),PAR(IP+9),
     +	IDEF,3,KHIT,POS,RNM,ISTAT)
	HIT=(KHIT.EQ.0)
        END
