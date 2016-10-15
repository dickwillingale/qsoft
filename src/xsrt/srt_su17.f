*+SRT_SU17      Intersection of ray with surface type 17
        SUBROUTINE SRT_SU17(IP,IDEF,ORG,DIR,HIT,POS,RNM,ISTAT)
        INTEGER IP,IDEF(2),ISTAT
	LOGICAL HIT
        DOUBLE PRECISION ORG(3),DIR(3),POS(3),RNM(3)
*IP   	input   index of surface parameters
*IDEF	input	deformation indices, IDEF(2)=number of apertures
*ORG    input   origin of ray
*DIR    input   direction of ray
*HIT 	output  hit surface if .true.
*POS    output  position on surface
*RNM    output  normal to surface
*ISTAT  in/out  returned status, 0 is OK
*Type 17 is a plane aperture with normal deformation and azimuthal limits
* Parameters in common array PAR(IP)
* (3) normal to surface
* (3) reference axis on surface
* (3) reference position on surface
* (1) number of sectors
* (1) width of sector arms (cw)
* (1) angular width of sector arms (aw)
* width of arm is cw+aw*r
*-Author Dick Willingale 1996-Nov-10
	INCLUDE 'SRT_COM'
	INTEGER KHIT
C Check status
        IF(ISTAT.NE.0) RETURN
C Use configuration 4
	CALL SRT_PLNA(DIR,ORG,PAR(IP+6),PAR(IP),PAR(IP+3),PAR(IP+9),
     +	IDEF,4,KHIT,POS,RNM,ISTAT)
	HIT=(KHIT.EQ.0)
        END
