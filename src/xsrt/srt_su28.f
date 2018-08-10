*+SRT_SU28      Intersection of ray with surface type 28
        SUBROUTINE SRT_SU28(IP,IDEF,ORG,DIR,HIT,POS,RNM,ISTAT)
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
*Type 28 is updated square pore MPO using the code for MCP with square pores
* and spherical slump
*The parameter sequence in common array PAR(IP) has been updated
* (3) axis of surface
* (3) tangential reference axis
* (3) centre of sphere
* (1) radius of sphere
* (1) full aperture half width of array of MPOs
* (1) void (spare)
* (1) void (spare)
* (1) void (spare)
* (1) void (spare)
* (1) void (spare)
* (1) local x (returned by SRT_SPHR with IKON=0)
* (1) local y (returned by SRT_SPHR with IKON=0)
* (1) this surface index
* This surface type consists of a series of components:
*	a spherical surface which is the front of the plate
*	a plane square aperture at the top of a pore
*	4 plane sides to the pore which reflect
* The pore aperture and sides are set dynamically for each ray.
* The deformation components are specified for each MPO in the array
* using the deformation arrays pointed to by indices IDEF.
*-Author Dick Willingale 2017-Oct-23
	INCLUDE 'SRT_COM'
	INTEGER IPACK
	DOUBLE PRECISION COST
C Check status
        IF(ISTAT.NE.0) RETURN
C Spherical surface with cartesian limits for all packing types
	CALL SRT_SPHR(DIR,ORG,PAR(IP),PAR(IP+3),PAR(IP+6),PAR(IP+9),
     +	0,5,HIT,POS,RNM,ISTAT)
	IF(.NOT.HIT) THEN
C Set up pore aperture and reflecting walls
		CALL SRT_SQMPOARR(DIR,POS,PAR(IP),PAR(IP+3),PAR(IP+6),
     +          PAR(IP+9),IDEF,ISTAT)
	ENDIF
        END
