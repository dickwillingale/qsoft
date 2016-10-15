*+SRT_SU21      Intersection of ray with surface type 21
        SUBROUTINE SRT_SU21(IP,IDEF,ORG,DIR,HIT,POS,RNM,ISTAT)
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
*Type 21 is a MCP with square pores and spherical slump
*The packing geometry depends on parameter in common below
* Parameters in common array PAR(IP)
* (3) axis of surface
* (3) tangential reference axis
* (3) centre of sphere
* (1) radius of sphere
* (1) radius of aperture
* (1) pitch x
* (1) pitch y
* (1) rib x
* (1) rib y
* (1) pore length
* (1) local x (returned by SRT_SPHR with IKON=0)
* (1) local y (returned by SRT_SPHR with IKON=0)
* (1) this surface index
* (1) surface quality
* (1) packing geometry
* (1) minimum pore length
* (1) maximum pore length
* This surface type consists of a series of components:
*	a spherical surface which is the front of the plate
*	a plane square aperture at the top of a pore
*	4 plane sides to the pore which reflect
* The pore aperture and sides are set dynamically for each ray.
* The deformation is applied to the normal and reference axes
* in this dynamic set up. IDEF is not passed to the components.
* 4 deformation matrices are expected by SRT_PORE if IDEF(1)>0:
*	1 to shift effective centre of sphere along CAR reference
*	2 to shift effective centre of sphere along NORM X CAR reference
*	3 to rotate pore about local normal
*	4 to distort pore shape
*       5 figure errors on reflecting surfaces
* If IPACK=6 then use aperture index and matrices are 1 dimensional
*-Author Dick Willingale 2004-Feb-13
	INCLUDE 'SRT_COM'
	INTEGER IPACK
	DOUBLE PRECISION COST
C Check status
        IF(ISTAT.NE.0) RETURN
C	IPACK=PAR(IP+20)
C	IF(IPACK.EQ.1.OR.IPACK.EQ.3) THEN
C Use configuration 6 for spherical surface with grid and cartesian limits
C		CALL SRT_SPHR(DIR,ORG,PAR(IP),PAR(IP+3),PAR(IP+6),PAR(IP+9),
C     +		0,5,HIT,POS,RNM,ISTAT)
C	ELSE
C Use configuration 5 for spherical surface with grid and radial limits
C	ENDIF
C Spherical surface with radial limits for all packing types
	CALL SRT_SPHR(DIR,ORG,PAR(IP),PAR(IP+3),PAR(IP+6),PAR(IP+9),
     +	0,5,HIT,POS,RNM,ISTAT)
	IF(.NOT.HIT) THEN
C Set up pore aperture and reflecting walls
		CALL SRT_PORE(DIR,POS,PAR(IP),PAR(IP+3),PAR(IP+6),PAR(IP+9),
     +		IDEF,ISTAT)
	ENDIF
        END
