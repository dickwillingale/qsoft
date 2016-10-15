*+SRT_SU23      Intersection of ray with surface type 23
        SUBROUTINE SRT_SU23(IP,IDEF,ORG,DIR,HIT,POS,RNM,ISTAT)
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
*Type 23 is a MOA with a cylindrical surface profile and rectangular slots
* Parameters in common array PAR(IP)
* (3) axis of cylindrical surface (also x axis on surface of MOA)
* (3) tangential reference axis (also y axis on surface of MOA at centre)
* (3) centre of cylinder below aperture
* (1) t**2 coefficient=0.0
* (1) t coefficient=0.0
* (1) const coefficient cylinder radius squared 
* (1) minimum x
* (1) minimum y
* (1) maximum x
* (1) maximum y
* (1) pitch of slots in y
* (1) rib width between slots in y
* (1) slot depth (thickness of MOA)
* (1) this surface index
* (1) surface quality
* This surface type consists of a series of components:
*	a cylindrical surface which is the front of the MOA
*	a plane rectanglar aperture at the top of each slot
*	4 plane sides for the slot, all reflecting
* The slot aperture and sides are set dynamically for each ray.
* The deformation is applied to the normal and reference axes
* in this dynamic set up. IDEF is not passed to the components.
*-Author Dick Willingale 2007-Feb-01
	INCLUDE 'SRT_COM'
C Check status
        IF(ISTAT.NE.0) RETURN
C Use configuration 3 for conic surface radial deformation axial/circum limits
	CALL SRT_CNIC(DIR,ORG,PAR(IP),PAR(IP+3),PAR(IP+6),PAR(IP+9),
     +	0,3,HIT,POS,RNM,ISTAT)
	IF(HIT) THEN
C Set up slot aperture and reflecting walls
		CALL SRT_MOA(DIR,POS,PAR(IP),PAR(IP+3),PAR(IP+6),PAR(IP+9),
     +		IDEF,ISTAT)
C Note that HIT is returned as compliment of value from SRT_CNIC
		HIT=.false.
	ELSE
		HIT=.true.
	ENDIF
        END
