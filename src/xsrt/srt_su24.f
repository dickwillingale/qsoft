*+SRT_SU24      Intersection of ray with surface type 24
        SUBROUTINE SRT_SU24(IP,IDEF,ORG,DIR,HIT,POS,RNM,ISTAT)
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
*Type 24 is a plane aperture with normal deformation and Si pores
* Parameters in common array PAR(IP)
* (3) reference position on surface (centre of aperture above join plane)
* (3) normal to surface
* (3) reference axis on surface
* (1) minimum radius
* (1) maximum radius
* (1) focal length
* (1) radial pitch of pores
* (1) azimuthal pitch of pores
* (1) wall thickness of pores
* (1) place marker for x local (value set by SRT_PLNA)
* (1) place marker for y local (value set by SRT_PLNA)
* (1) surface quality
* (1) index of surface
*-Author Dick Willingale 2008-Oct-2
	INCLUDE 'SRT_COM'
C Check status
        IF(ISTAT.NE.0) RETURN
C Use plane aperture with configuration 9
	CALL SRT_PLNA(DIR,ORG,PAR(IP),PAR(IP+3),PAR(IP+6),PAR(IP+9),
     +	IDEF,9,HIT,POS,RNM,ISTAT)
	IF(.NOT.HIT) THEN
C Set up Si pore aperture and reflecting walls
		CALL SRT_SIPORE(DIR,POS,PAR(IP),PAR(IP+3),PAR(IP+6),PAR(IP+9),
     +		IDEF,ISTAT)
	ENDIF
        END
