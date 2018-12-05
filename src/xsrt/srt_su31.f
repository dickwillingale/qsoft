*+SRT_SU31     Intersection of ray with surface type 31,  derived from SRT_SU25
        SUBROUTINE SRT_SU31(IP,IDEF,ORG,DIR,HIT,POS,RNM,ISTAT)
        INTEGER IP,IDEF(2),ISTAT
        DOUBLE PRECISION ORG(3),DIR(3),POS(3),RNM(3)
        LOGICAL HIT
*IP           input   index of surface parameters
*IDEF        input        deformation indices
*ORG    input   origin of ray
*DIR    input   direction of ray
*HIT    output  true if hit surface
*POS    output  position on surface
*RNM    output  normal to surface
*ISTAT  in/out  returned status, 0 is OK
*Type 31 is a spherical aperture plane and Schmidt Lobster stack
* Parameters in common array PAR(IP)
* (3) axis of surface
* (3) tangential reference axis
* (3) centre of curvature of spherical aperture 
* (1) radius of sphere (2*FLEN+PLMAX)
* (1) maximum radius
* (1) minimum radius
* (1) pitch of slots
* (1) wall thickness of wafers
* (1) minimum axial slot length
* (1) maximum axial slot length
* (1) local x (returned by SRT_SPHR)
* (1) local y (returned by SRT_SPHR)
* (1) index of surface
* (1) surface quality
*-
* Derived from srt_su25 by Vladimir Tichy
*
* last modify 24 Jan 2018
*
        INCLUDE 'SRT_COM'
C Check status
        IF(ISTAT.NE.0) RETURN
C Use configuration 5 for spherical surface with grid and radial limits
        CALL SRT_SPHR(DIR,ORG,PAR(IP),PAR(IP+3),PAR(IP+6),PAR(IP+9),
     +        0,5,HIT,POS,RNM,ISTAT)
C
C this is exexuted if aperture is not hit, i.e. HIT=.FALSE.
        IF(.NOT.HIT) THEN
C                WRITE (*,*) 'ORG:',ORG(1), ORG(2), ORG(3)
C                WRITE (*,*) 'POS:',POS(1), POS(2), POS(3)
C                WRITE (*,*)
C Set up SLE stack aperture and reflecting walls
                CALL SRT_SLE(DIR,ORG,POS,PAR(IP),PAR(IP+3),PAR(IP+6),PAR(IP+9),
     +                IDEF,HIT,ISTAT)
C
C pokus:
                POS(1)=ORG(1)
                POS(2)=ORG(2)
                POS(3)=ORG(3)
        ENDIF
        END
