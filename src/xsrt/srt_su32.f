*+SRT_SU32      Intersection of ray with surface type 32
        SUBROUTINE SRT_SU32(IP,IDEF,ORG,DIR,HIT,POS,RNM,ISTAT)
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
*Type 32 is an array of plane aperture with normal deformation and cartesian limits
* Parameters in common array PAR(IP) for each aperture in turn
* (3) normal to surface
* (3) reference axis on surface
* (3) reference position on surface
* (1) minimum x
* (1) minimum y
* (1) maximum x
* (1) maximum y
* (1) number of apertures
*-Author Dick Willingale 2019-Mar-03
        INCLUDE 'SRT_COM'
        INTEGER KHIT,NA,I,K
C Check status
        IF(ISTAT.NE.0) RETURN
C Loop for all apertures
        NA=int(PAR(IP+13))
        I=0
        KHIT=0
        DO WHILE((KHIT.EQ.0).AND.(I.LT.NA))
               K=I*14
C Use plane aperture with configuration 1
               CALL SRT_PLNA(DIR,ORG,PAR(IP+6+K),PAR(IP+K),PAR(IP+3+K),
     +         PAR(IP+9+K),IDEF,1,KHIT,POS,RNM,ISTAT)
               I=I+1
        ENDDO
        HIT=(KHIT.EQ.0)
        END
