*+SRT_SU32      Intersection of ray with surface type 32
* The plane with specific deformations of Schmidt lobster
* derived from SRT_SU5
        SUBROUTINE SRT_SU32(IP,IDEF,ORG,DIR,HIT,POS,RNM,ISTAT)
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
*Type 5 is a plane with normal deformation and cartesian limits
* Parameters in common array PAR(IP)
* (3) normal to surface
* (3) reference axis on surface
* (3) reference position on surface
* (1) minimum x
* (1) minimum y
* (1) maximum x
* (1) maximum y
* Derived from srt_su5 by Vladimir Tichy
	INCLUDE 'SRT_COM'
	INTEGER KHIT
C Check status
        IF(ISTAT.NE.0) RETURN
C Use plane with specific deformations of SLE
	CALL SRT_PLNS(DIR,ORG,PAR(IP+6),PAR(IP),PAR(IP+3),PAR(IP+9),
     +	IDEF,1,KHIT,POS,RNM,ISTAT)
	HIT=(KHIT.GT.0)
        END
