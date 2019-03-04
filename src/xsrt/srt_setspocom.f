*+SRT_SETSPOCOM   Set values in constellation common block for SPO array
	SUBROUTINE SRT_SETSPOCOM(NMOD,RM,PM,TM,WM,HM,AM,
     +	CM,GM,RPITCH,WALL,APITCH,RWI,WFR,SIQ,BIQ,ISTAT)
*NMOD	   input   number of MPOs in full array
*RM        input   array of module radial positions (mm)
*PM        input   array of module azimuthal positions (radians)
*TM        input   array of module rotations (radians normally 0)
*WM        input   array of module widths (radial mm)
*HM        input   array of module heights (azimuthal mm)
*AM        input   array of module lengths (axial pore length mm)
*CM        input   array of module curvature indices
*GM        input   array of module grazing angle ratios
*RPITCH    input   array of module pore radial pitch (mm)
*WALL      input   array of module wall thickness (mm)
*APITCH    input   array of module pore azimuthal pitch (rib spacing mm)
*RWI       input   array of module pore rib thickness (mm)
*WFR       input   array of module frame widths (mm)
*SIQ       input   array of module reflecting surface quality indices
*BIQ       input   array of module non-reflecting surface quality indices
*ISTAT	in/out	status 0 OK
      	IMPLICIT NONE
	INTEGER NMOD,ISTAT
	DOUBLE PRECISION RM(NMOD),PM(NMOD),TM(NMOD)
	DOUBLE PRECISION WM(NMOD),HM(NMOD),AM(NMOD)
	DOUBLE PRECISION CM(NMOD),GM(NMOD),RPITCH(NMOD)
	DOUBLE PRECISION WALL(NMOD),APITCH(NMOD),SIQ(NMOD)
	DOUBLE PRECISION WFR(NMOD),RWI(NMOD),BIQ(NMOD)
*-Dick Willingale 2018-Mar-26
	INCLUDE 'SRT_COM'
	INTEGER J
C
	IF(ISTAT.NE.0) RETURN
C
	IF(NMOD.GT.maxlist) THEN
		WRITE(*,*) 'SRT_SETSPOCOM error -# SPOs >',maxlist
		ISTAT=1
		RETURN
	ENDIF
C Set scalars in common
	iptype=1	
	inrad=0.0
	outrad=0.0
	nrofelts=NMOD
C Now run through and set SPO parameters in common
	DO J=1,nrofelts
		rlist(J)=RM(J)
		philist(J)=PM(J)
		thlist(J)=TM(J)
		wlist(J)=WM(J)
		hlist(J)=HM(J)
	        llist(J)=AM(J)
	        clist(J)=CM(J)
	        glist(J)=GM(J)
		olist(J)=RPITCH(J)
		plist(J)=WALL(J)
		qlist(J)=APITCH(J)
		ulist(J)=RWI(J)
		vlist(J)=WFR(J)
		zlist(J)=SIQ(J)
		slist(J)=BIQ(J)
	ENDDO
	END
