*+SRT_SETMPOARR   Set values in constellation common block for MPO array
	SUBROUTINE SRT_SETMPOARR(HWID,NMOD,XP,YP,TC,WC,HC,AC,
     +	CU,MF,PP,WA,SQ,BU,BV,BZ,SP,ISTAT)
*HWID	input	half width of full array aperture
*NMOD	input	number of MPOs in full array
*XP	input	x position of each MPO
*YP	input	y position of each MPO
*TC	input	rotation angle of each MPO
*WC	input	width of each MPO delx 
*HC	input	height of each MPO dely 
*AC	input	axial length of each MPO (thickness)
*CU	input	radius of curvature of each MPO
*MF	input	size of multifibres in MPO
*PP	input	pitch of pores (pore size + wall thickness)
*WA	input	wall thickness between pores
*SQ	input	reflecting surface quality index for pores
*BU	input	bias angle x radians
*BV	input	bias angle y radians
*BZ	input	efficiency of reflection wrt theory (1 reflection in pore)
*SP	input	spare parameter!
*ISTAT	in/out	status 0 OK
      	IMPLICIT NONE
	INTEGER NMOD,ISTAT
	DOUBLE PRECISION HWID
	DOUBLE PRECISION XP(NMOD),YP(NMOD),TC(NMOD)
	DOUBLE PRECISION WC(NMOD),HC(NMOD),AC(NMOD)
	DOUBLE PRECISION CU(NMOD),MF(NMOD),PP(NMOD)
	DOUBLE PRECISION WA(NMOD),SQ(NMOD),SP(NMOD)
	DOUBLE PRECISION BU(NMOD),BV(NMOD),BZ(NMOD)
*-Dick Willingale 2017-Oct-23
	INCLUDE 'SRT_COM'
	INTEGER J
C
	IF(ISTAT.NE.0) RETURN
C
	IF(NMOD.GT.maxlist) THEN
		WRITE(*,*) 'SRT_SETMPOARR error -# MPOs >',maxlist
		ISTAT=1
		RETURN
	ENDIF
C Set scalars in common
	iptype=1	
	inrad=0.0
	outrad=HWID
	nrofelts=NMOD
C Now run through and set MPO parameters in common
	DO J=1,nrofelts
		rlist(J)=XP(J)
		philist(J)=YP(J)
		thlist(J)=TC(J)
		wlist(J)=WC(J)
		hlist(J)=HC(J)
	        llist(J)=AC(J)
	        clist(J)=CU(J)
	        glist(J)=MF(J)
		olist(J)=PP(J)
		plist(J)=WA(J)
		qlist(J)=SQ(J)
		ulist(J)=BU(J)
		vlist(J)=BV(J)
		zlist(J)=BZ(J)
		slist(J)=SP(J)
	ENDDO
	END
