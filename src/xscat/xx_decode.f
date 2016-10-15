*+XX_DECODE	Decode material specification string
	SUBROUTINE XX_DECODE(SPEC,MAT,ATOM,ICOMP,NAT,ISTAT)
	CHARACTER SPEC*(*)
	INTEGER MAT,ICOMP(MAT),NAT,ISTAT
	CHARACTER ATOM(MAT)*(2)
        EXTERNAL SCAN
	INTEGER SCAN,LEN_TRIM
*SPEC	input	material specification e.g. 'C8 H8'
*MAT	input	maximum number of atomic types
*ATOM	output	atomic types
*ICOMP	output	number of each atomic type
*NAT	output	number of atomic types
*ISTAT	in/out	returned status
*-Author Dick Willingale 1992-Dec-16
C
	IF(ISTAT.NE.0) RETURN
C
	LS=LEN_TRIM(SPEC)
C
	JS=0
	JE=0
	NAT=0
	DO J=1,LS+1
		IF(J.GT.LS.OR.(JE.EQ.0.AND.(SPEC(J:J).LE.' '))) THEN
			JE=J
			JL=JS+1
			JH=J-1
			NAT=NAT+1
			IF(NAT.GT.MAT) THEN
				WRITE(*,*) 'XX_DECODE too many atomic types'
				ISTAT=1
				RETURN
			ENDIF
			IN=SCAN(SPEC(JL:JH),'0123456789',.FALSE.)
			IF(IN.GT.0) THEN
				IN=IN+JL-1
				IC=JH-IN+1
			ELSE
				IC=0
			ENDIF
			IF(IC.EQ.0) THEN
				ICOMP(NAT)=1
				ATOM(NAT)=SPEC(JL:JH)
				IERR=0
			ELSEIF(IC.EQ.1) THEN
				READ(SPEC(IN:JH),1,IOSTAT=IERR) ICOMP(NAT)
				ATOM(NAT)=SPEC(JL:IN-1)
			ELSEIF(IC.EQ.2) THEN
				READ(SPEC(IN:JH),2,IOSTAT=IERR) ICOMP(NAT)
				ATOM(NAT)=SPEC(JL:IN-1)
			ELSEIF(IC.EQ.3) THEN
				READ(SPEC(IN:JH),3,IOSTAT=IERR) ICOMP(NAT)
				ATOM(NAT)=SPEC(JL:IN-1)
			ELSEIF(IC.EQ.4) THEN
				READ(SPEC(IN:JH),4,IOSTAT=IERR) ICOMP(NAT)
				ATOM(NAT)=SPEC(JL:IN-1)
			ENDIF
			IF(IERR.NE.0) THEN
				WRITE(*,*) 'XX_DECODE failed to read natom'
				ISTAT=1
				RETURN
			ENDIF
    1			FORMAT(I1)
    2			FORMAT(I2)
    3			FORMAT(I3)
    4			FORMAT(I4)
			CALL SYS_UPCASE(ATOM(NAT))
		ELSEIF(JE.NE.0.AND.SPEC(J:J).GT.' ') THEN
			JE=0
			JS=J-1
		ENDIF
	ENDDO
	END
