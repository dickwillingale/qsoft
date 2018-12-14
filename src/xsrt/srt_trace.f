*+SRT_TRACE	trace and detect rays and write results to ASCII table files
	SUBROUTINE SRT_TRACE(RIRIS,IOPT,AREA,DSHFT,YBAR,ZBAR,RMS,ISTAT)
*RIRIS	input	radius of iris for blur calculation
*		if 0.0 then no analysis performed
*IOPT	input	-2 save traced.dat and detected.dat files
*		-1 save detected.dat
*		0 don't save files or adjust focus
*		1 adjust focus and save detected.dat
*		2 adjust focus and save detected.dat and traced
*		   only rays with IOPT reflections are used for adjustment
*AREA	output	total detected area (0 if no detected rays)
*DSHFT 	output	detector axial shift to minimize rms (0.0 if IOPT<=0)
*YBAR	output	mean detector y position
*ZBAR	output	mean detector z position
*RMS	output	root mean square radius of detected positions
*ISTAT	in/out	returned status
	IMPLICIT NONE
	INTEGER IOPT,ISTAT
	DOUBLE PRECISION RIRIS,DSHFT,YBAR,ZBAR,RMS,AREA
* Added SU30 and SU31 21st November 2018
*-Author Dick Willingale 2012-May-1
	INCLUDE 'SRT_COM'
	INTEGER MAXNPR
        PARAMETER (MAXNPR=20)
        INTEGER NPR,IQU(MAXNPR),KSUR,KOLD,IDTYPE
        DOUBLE PRECISION RNM(3),RPOS(3,MAXNPR),RDIR(3),RQUA(3),RSPHR
	DOUBLE PRECISION ARAY(MAXNPR)
        LOGICAL HIT
C
        INTEGER IRU,IDU,ITU
	DOUBLE PRECISION DNML(3),DRXL(3),DPOS(3),YAX(3)
	DOUBLE PRECISION X,Y,Z,XR,YR,ZR,XD,YD,ZD,XC,YC,ZC,WL,AD
	INTEGER IR,IB,ID,IT,I,J,K,IERR,NTR
	DOUBLE PRECISION DL,DX,DY,DZ,DR(3)
	DOUBLE PRECISION YB,ZB,Y2,Z2,YPB,ZPB,YYB,ZZB,YP2,ZP2,GX,GY,GZ
	DOUBLE PRECISION RR,SS
	DOUBLE PRECISION A
C Check status
        IF(ISTAT.NE.0) RETURN
C Get detector parameters
	DNML(1)=PAR(IPAR(IDET))
	DNML(2)=PAR(IPAR(IDET)+1)
	DNML(3)=PAR(IPAR(IDET)+2)
	DRXL(1)=PAR(IPAR(IDET)+3)
	DRXL(2)=PAR(IPAR(IDET)+4)
	DRXL(3)=PAR(IPAR(IDET)+5)
	DPOS(1)=PAR(IPAR(IDET)+6)
	DPOS(2)=PAR(IPAR(IDET)+7)
	DPOS(3)=PAR(IPAR(IDET)+8)
	IDTYPE=ISURS(IDET)
	IF(IDTYPE.EQ.14.OR.IDTYPE.EQ.15) THEN
C If spherical detector shift position to spherical surface
		RSPHR=PAR(IPAR(IDET)+9)
		DPOS(1)=DPOS(1)+DNML(1)*RSPHR
		DPOS(2)=DPOS(2)+DNML(2)*RSPHR
		DPOS(3)=DPOS(3)+DNML(3)*RSPHR
	ENDIF
C Find third detector axis
	CALL SRT_VCRS(DNML,DRXL,YAX)
	CALL SRT_VNRM(YAX,ISTAT)
C Open ASCII tabulation files as required
	IF(ABS(IOPT).EQ.2) THEN
	 CALL SYS_GETLUN(IRU,ISTAT)
	 OPEN(UNIT=IRU,FILE='traced.dat',STATUS='UNKNOWN',IOSTAT=IERR)
	 WRITE(IRU,'(A)')
     +	'      RXP          RYP          RZP             AREA           IQU'
	ENDIF
	IF(IOPT.NE.0) THEN
	 CALL SYS_GETLUN(IDU,ISTAT)
	 OPEN(UNIT=IDU,FILE='detected.dat',STATUS='UNKNOWN',IOSTAT=IERR)
	 WRITE(IDU,'(A,A,A,A)')
     +   '         XD            YD            ZD',
     +   '         XC            YC            ZC',
     +   '         XR            YR            ZR',
     +   '        YDET         ZDET          AREA     IREF'
	 IF(IOPT.GT.0) THEN
		CALL SYS_GETLUN(ITU,ISTAT)
	 	OPEN(UNIT=ITU,FILE='temp.dat',STATUS='UNKNOWN',IOSTAT=IERR)
	 ENDIF
	ENDIF
1001	FORMAT(4(1X,F13.5),1X,I10)
1002	FORMAT(3(1X,F13.5),3(1X,F13.9),3(1X,F13.5),2(1X,F13.5),1X,F20.10,1X,I2)
C Initialize sums etc.
	YB=0.0
	ZB=0.0
	Y2=0.0
	Z2=0.0
	YPB=0.0
	ZPB=0.0
	YYB=0.0
	ZZB=0.0
	YP2=0.0
	ZP2=0.0
	SS=0.0
C Trace rays
	NDET=0
	NPOS=0
	NTR=0
	DO WHILE(ISTAT.EQ.0.AND.NTR.LT.NRAYS)
                CALL SRT_GRAY(0,MAXNPR,NPR,RPOS,RDIR,RQUA,IQU,ISTAT)
		ARAY(NPR)=RQUA(2)
		IF(IDEBUG.EQ.1) THEN
			WRITE(*,*) 'debug surface position direction rqua'
			WRITE(*,990) 0,(RPOS(I,NPR),I=1,3),RDIR,RQUA(2)
990			FORMAT(' debug',I5,1X,3F9.3,1X,3F9.5,1X,F9.1)
991			format(' debug',I5,'  missed')
		ENDIF
                KSUR=1
                DO WHILE(ISTAT.EQ.0.AND.KSUR.GT.0)
C Find intersection of ray with surface
                        IF(ISURS(KSUR).EQ.1) THEN
                                CALL SRT_SU1(IPAR(KSUR),ISDF(1,KSUR),
     +				RPOS(1,NPR),RDIR,HIT,RPOS(1,NPR+1),RNM,ISTAT)
                        ELSEIF(ISURS(KSUR).EQ.2) THEN
                                CALL SRT_SU2(IPAR(KSUR),ISDF(1,KSUR),
     +				RPOS(1,NPR),RDIR,HIT,RPOS(1,NPR+1),RNM,ISTAT)
                        ELSEIF(ISURS(KSUR).EQ.3) THEN
                                CALL SRT_SU3(IPAR(KSUR),ISDF(1,KSUR),
     +				RPOS(1,NPR),RDIR,HIT,RPOS(1,NPR+1),RNM,ISTAT)
                        ELSEIF(ISURS(KSUR).EQ.4) THEN
                                CALL SRT_SU4(IPAR(KSUR),ISDF(1,KSUR),
     +				RPOS(1,NPR),RDIR,HIT,RPOS(1,NPR+1),RNM,ISTAT)
                        ELSEIF(ISURS(KSUR).EQ.5) THEN
                                CALL SRT_SU5(IPAR(KSUR),ISDF(1,KSUR),
     +				RPOS(1,NPR),RDIR,HIT,RPOS(1,NPR+1),RNM,ISTAT)
                        ELSEIF(ISURS(KSUR).EQ.6) THEN
                                CALL SRT_SU6(IPAR(KSUR),ISDF(1,KSUR),
     +				RPOS(1,NPR),RDIR,HIT,RPOS(1,NPR+1),RNM,ISTAT)
                        ELSEIF(ISURS(KSUR).EQ.7) THEN
                                CALL SRT_SU7(IPAR(KSUR),ISDF(1,KSUR),
     +				RPOS(1,NPR),RDIR,HIT,RPOS(1,NPR+1),RNM,ISTAT)
                        ELSEIF(ISURS(KSUR).EQ.8) THEN
                                CALL SRT_SU8(IPAR(KSUR),ISDF(1,KSUR),
     +				RPOS(1,NPR),RDIR,HIT,RPOS(1,NPR+1),RNM,ISTAT)
                        ELSEIF(ISURS(KSUR).EQ.9) THEN
                                CALL SRT_SU9(IPAR(KSUR),ISDF(1,KSUR),
     +				RPOS(1,NPR),RDIR,HIT,RPOS(1,NPR+1),RNM,ISTAT)
                        ELSEIF(ISURS(KSUR).EQ.10) THEN
                                CALL SRT_SU10(IPAR(KSUR),ISDF(1,KSUR),
     +				RPOS(1,NPR),RDIR,HIT,RPOS(1,NPR+1),RNM,ISTAT)
                        ELSEIF(ISURS(KSUR).EQ.11) THEN
                                CALL SRT_SU11(IPAR(KSUR),ISDF(1,KSUR),
     +				RPOS(1,NPR),RDIR,HIT,RPOS(1,NPR+1),RNM,ISTAT)
                        ELSEIF(ISURS(KSUR).EQ.12) THEN
                                CALL SRT_SU12(IPAR(KSUR),ISDF(1,KSUR),
     +				RPOS(1,NPR),RDIR,HIT,RPOS(1,NPR+1),RNM,ISTAT)
                        ELSEIF(ISURS(KSUR).EQ.13) THEN
                                CALL SRT_SU13(IPAR(KSUR),ISDF(1,KSUR),
     +				RPOS(1,NPR),RDIR,HIT,RPOS(1,NPR+1),RNM,ISTAT)
                        ELSEIF(ISURS(KSUR).EQ.14) THEN
                                CALL SRT_SU14(IPAR(KSUR),ISDF(1,KSUR),
     +				RPOS(1,NPR),RDIR,HIT,RPOS(1,NPR+1),RNM,ISTAT)
                        ELSEIF(ISURS(KSUR).EQ.15) THEN
                                CALL SRT_SU15(IPAR(KSUR),ISDF(1,KSUR),
     +				RPOS(1,NPR),RDIR,HIT,RPOS(1,NPR+1),RNM,ISTAT)
                        ELSEIF(ISURS(KSUR).EQ.16) THEN
                                CALL SRT_SU16(IPAR(KSUR),ISDF(1,KSUR),
     +				RPOS(1,NPR),RDIR,HIT,RPOS(1,NPR+1),RNM,ISTAT)
                        ELSEIF(ISURS(KSUR).EQ.17) THEN
                                CALL SRT_SU17(IPAR(KSUR),ISDF(1,KSUR),
     +				RPOS(1,NPR),RDIR,HIT,RPOS(1,NPR+1),RNM,ISTAT)
                        ELSEIF(ISURS(KSUR).EQ.18) THEN
                                CALL SRT_SU18(IPAR(KSUR),ISDF(1,KSUR),
     +				RPOS(1,NPR),RDIR,HIT,RPOS(1,NPR+1),RNM,ISTAT)
                        ELSEIF(ISURS(KSUR).EQ.19) THEN
                                CALL SRT_SU19(IPAR(KSUR),ISDF(1,KSUR),
     +				RPOS(1,NPR),RDIR,HIT,RPOS(1,NPR+1),RNM,ISTAT)
                        ELSEIF(ISURS(KSUR).EQ.20) THEN
                                CALL SRT_SU20(IPAR(KSUR),ISDF(1,KSUR),
     +				RPOS(1,NPR),RDIR,HIT,RPOS(1,NPR+1),RNM,ISTAT)
                        ELSEIF(ISURS(KSUR).EQ.21) THEN
                                CALL SRT_SU21(IPAR(KSUR),ISDF(1,KSUR),
     +				RPOS(1,NPR),RDIR,HIT,RPOS(1,NPR+1),RNM,ISTAT)
                        ELSEIF(ISURS(KSUR).EQ.22) THEN
                                CALL SRT_SU22(IPAR(KSUR),ISDF(1,KSUR),
     +				RPOS(1,NPR),RDIR,HIT,RPOS(1,NPR+1),RNM,ISTAT)
                        ELSEIF(ISURS(KSUR).EQ.23) THEN
                                CALL SRT_SU23(IPAR(KSUR),ISDF(1,KSUR),
     +				RPOS(1,NPR),RDIR,HIT,RPOS(1,NPR+1),RNM,ISTAT)
                        ELSEIF(ISURS(KSUR).EQ.24) THEN
                                CALL SRT_SU24(IPAR(KSUR),ISDF(1,KSUR),
     +				RPOS(1,NPR),RDIR,HIT,RPOS(1,NPR+1),RNM,ISTAT)
                        ELSEIF(ISURS(KSUR).EQ.25) THEN
                                CALL SRT_SU25(IPAR(KSUR),ISDF(1,KSUR),
     +				RPOS(1,NPR),RDIR,HIT,RPOS(1,NPR+1),RNM,ISTAT)
                        ELSEIF(ISURS(KSUR).EQ.26) THEN
                                CALL SRT_SU26(IPAR(KSUR),ISDF(1,KSUR),
     +				RPOS(1,NPR),RDIR,HIT,RPOS(1,NPR+1),RNM,ISTAT)
                        ELSEIF(ISURS(KSUR).EQ.27) THEN
                                CALL SRT_SU27(IPAR(KSUR),ISDF(1,KSUR),
     +				RPOS(1,NPR),RDIR,HIT,RPOS(1,NPR+1),RNM,ISTAT)
                        ELSEIF(ISURS(KSUR).EQ.28) THEN
                                CALL SRT_SU28(IPAR(KSUR),ISDF(1,KSUR),
     +				RPOS(1,NPR),RDIR,HIT,RPOS(1,NPR+1),RNM,ISTAT)
                        ELSEIF(ISURS(KSUR).EQ.29) THEN
                                CALL SRT_SU29(IPAR(KSUR),ISDF(1,KSUR),
     +				RPOS(1,NPR),RDIR,HIT,RPOS(1,NPR+1),RNM,ISTAT)
                        ELSEIF(ISURS(KSUR).EQ.30) THEN
                                CALL SRT_SU30(IPAR(KSUR),ISDF(1,KSUR),
     +				RPOS(1,NPR),RDIR,HIT,RPOS(1,NPR+1),RNM,ISTAT)
                        ELSEIF(ISURS(KSUR).EQ.31) THEN
                                CALL SRT_SU31(IPAR(KSUR),ISDF(1,KSUR),
     +				RPOS(1,NPR),RDIR,HIT,RPOS(1,NPR+1),RNM,ISTAT)
			ELSEIF(ISURS(KSUR).EQ.0) THEN
C Null surface
				HIT=.FALSE.
			ELSE
C Unknown surface
				WRITE(*,*) 'SRT_TRC error - unknown surface'
				ISTAT=1
				HIT=.FALSE.
                        ENDIF
			KOLD=KSUR
                        IF(HIT) THEN
C Update position count and history of ray
                                NPR=NPR+1
				ARAY(NPR)=RQUA(2)
                                IQU(NPR)=ISTY(KSUR)
				IF(ISTY(KSUR).GT.0) THEN
C Reflect, refract or diffract the ray
                                	CALL SRT_REFS(RPOS(1,NPR),
     +					PAR(IPAR(KSUR)),
     +					RDIR,RNM,ISTY(KSUR),RDIR,RQUA,ISTAT)
					KSUR=IHIT(KSUR)
					ARAY(NPR)=RQUA(2)
				ELSEIF(ISTY(KSUR).LT.0) THEN
C Ray detected, ray stops
					KSUR=0
				ELSE
C Ray absorbed, ray stops
					KSUR=0
				ENDIF
                        ELSE
C Decide where to go next
                                KSUR=IMISS(KSUR)
                        ENDIF
			IF(IDEBUG.EQ.1) THEN
				IF(HIT) THEN
					WRITE(*,990) KOLD,(RPOS(I,NPR),I=1,3),
     +					RDIR,RQUA(2)
				ELSE
					WRITE(*,991) KOLD
				ENDIF
			ENDIF
C Check number of rays positions
                        IF(NPR.EQ.MAXNPR) THEN
C WRITE(*,*) 'srt_trace - ray hit max positions ',MAXNPR
			 KSUR=0
                        ENDIF
                ENDDO
		NPOS=NPOS+NPR
		NTR=NTR+1
		IF(IDEBUG.EQ.1.AND.MOD(NTR,1000).EQ.0) THEN
			WRITE(*,*) 'srt_trace - traced ',NTR
		ENDIF
C Write ray to file if required
		IF(ABS(IOPT).EQ.2) THEN
			DO J=1,NPR
			WRITE(IRU,1001) (RPOS(I,J),I=1,3),ARAY(J),IQU(J)
			ENDDO
		ENDIF
		IF(IQU(NPR).EQ.-1) THEN
C Detected ray
			NDET=NDET+1
			XR=RPOS(1,NPR-1)
			YR=RPOS(2,NPR-1)
			ZR=RPOS(3,NPR-1)
			XD=RPOS(1,NPR)
			YD=RPOS(2,NPR)
			ZD=RPOS(3,NPR)
			XC=RDIR(1)
			YC=RDIR(2)
			ZC=RDIR(3)
			WL=RQUA(1)
			AD=RQUA(2)
C Count up history of ray
			IR=0
			IB=0
			ID=0
			DO I=1,NPR-1
				IF(IQU(I).EQ.1) THEN
					IR=IR+1
				ELSEIF(IQU(I).EQ.2) THEN
					IB=IB+1
				ELSEIF(IQU(I).EQ.3) THEN
					ID=ID+1
				ENDIF
			ENDDO
C Find length along ray
			DL=0.0
			DO I=1,NPR
				DX=RPOS(1,I)-RPOS(1,I-1)
				DY=RPOS(2,I)-RPOS(2,I-1)
				DZ=RPOS(3,I)-RPOS(3,I-1)
				DL=DL+SQRT(DX*DX+DY*DY+DZ*DZ)
			ENDDO
C Find position in detector coordinates
			DR(1)=XD-DPOS(1)
			DR(2)=YD-DPOS(2)
			DR(3)=ZD-DPOS(3)
			CALL SRT_VDOT(DNML,DR,X)
			CALL SRT_VDOT(DRXL,DR,Y)
			CALL SRT_VDOT(YAX,DR,Z)
C Find gradients with respect to x-shift
			DR(1)=XC
			DR(2)=YC
			DR(3)=ZC
			CALL SRT_VDOT(DR,DNML,GX)
			IF(GX.NE.0.0) THEN
				CALL SRT_VDOT(DR,DRXL,GY)
				CALL SRT_VDOT(DR,YAX,GZ)
				GY=GY/GX
				GZ=GZ/GX
			ELSE
				GY=0.0
				GZ=0.0
			ENDIF
C Only rays with IOPT reflections within iris radius count
			IF(IOPT.GT.0) THEN
				IT=IR+IB
			ELSE
				IT=IOPT
			ENDIF
			IF(SQRT(Y**2+Z**2).LT.RIRIS.AND.IT.EQ.IOPT) THEN
				A=AD
				SS=SS+A
				YB=YB+Y*A
				ZB=ZB+Z*A
				Y2=Y2+A*Y**2
				Z2=Z2+A*Z**2
				YPB=YPB+GY*A
				ZPB=YPB+GZ*A
				YYB=YYB+Y*GY*A
				ZZB=ZZB+Z*GZ*A
				YP2=YP2+GY*GY*A
				ZP2=ZP2+GZ*GZ*A
			ENDIF
			IF(IOPT.GT.0) THEN
C Put detected info on temporary file
				WRITE(ITU,*) XD,YD,ZD,XC,YC,ZC,XR,YR,ZR,
     +				X,Y,Z,GX,GY,GZ,AD,IR
			ELSEIF(IOPT.LT.0) THEN
C Put detected info on detected.dat file
				WRITE(IDU,1002) XD,YD,ZD,XC,YC,ZC,
     +				XR,YR,ZR,Y,Z,AD,IR
			ENDIF
		ENDIF
	ENDDO
	IF(ABS(IOPT).EQ.2) THEN
		CLOSE(IRU)
	ENDIF
C Calculate area weighted mean and rms radius of detected positions
	DSHFT=0.0
	IF(SS.GT.0.0) THEN
		YB=YB/SS
		ZB=ZB/SS
		Y2=Y2/SS
		Z2=Z2/SS
		YPB=YPB/SS
		ZPB=ZPB/SS
		YYB=YYB/SS
		ZZB=ZZB/SS
		YP2=YP2/SS
		ZP2=ZP2/SS
		RR=YP2+ZP2-YPB**2-ZPB**2
		IF(IOPT.GT.0.AND.RR.GT.0.0) THEN
C Shift detector along normal to position of minimum rms radius
			DSHFT=(YPB*YB+ZPB*ZB-YYB-ZZB)/RR
C Shift detector along normal to position of minimum rms y
C			DSHFT=(YPB*YB-YYB)/(YP2-YPB**2)
C Shift detector along normal to position of minimum rms z
C			DSHFT=(ZPB*ZB-ZZB)/(ZP2-ZPB**2)
		ENDIF
	ENDIF
C
	AREA=SS
	YBAR=DSHFT*YPB+YB
	ZBAR=DSHFT*ZPB+ZB
	RMS=(YP2+ZP2-YPB**2-ZPB**2)*DSHFT**2+
     +	(YYB+ZZB-YPB*YB-ZPB*ZB)*DSHFT*2.D0+Y2+Z2-YB**2-ZB**2
	RMS=SQRT(ABS(RMS))
C 
	IF(IOPT.GT.0) THEN
		REWIND ITU
		DO J=1,NDET
			READ(ITU,*) XD,YD,ZD,XC,YC,ZC,XR,YR,ZR,
     +			X,Y,Z,GX,GY,GZ,AD,IR
			DR(1)=DSHFT
			DR(2)=DSHFT*GY
			DR(3)=DSHFT*GZ
			Y=Y+DR(2)
			Z=Z+DR(3)
			XD=XD+DR(1)*DNML(1)+DR(2)*DRXL(1)+DR(3)*YAX(1)
			YD=YD+DR(1)*DNML(2)+DR(2)*DRXL(2)+DR(3)*YAX(2)
			ZD=ZD+DR(1)*DNML(3)+DR(2)*DRXL(3)+DR(3)*YAX(3)
			WRITE(IDU,1002) XD,YD,ZD,XC,YC,ZC,
     +			XR,YR,ZR,Y,Z,AD,IR
		ENDDO
		CLOSE(ITU)
	ENDIF
	IF(IOPT.NE.0) THEN
		CLOSE(IDU)
	ENDIF
	END
