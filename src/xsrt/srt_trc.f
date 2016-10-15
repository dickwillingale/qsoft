*+SRT_TRC       Trace rays
        SUBROUTINE SRT_TRC(ISTAT)
	IMPLICIT NONE
        INTEGER ISTAT
*ISTAT  in/out  returned status, 0 is OK
*Author Dick Willingale 1996-Nov-10
	INCLUDE 'SRT_COM'
	INTEGER MAXNPR
        PARAMETER (MAXNPR=30)
        INTEGER NPR,IQU(MAXNPR),J,KSUR,I,KOLD,NFAIL
        DOUBLE PRECISION RNM(3),RPOS(3,MAXNPR),RDIR(3),RQUA(3)
        LOGICAL HIT
C Check status
        IF(ISTAT.NE.0) RETURN
C Trace rays
	REWIND ITRA
	NFAIL=0
	NDET=0
	NPOS=0
	J=0
	DO WHILE(ISTAT.EQ.0.AND.J.LT.NRAYS)
		J=J+1
                CALL SRT_GRAY(0,MAXNPR,NPR,RPOS,RDIR,RQUA,IQU,ISTAT)
		IF(IDEBUG.EQ.1) THEN
			WRITE(*,*) 'debug surface position direction rqua'
			WRITE(*,990) 0,(RPOS(I,NPR),I=1,3),RDIR,RQUA(2)
990			FORMAT(' debug',I5,1X,3F9.3,1X,3F9.5,F9.1)
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
                                IQU(NPR)=ISTY(KSUR)
				IF(ISTY(KSUR).GT.0) THEN
C Reflect, refract or diffract the ray
                                	CALL SRT_REFS(RPOS(1,NPR),
     +					PAR(IPAR(KSUR)),
     +					RDIR,RNM,ISTY(KSUR),RDIR,RQUA,ISTAT)
					KSUR=IHIT(KSUR)
				ELSEIF(ISTY(KSUR).LT.0) THEN
C Ray detected, ray stops
					NDET=NDET+1
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
				NFAIL=NFAIL+1
                                KSUR=0
                        ENDIF
                ENDDO
                CALL SRT_PRAY(ITRA,MAXNPR,NPR,RPOS,RDIR,RQUA,IQU,ISTAT)
		NPOS=NPOS+NPR
		IF(IDEBUG.EQ.1.AND.MOD(J,1000).EQ.0) THEN
			WRITE(*,*) 'srt_trc - traced ',J
		ENDIF
        ENDDO
        IF(NFAIL.GT.0) THEN
	        WRITE(*,*) 'srt_trc - ',NFAIL,' rays hit max ',MAXNPR
	ENDIF
        END
