*+QRT_MIRROR Set up plane mirror
        SUBROUTINE QRT_MIRROR(ID,IDF,IQ,ANML,ARFX,APOS,N,ALIM,NSURF)
        IMPLICIT NONE
        INTEGER ID,IDF,IQ,N,NSURF
        DOUBLE PRECISION ANML(3),ARFX(3),APOS(3),ALIM(N)
*ID        input        aperture type
*IDF       input        deformation index
*IQ        input        surface quality index
*ANML      input        surface normal
*ARFX      input        surface reference axis
*APOS      input        surface reference position
*N         input        number of limits
*ALIM      input        limits array
*NSURF     input        number of subsequent surfaces ID=2
Cf2py  intent(in) ID,IDF,IQ,ANML,ARFX,APOS,N,ALIM,NSURF
*-Author Dick Willingale 2012-Jun-28
        INCLUDE 'SRT_COM'
        INTEGER IT,IDEF(2),NP,MAXRS,ISS,KSUR,I
        PARAMETER (MAXRS=101)
        DOUBLE PRECISION PL(MAXRS+11)
        INCLUDE 'QR_COM'
C
        IF(ISTAT.NE.0) RETURN
C
        IDEF(1)=IDF
C Get surface parameters
        PL(1)=ANML(1)
        PL(2)=ANML(2)
        PL(3)=ANML(3)
        PL(4)=ARFX(1)
        PL(5)=ARFX(2)
        PL(6)=ARFX(3)
        PL(7)=APOS(1)
        PL(8)=APOS(2)
        PL(9)=APOS(3)
        DO I=1,N
                PL(I+9)=ALIM(I)
        ENDDO
C
        IF(ID.EQ.1) THEN
C Single annulus, radial limits
                IF(N.NE.3.OR.PL(11).GT.PL(12)) THEN
                    WRITE(*,*) 'QRT_MIRROR error -  bad radial limits'
                    ISTAT=1
                    RETURN
                ENDIF
                IT=7
                NP=9+N
                IDEF(2)=1
        ELSEIF(ID.EQ.2) THEN
C Nested annuli, radial limits
                IF(N.EQ.0.OR.MOD(N,2).NE.1) THEN
                    WRITE(*,*) 'QRT_MIRROR error -  bad radial limits'
                    ISTAT=1
                    RETURN
                ENDIF
                IT=8
                CALL SRT_NSUR(KSUR,ISTAT)
                ISS=NSURF
                PL(N+10)=KSUR
                PL(N+11)=ISS
                NP=11+N
                IDEF(2)=(N-1)/2
        ELSEIF(ID.EQ.3.OR.ID.EQ.4) THEN
C Cartesian limits
                IF(N.NE.4.OR.PL(10).GT.PL(12).OR.PL(11).GT.PL(13)) THEN
                  WRITE(*,*) 'QRT_MIRROR error -  bad cartesian limits'
                  ISTAT=1
                  RETURN
                ENDIF
                IF(ID.EQ.3) THEN
                        IT=5
                ELSE
                        IT=5
                ENDIF
                NP=9+N
                IDEF(2)=0
        ELSEIF(ID.EQ.5) THEN
C Cartesian grid limits pitchx pitchy ribx riby
                IF(N.NE.4.OR.PL(10).LT.PL(12).OR.PL(11).LT.PL(13)) THEN
                  WRITE(*,*)'QRT_MIRROR error - bad cart. grid limits'
                  ISTAT=1
                  RETURN
                ENDIF
C surface type 18
                IT=18
                NP=9+N
                IDEF(2)=0
        ELSEIF(ID.EQ.7) THEN
C Cartesian slats limits as above plus pitchx ribx
                IF(N.NE.6.OR.PL(10).GT.PL(12).OR.PL(11).GT.PL(13)) THEN
                  WRITE(*,*)'QRT_MIRROR error - bad cart. slat limits'
                  ISTAT=1
                  RETURN
                ENDIF
C surface type 20
                IT=20
                NP=9+N
                IDEF(2)=0
        ELSE 
                WRITE(*,*) 'QRT_MIRROR error - aperture type',ID
                ISTAT=1
                RETURN
        ENDIF
C Finally set parameters in common for mirror
        CALL SRT_SETF(0,IT,NP,PL,IDEF,IQ,-1,-1,ISTAT)
        END
