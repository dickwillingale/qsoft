*+QRT_OPGRAT Set up single off-plane grating
        SUBROUTINE QRT_OPGRAT(ID,DEFI,IQ,AL,FPOS,ZPOS,GPOS,NLIM,ALIM,
     +  ADIR,RDIR,HPOS,GRAZ,GAM,DHUB)
*ID     input        1 radial, 2 nested rad., 4 cart., 6 slats
*DEFI   input        deformation index
*IQ     input        surface quality index
*AL     input         azimuthal exit angle of zeroth order on cone degrees
*FPOS   input        position of primary focus (not reflected)
*ZPOS   input        position of zeroth order focus (reflected)
*GPOS   input        position of centre of grating
*NLIM   input   number of limits
*ALIM   input        limits on surface
*ADIR   output  grating normal
*RDIR   output  grating ruling
*HPOS   output  grating hub
*GRAZ   output  grating grazing angle radians
*GAM    output  grating cone angle radians
*DHUB   output  grating hub distance
Cf2py  intent(in) ID,DEFI,IQ,AL,FPOS,ZPOS,GPOS,NLIM,ALIM
Cf2py  intent(out) ADIR,RDIR,HPOS,GRAZ,GAM,DHUB
        IMPLICIT NONE
        INTEGER ID, DEFI, IQ, NLIM
        DOUBLE PRECISION AL,FPOS(3),ZPOS(3),GPOS(3),ALIM(6)
        DOUBLE PRECISION ADIR(3),RDIR(3),HPOS(3),GRAZ,GAM,DHUB
*-Author Dick Willingale 2010-Apr-29
        INCLUDE 'QR_COM'
        INTEGER IT,IDEF(2),NP,KSUR,J
        DOUBLE PRECISION SDIR(3),ZDIR(3),COSTG,BDIR(3),CDIR(3),ZL
        DOUBLE PRECISION ALPHA,THETA,PI
        DOUBLE PRECISION PL(16)
C Check status
        IF(ISTAT.NE.0) RETURN
C
        PI=ASIN(1.D0)*2.D0
C
        IDEF(1)=DEFI
C Set grating parameters
C convert angles to radians
        ALPHA=AL*PI/180.0
C Find direction SDIR from grating to primary focus
C Find direction ZDIR from grating to zeroth order focus
C ZL is distance to focus
        ZL=0.0
        do j=1,3
                SDIR(j)=FPOS(j)-GPOS(j)
                ZL=ZL+SDIR(J)**2
                ZDIR(j)=ZPOS(j)-GPOS(j)
        enddo
        ZL=SQRT(ZL)
        CALL SRT_VNRM(ZDIR,ISTAT)
        CALL SRT_VNRM(SDIR,ISTAT)
C Find ADIR grating surface normal
        do j=1,3
                ADIR(j)=ZPOS(j)-FPOS(j)
        enddo
        CALL SRT_VNRM(ADIR,ISTAT)
C find grazing angle
        CALL SRT_VDOT(ZDIR,SDIR,COSTG)
        GRAZ=ACOS(COSTG)*0.5
C axis perpenicular to ADIR outof plane
        CALL SRT_VCRS(ZDIR,ADIR,BDIR)
        CALL SRT_VNRM(BDIR,ISTAT)
C other axis in grating plane surface
        CALL SRT_VCRS(ADIR,BDIR,CDIR)
        CALL SRT_VNRM(CDIR,ISTAT)
C Cone angle gamma
        GAM=asin(sin(GRAZ)/COS(ALPHA))
C 3rd axis in plane
        CALL SRT_VCRS(ADIR,BDIR,CDIR)
        CALL SRT_VNRM(CDIR,ISTAT)
C Angle on plane between CDIR and ruling direction RDIR 
        THETA=ACOS(COS(GAM)/COS(GRAZ))
C Find ruling direction +ve from hub to grating
        do j=1,3
                RDIR(j)=-CDIR(J)*COS(THETA)+BDIR(j)*SIN(THETA)
        enddo
        CALL SRT_VNRM(RDIR,ISTAT)
C Find hub position
        DHUB=ZL*COS(GAM)
        do j=1,3
                HPOS(J)=GPOS(j)-RDIR(j)*DHUB
        enddo
C set normal, reference axis in surface and position of grating
        DO J=1,3
                PL(J)=ADIR(J)
                PL(J+3)=RDIR(J)
                PL(J+6)=GPOS(J)
        ENDDO
C set aperture limits
        DO J=1,NLIM
                PL(J+9)=ALIM(J)
        ENDDO
        NP=9+NLIM
        IF(ID.EQ.1) THEN
C Single annulus, radial limits
                IF(PL(11).GT.PL(12)) THEN
                   WRITE(*,*) 'QRT_OPGRAT error -  bad radial limits'
                   ISTAT=1
                   RETURN
                ENDIF
                IT=7
                IDEF(2)=1
        ELSEIF(ID.EQ.2) THEN
C Nested annuli, radial limits
                IF(NLIM.EQ.0.OR.MOD(NLIM,2).NE.1) THEN
                   WRITE(*,*) 'QRT_OPGRAT error -  bad radial limits'
                   ISTAT=1
                   RETURN
                ENDIF
                IT=8
                CALL SRT_NSUR(KSUR,ISTAT)
                PL(NP+1)=KSUR
                PL(NP+2)=1
                NP=NP+2
                IDEF(2)=(NP-1)/2
        ELSEIF(ID.EQ.3.OR.ID.EQ.4) THEN
C Cartesian limits
                IF(NLIM.NE.4.OR.PL(10).GT.PL(12).OR.
     +          PL(11).GT.PL(13)) THEN
                  WRITE(*,*) 'QRT_OPGRAT error -  bad cartesian limits'
                  ISTAT=1
                  RETURN
                ENDIF
                IF(ID.EQ.3) THEN
                        IT=5
                ELSE
                        IT=5
                ENDIF
                IDEF(2)=0
        ELSEIF(ID.EQ.5) THEN
C Cartesian grid limits pitchx pitchy ribx riby
                IF(NLIM.NE.4.OR.PL(10).LT.PL(12).OR.
     +          PL(11).LT.PL(13)) THEN
                 WRITE(*,*)'QRT_OPGRAT error - bad cart. grid limits'
                 ISTAT=1
                 RETURN
                ENDIF
C surface type 18
                IT=18
                IDEF(2)=0
        ELSEIF(ID.EQ.6) THEN
C Cartesian slats limits as above plus pitchx ribx
                IF(NLIM.NE.6.OR.PL(10).GT.PL(12).OR.
     +          PL(11).GT.PL(13)) THEN
                   WRITE(*,*)'QRT_OPGRAT error - bad cart. slat limits'
                   ISTAT=1
                   RETURN
                ENDIF
C surface type 20
                IT=20
                IDEF(2)=0
        ENDIF
C Finally set parameters in common for grating
        CALL SRT_SETF(0,IT,NP,PL,IDEF,IQ,-1,-1,ISTAT)
        END
