*+QRT_APERTURE Set aperture stop parameters
        SUBROUTINE QRT_APERTURE(ID,IDF,AP,AN,AR,NL,ALIM,NSURF)
*ID        input        aperture type
*   1 Single annulus, radial limits (aref,rmin1,rmax1,rmin2...)
*   2 Nested annuli, radial limits        (aref,rmin1,rmax1,rmin2...)
*   3 Hole Cartesian limits (xmin,ymin,xmax,ymax)
*   4 Block Cartesian limits (xmin,ymin,xmax,ymax)
*   5 Cartesian grid limits (pitchx pitchy ribx riby)
*   6 Radial/azimuthal sector limits (rmin,rmax,amin,amax)
*     amin and amax in radians range 0-2pi
*   7 Parallelogram limits (xmin,ymin,xmax,ymax,dx)
*   8 Aperture for MCO test station limits (hsize dols)
*IDF      input        deformation index
*AP       input        position of aperture
*AN       input        normal to aperture plane
*AR       input        reference axis in aperture plane
*NL       input        number of limit values
*ALIM     input        limit values
*NSURF    input        number of subsequent surfaces per aperture (ID=2)
Cf2py  intent(in) ID,IDF,AP,AN,AR,NL,ALIM,NSURF
        IMPLICIT NONE
        INTEGER ID,IDF,NL,NSURF
        DOUBLE PRECISION AP(3),AN(3),AR(3),ALIM(NL)
*-Author Dick Willingale 2012-May-16
        INTEGER IT,IDEF(2),NP,MAXRS,KSUR,J
        PARAMETER (MAXRS=101)
        DOUBLE PRECISION PL(MAXRS+11)
        INCLUDE 'QR_COM'
C
        IF(ISTAT.NE.0) RETURN
C
        IDEF(1)=IDF
C Set surface parameters
        PL(1)=AN(1)
        PL(2)=AN(2)
        PL(3)=AN(3)
        PL(4)=AR(1)
        PL(5)=AR(2)
        PL(6)=AR(3)
        PL(7)=AP(1)
        PL(8)=AP(2)
        PL(9)=AP(3)
        IF(ID.EQ.1) THEN
C Single annulus, radial limits
                IF(NL.NE.3.OR.ALIM(2).GT.ALIM(3)) THEN
                   WRITE(*,*) 'QRT_APERTURE error -  bad radial limits'
                   ISTAT=1
                   RETURN
                ENDIF
                IT=3
                NP=12
                PL(10)=ALIM(1)
                PL(11)=ALIM(2)
                PL(12)=ALIM(3)
                IDEF(2)=1
        ELSEIF(ID.EQ.2) THEN
C Nested annuli, radial limits
                IF(NL.EQ.0.OR.MOD(NL,2).NE.1) THEN
                    WRITE(*,*) 'QRT_APERTURE error -  bad radial limits'
                    ISTAT=1
                    RETURN
                ENDIF
                IF(NL.GT.MAXRS) THEN
                    WRITE(*,*) 'QRT_APERTURE error -  too many annuli'
                    ISTAT=1
                    RETURN
                ENDIF
                IT=4
                CALL SRT_NSUR(KSUR,ISTAT)
                DO J=1,NL
                        PL(9+J)=ALIM(J)
                ENDDO
                PL(NL+10)=KSUR
                PL(NL+11)=NSURF
                NP=11+NL
                IDEF(2)=(NL-1)/2
        ELSEIF(ID.EQ.3.OR.ID.EQ.4) THEN
C Cartesian limits
                IF(NL.NE.4.OR.ALIM(1).GT.ALIM(3).OR.
     +          ALIM(2).GT.ALIM(4)) THEN
                 WRITE(*,*) 'QRT_APERTURE error-  bad cartesian limits'
                 ISTAT=1
                 RETURN
                ENDIF
                IF(ID.EQ.3) THEN
                        IT=1
                ELSE
                        IT=5
                ENDIF
                NP=9+NL
                DO J=1,NL
                        PL(9+J)=ALIM(J)
                ENDDO
                IDEF(2)=0
        ELSEIF(ID.EQ.5) THEN
C Cartesian grid limits pitchx pitchy ribx riby
                IF(NL.NE.4.OR.ALIM(1).LT.ALIM(3).OR.
     +          ALIM(2).LT.ALIM(4)) THEN
                   WRITE(*,*) 'QRT_APERTURE - bad cartesian grid limits'
                   ISTAT=1
                   RETURN
                ENDIF
C surface type 18
                IT=18
                NP=9+NL
                DO J=1,NL
                        PL(9+J)=ALIM(J)
                ENDDO
                IDEF(2)=0
        ELSEIF(ID.EQ.6) THEN
C Radial/azimuthal sector limits rmin,rmax,amin,amax
C amin and amax in radians
                IF(NL.NE.4.OR.ALIM(1).GT.ALIM(2).OR.
     +          ALIM(3).GT.ALIM(4)) THEN
                    WRITE(*,*) 'QRT_APERTURE error - bad sector limits'
                    ISTAT=1
                    RETURN
                ENDIF
C surface type 22
                IT=22
                NP=9+NL
                DO J=1,NL
                        PL(9+J)=ALIM(J)
                ENDDO
                IDEF(2)=0
        ELSEIF(ID.EQ.7) THEN
C Cartesian limits
                IF(NL.NE.5.OR.ALIM(1).GT.ALIM(3).OR.
     +          ALIM(2).GT.ALIM(4)) THEN
                 WRITE(*,*) 'QRT_APERTURE error -  bad parallel limits'
                 ISTAT=1
                 RETURN
                ENDIF
                IT=27
                NP=9+NL
                DO J=1,NL
                        PL(9+J)=ALIM(J)
                ENDDO
                IDEF(2)=0
        ELSEIF(ID.EQ.8) THEN
C MCO test station aperture
                IF(NL.NE.2.OR.ALIM(1).GT.ALIM(2)) THEN
                 WRITE(*,*) 'QRT_APERTURE error -  ID=8 bad paramters'
                 ISTAT=1
                 RETURN
                ENDIF
                IT=29
                NP=9+NL
                DO J=1,NL
                        PL(9+J)=ALIM(J)
                ENDDO
                IDEF(2)=0
        ENDIF
C Finally set parameters in common for stop
        CALL SRT_SETF(0,IT,NP,PL,IDEF,0,0,-1,ISTAT)
        END        
