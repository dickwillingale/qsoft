*+QRT_DETECTOR Set detector parameters
        SUBROUTINE QRT_DETECTOR(ID,DPOS,DNML,DRFX,DLIM,RADET)
*ID        input        detector type
*DPOS      input        detector position
*DNML      input        detector normal
*DRFX      input        detector reference axis
*DLIM      input        detector limits
*RADET     input        radius of curvature of spherical detector
Cf2py  intent(in) ID,DPOS,DNML,DRFX,DLIM,RADET
        IMPLICIT NONE
        INTEGER ID
        DOUBLE PRECISION DPOS(3),DNML(3),DRFX(3),DLIM(4),RADET
*-Author Dick Willingale 2012-Apr-30
        INTEGER IT,IDEF(2),NP,IQ
        DOUBLE PRECISION PL(14)
        INCLUDE 'QR_COM'
C
        IF(ISTAT.NE.0) RETURN
C Get detector surface parameters
        PL(1)=DNML(1)
        PL(2)=DNML(2)
        PL(3)=DNML(3)
              PL(4)=DRFX(1)
              PL(5)=DRFX(2)
              PL(6)=DRFX(3)
        PL(7)=DPOS(1)
        PL(8)=DPOS(2)
        PL(9)=DPOS(3)
        IF(ID.EQ.1) THEN
C Planar detector, radial limits
                IF(DLIM(1).GT.DLIM(2)) THEN
                      WRITE(*,*) 'TRS_DETR error -  bad radial limits'
                      ISTAT=1
                      RETURN
                ENDIF
                PL(10)=DLIM(1)
                PL(11)=DLIM(2)
                IT=6
                NP=11
        ELSEIF(ID.EQ.2) THEN
C Planar detector cartesian limits
                IF(DLIM(1).GT.DLIM(3).OR.DLIM(2).GT.DLIM(4)) THEN
                    WRITE(*,*) 'TRS_DETR error -  bad cartesian limits'
                    ISTAT=1
                    RETURN
                ENDIF
                PL(10)=DLIM(1)
                PL(11)=DLIM(2)
                PL(12)=DLIM(3)
                PL(13)=DLIM(4)
                IT=5
                NP=13
        ELSEIF(ID.EQ.3) THEN
C Spherical detector, radial limits
                IF(DLIM(1).GT.DLIM(2)) THEN
                    WRITE(*,*) 'TRS_DETR error -  bad radial limits'
                    ISTAT=1
                    RETURN
                ENDIF
C Get radius and set centre of sphere
                PL(7)=DPOS(1)-DNML(1)*RADET
                PL(8)=DPOS(2)-DNML(2)*RADET
                PL(9)=DPOS(3)-DNML(3)*RADET
C Move limits and put radius into parameter array
                PL(10)=RADET
                PL(11)=DLIM(1)
                PL(12)=DLIM(2)
                IT=15
                NP=12
        ELSEIF(ID.EQ.4) THEN
C Spherical detector cartesian limits
                IF(DLIM(1).GT.DLIM(3).OR.DLIM(2).GT.DLIM(4)) THEN
                    WRITE(*,*) 'TRS_DETR error -  bad cartesian limits'
                    ISTAT=1
                    RETURN
                ENDIF
C Get radius and set centre of sphere
                PL(7)=DPOS(1)-DNML(1)*RADET
                PL(8)=DPOS(2)-DNML(2)*RADET
                PL(9)=DPOS(3)-DNML(3)*RADET
C Move limits and put radius into parameter array
                PL(10)=RADET
                PL(11)=DLIM(1)
                PL(12)=DLIM(2)
                PL(13)=DLIM(3)
                PL(14)=DLIM(4)
                IT=14
                NP=14
        ENDIF
C Finally set parameters in common for detector
C Note that the surface quality is -ve for a detector
        IDEF(1)=0
        IDEF(2)=0
        IQ=-1
        CALL SRT_SETF(0,IT,NP,PL,IDEF,IQ,0,-1,ISTAT)
        END
