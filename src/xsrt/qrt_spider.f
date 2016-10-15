*+QRT_SPIDER Set support spider parameters
        SUBROUTINE QRT_SPIDER(CONE,APOS,ANML,ARFX,NSEC,CWID,AWID)
*CONE        input        90-half cone angle degrees (0.0 for plane)
*APOS        input        axial position of vertex (centre)
*ANML        input        direction of normal to aperture (optic axis)
*ARFX        input        reference axis in aperture
*NSEC        input        number of sectors (number of arms)
*CWID        input        constant arm width
*AWID        input        angular arm width degrees
Cf2py  intent(in) CONE,APOS,ANML,ARFX,NSEC,CWID,AWID
        IMPLICIT NONE
        INTEGER NSEC
        DOUBLE PRECISION CONE,APOS(3),ANML(3),ARFX(3),CWID,AWID
*-Author Dick Willingale 2012_May-01
        INCLUDE 'SRT_COM'
        INTEGER IDEF(2)
        DOUBLE PRECISION PL(15)
        INCLUDE 'QR_COM'
C
        IF(ISTAT.NE.0) RETURN
C
        PL(1)=ANML(1)
        PL(2)=ANML(2)
        PL(3)=ANML(3)
        PL(4)=ARFX(1)
        PL(5)=ARFX(2)
        PL(6)=ARFX(3)
        PL(7)=APOS(1)
        PL(8)=APOS(2)
        PL(9)=APOS(3)
        IDEF(1)=0
        IDEF(2)=0
        IF(CONE.EQ.0.0) THEN
C Plane spider
                PL(10)=NSEC
                PL(11)=CWID
                PL(12)=AWID*PI/180.D0
                CALL SRT_SETF(0,17,12,PL,IDEF,0,0,-1,ISTAT)
        ELSE
C Conical spider
                PL(10)=TAN((90.D0-CONE)*PI/180.D0)**2
                PL(11)=0.0
                PL(12)=0.0
                PL(13)=NSEC
                PL(14)=CWID
                PL(15)=AWID*PI/180.D0
                CALL SRT_SETF(0,16,15,PL,IDEF,0,0,-1,ISTAT)
        ENDIF
        END        
