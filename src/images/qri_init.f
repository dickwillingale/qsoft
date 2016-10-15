*+QRI_INIT Initialisation for image processing
        SUBROUTINE QRI_INIT()
*-Author Dick Willingale 2012-May-11
        INCLUDE 'QRI_TRANSCOM'
C The following from SLA_EQGAL for 2000 epoch
        DATA CTOG/-0.054875539726D0,-0.873437108010D0,-0.483834985808D0,
     :            +0.494109453312D0,-0.444829589425D0,+0.746982251810D0,
     :            -0.867666135858D0,-0.198076386122D0,+0.455983795705D0/
        DATA GTOC/-0.054875539726D0,+0.494109453312D0,-0.867666135858D0,
     :            -0.873437108010D0,-0.444829589425D0,-0.198076386122D0,
     :            -0.483834985808D0,+0.746982251810D0,+0.455983795705D0/
        DOUBLE PRECISION DAXIS(2),OBLIQECL,DCEL(2),DTOR,T
C
        CALL QR_INIT()
C
        DTOR=ASIN(1.0D0)/90.D0
C Default is to hide sky transformations although they are initialised
        SKYTRANS=.FALSE.
C Set Ra, Dec and Roll
        DAXIS(1)=0.D0
        DAXIS(2)=0.D0
        QRA=DAXIS(1)
        QDEC=DAXIS(2)
        QROLL=0.D0
        QMJD=45000.0D0
C Obliquity of Ecliptic from SLA_ECMAT
       T=(QMJD-51544.5D0)/36525D0
       OBLIQECL=(84381.448D0+(-46.8150D0+(-0.00059D0+0.001813D0*T)*T)*T)
       OBLIQECL=DTOR*OBLIQECL/3600.D0
       DCEL(1)=0.D0
       DCEL(2)=0.D0
        CALL AX_DMAT(DCEL,OBLIQECL,CTOE,ETOC)
        CALL AX_DMAT(DAXIS,QROLL,CTOS,STOC)
        CALL AX_DONMXM(GTOC,CTOS,GTOS)
        CALL AX_DONMXM(STOC,CTOG,STOG)
        CALL AX_DONMXM(ETOC,CTOS,ETOS)
        CALL AX_DONMXM(STOC,CTOE,STOE)
C Set default pixel sizes etc.
        XYTORAD(1)=-DTOR
        XYTORAD(2)=DTOR
        BLHC(1)=0.0
        BLHC(2)=0.0
        XYPIXEL(1)=DTOR
        XYPIXEL(2)=DTOR
        XYNOW(1)=DCEL(1)
        XYNOW(2)=DCEL(2)
        IPROJ=0
        CALL QRI_LTOS(DCEL,DAXIS)
        AENOW(1)=DAXIS(1)
        AENOW(2)=DAXIS(2)
        END
*+QRI_PTOL        Pixel to local XY coordinate transformation
        SUBROUTINE QRI_PTOL(PIX,XYP)
*PIX        input        PIX position 0 to NX-1, 0 to NY-1
*XYP        outout        local XY position
        IMPLICIT NONE
        DOUBLE PRECISION PIX(2),XYP(2)
*-Author Dick Willingale 2012-May-11
        INCLUDE 'QRI_TRANSCOM'
        INCLUDE 'QR_COM'
C
        IF(ISTAT.NE.0) RETURN
C
        XYP(1)=BLHC(1)+PIX(1)*XYPIXEL(1)
        XYP(2)=BLHC(2)+PIX(2)*XYPIXEL(2)
        END
*+QRI_LTOP        local XY to pixel coordinate transformation
        SUBROUTINE QRI_LTOP(XYP,PIX)
*XYP        input        local XY position
*PIX        output        PIX position 0 to NX-1, 0 to NY-1
        IMPLICIT NONE
        DOUBLE PRECISION PIX(2),XYP(2)
*-Author Dick Willingale 2012-May-11
        INCLUDE 'QRI_TRANSCOM'
        INCLUDE 'QR_COM'
C
        IF(ISTAT.NE.0) RETURN
C
        PIX(1)=(XYP(1)-BLHC(1))/XYPIXEL(1)
        PIX(2)=(XYP(2)-BLHC(2))/XYPIXEL(2)
        END
*+QRI_LTOS        local XY to local radians coordinate transformation
        SUBROUTINE QRI_LTOS(XYP,LOC)
*XYP        input        local XY position
*LOC        output        local position Azimuth,Elevation radians
        IMPLICIT NONE
        DOUBLE PRECISION XYP(2),LOC(2)
*-Author Dick Willingale 2012-May-11
        INCLUDE 'QRI_TRANSCOM'
        INCLUDE 'QR_COM'
        DOUBLE PRECISION PIBY2
        INTEGER STATUS
C
        IF(ISTAT.NE.0) RETURN
C
        PIBY2=ASIN(1.0D0)
        STATUS=0
C
        IF(IPROJ.EQ.0) THEN
                LOC(1)=XYP(1)*XYTORAD(1)
                LOC(2)=XYP(2)*XYTORAD(2)
                IF(ABS(LOC(1)).GT.PIBY2*2.D0.OR.
     +            ABS(LOC(2)).GT.PIBY2) THEN
                    STATUS=1
                ENDIF
        ELSEIF(IPROJ.EQ.1) THEN
                CALL AX_DAMMERINV(XYP(1),XYP(2),LOC(1),LOC(2),STATUS)
        ELSEIF(IPROJ.EQ.2) THEN
                CALL AX_LAMBERTINV(XYP(1),XYP(2),LOC(1),LOC(2),STATUS)
        ENDIF
        IF(STATUS.NE.0) THEN
                LOC(1)=0.0
                LOC(2)=0.0
        ENDIF
        END
*+QRI_STOL        local spherical radians to local XYs coordinate transformation
        SUBROUTINE QRI_STOL(LOC,XYP)
*LOC        input        local position Azimuth,Elevation radians
*XYP        output        local XY position
        IMPLICIT NONE
        DOUBLE PRECISION XYP(2),LOC(2)
*-Author Dick Willingale 2012-May-11
        INCLUDE 'QRI_TRANSCOM'
        INCLUDE 'QR_COM'
C
        IF(ISTAT.NE.0) RETURN
C
        IF(IPROJ.EQ.0) THEN
                XYP(1)=LOC(1)/XYTORAD(1)
                XYP(2)=LOC(2)/XYTORAD(2)
        ELSEIF(IPROJ.EQ.1) THEN
                CALL AX_DAMMER(LOC(1),LOC(2),XYP(1),XYP(2))
        ELSEIF(IPROJ.EQ.2) THEN
                CALL AX_LAMBERT(LOC(1),LOC(2),XYP(1),XYP(2))
        ENDIF
        END
*+QRI_STOC        Local to Celestial coordinate transformation
        SUBROUTINE QRI_STOC(LOC,CEL)
*LOC        input        Local position Azimuth,Elevation radians
*CEL        output        Celestial position RA,DEC radians
        IMPLICIT NONE
        DOUBLE PRECISION CEL(2),LOC(2)
*-Author Dick Willingale 2012-May-11
        INCLUDE 'QRI_TRANSCOM'
        INCLUDE 'QR_COM'
C
        IF(ISTAT.NE.0) RETURN
C
        CALL AX_DONVRT(LOC(1),LOC(2),STOC,CEL(1),CEL(2))
        END
*+QRI_CTOS        Celestial to local coordinate transformation
        SUBROUTINE QRI_CTOS(CEL,LOC)
*CEL        input        Celestial position RA,DEC radians
*LOC        output        Local position Azimuth,Elevation radians
        IMPLICIT NONE
        DOUBLE PRECISION CEL(2),LOC(2)
*-Author Dick Willingale 2012-May-11
        INCLUDE 'QRI_TRANSCOM'
        INCLUDE 'QR_COM'
C
        IF(ISTAT.NE.0) RETURN
C
        CALL AX_DONVRT(CEL(1),CEL(2),CTOS,LOC(1),LOC(2))
        END
*+QRI_CTOE        Celestial to Ecliptic coordinate transformation
        SUBROUTINE QRI_CTOE(CEL,ECL)
*CEL        input        Celestial position RA,DEC radians
*ECL        output        Ecliptic position ECL,ECD radians
        IMPLICIT NONE
        DOUBLE PRECISION CEL(2),ECL(2)
*-Author Dick Willingale 2012-May-11
        INCLUDE 'QRI_TRANSCOM'
        INCLUDE 'QR_COM'
C
        IF(ISTAT.NE.0) RETURN
C
        CALL AX_DONVRT(CEL(1),CEL(2),CTOE,ECL(1),ECL(2))
        END
*+QRI_ETOC        Ecliptic to Celestial coordinate transformation
        SUBROUTINE QRI_ETOC(ECL,CEL)
*ECL        output        Ecliptic position ECL,ECD radians
*CEL        input        Celestial position RA,DEC radians
        IMPLICIT NONE
        DOUBLE PRECISION CEL(2),ECL(2)
*-Author Dick Willingale 2012-May-11
        INCLUDE 'QRI_TRANSCOM'
        INCLUDE 'QR_COM'
C
        IF(ISTAT.NE.0) RETURN
C
        CALL AX_DONVRT(ECL(1),ECL(2),ETOC,CEL(1),CEL(2))
        END
*+QRI_CTOG        Celestial to Galactic coordinate transformation
        SUBROUTINE QRI_CTOG(CEL,GAL)
*CEL        input        Celestial position RA,DEC radians
*GAL        output        Galactic position LII,BII radians
        IMPLICIT NONE
        DOUBLE PRECISION CEL(2),GAL(2)
*-Author Dick Willingale 2012-May-11
        INCLUDE 'QRI_TRANSCOM'
        INCLUDE 'QR_COM'
C
        IF(ISTAT.NE.0) RETURN
C
        CALL AX_DONVRT(CEL(1),CEL(2),CTOG,GAL(1),GAL(2))
        END
*+QRI_GTOC        Galactic to Celestial coordinate transformation
        SUBROUTINE QRI_GTOC(GAL,CEL)
*GAL        output        Galactic position LII,BII radians
*CEL        input        Celestial position RA,DEC radians
        IMPLICIT NONE
        DOUBLE PRECISION CEL(2),GAL(2)
*-Author Dick Willingale 2012-May-11
        INCLUDE 'QRI_TRANSCOM'
        INCLUDE 'QR_COM'
C
        IF(ISTAT.NE.0) RETURN
C
        CALL AX_DONVRT(GAL(1),GAL(2),GTOC,CEL(1),CEL(2))
        END
