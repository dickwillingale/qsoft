*+QRI_SETSKY        Set sky coordinates for image field
        SUBROUTINE QRI_SETSKY(XTODEG,YTODEG,IPR,MJD,RA,DEC,ROLL)
        IMPLICIT NONE
        INTEGER IPR
        DOUBLE PRECISION XTODEG,YTODEG,MJD,RA,DEC,ROLL
Cf2py   intent(in) XTODEG,YTODEG,IPR,MJD,RA,DEC,ROLL
*XTODEG     input    scale from local X to degrees (usually -ve)
*YTODEG     input    scale from local Y to degrees
*IPR        input    projection between local XY and local spherical
*                0 Plate Carre
*                1 Aitoff (Hammer)
*MJD        input    Modified Julian date 
*R          input    Right Ascension (degrees J2000 at local origin)
*DEC        input    Declination (degrees J2000 at local origin)
*ROLL       input    Roll angle (degrees from North to +ve elev. +ve clockwise)
*-Author Dick Willingale 2012-May-11
        INCLUDE 'QRI_TRANSCOM'
        INCLUDE 'QR_COM'

        DOUBLE PRECISION DAXIS(2),OBLIQECL,DCEL(2),DTOR,T
C
        IF(ISTAT.NE.0) RETURN
C
        DTOR=ASIN(1.0D0)/90.D0
        XYTORAD(1)=XTODEG*DTOR
        XYTORAD(2)=YTODEG*DTOR
C Set Ra, Dec and Roll
        DAXIS(1)=RA*DTOR
        DAXIS(2)=DEC*DTOR
        QRA=DAXIS(1)
        QDEC=DAXIS(2)
        QMJD=MJD
        QROLL=ROLL*DTOR
C Obliquity of Ecliptic from SLA_ECMAT
        T=(QMJD-51544.5D0)/36525D0
        OBLIQECL=
     +  (84381.448D0+(-46.8150D0+(-0.00059D0+0.001813D0*T)*T)*T)
        OBLIQECL=DTOR*OBLIQECL/3600.D0
        DCEL(1)=0.D0
        DCEL(2)=DCEL(1)
        CALL AX_DMAT(DCEL,OBLIQECL,CTOE,ETOC)
        CALL AX_DMAT(DAXIS,QROLL,CTOS,STOC)
        CALL AX_DONMXM(STOC,CTOG,STOG)
        CALL AX_DONMXM(ETOC,CTOS,ETOS)
        CALL AX_DONMXM(STOC,CTOE,STOE)
C
        IPROJ=IPR
        CALL QRI_STOL(DCEL,DAXIS)
        AENOW(1)=DCEL(1)
        AENOW(2)=DCEL(2)
        XYNOW(1)=DAXIS(1)
        XYNOW(2)=DAXIS(2)
C
        SKYTRANS=.TRUE.
        END
