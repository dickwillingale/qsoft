*+QRT_BAFFLEARRAY Set up an array of baffle elements
        SUBROUTINE QRT_BAFFLEARRAY(NA,AP,AN,AR,RMIN,RMAX,AMIN,AMAX)
*NA       input        number of baffle elements
*AP       input        positions of baffle elements
*AN       input        normals to baffle elements
*AR       input        reference axes in baffle element planes
*RMIN     input        baffle element minimum radius
*RMAX     input        baffle element maximum radius
*AMIN     input        baffle element minimum azimuth radians
*AMAX     input        baffle element maximum azimuth radians
Cf2py  intent(in) NA,AP,AN,AR,RMIN,RMAX,AMIN,AMAX
        IMPLICIT NONE
        INTEGER NA
        DOUBLE PRECISION AP(3,NA),AN(3,NA),AR(3,NA),RMIN(NA),RMAX(NA)
        DOUBLE PRECISION AMIN(NA),AMAX(NA)
*-Author Dick Willingale 2019-Jun-03
        INTEGER IDEF(2),J,NEXT
        DOUBLE PRECISION PL(14)
        INCLUDE 'QR_COM'
C
        IF(ISTAT.NE.0) RETURN
C Set deformation
        IDEF(1)=0
        IDEF(2)=0
C Find surface index after array
        CALL SRT_NSUR(NEXT,ISTAT)
        NEXT=NEXT+NA
C Loop for all apertures
        DO J=1,NA
C Set surface parameters
                PL(1)=AN(1,J)
                PL(2)=AN(2,J)
                PL(3)=AN(3,J)
                PL(4)=AR(1,J)
                PL(5)=AR(2,J)
                PL(6)=AR(3,J)
                PL(7)=AP(1,J)
                PL(8)=AP(2,J)
                PL(9)=AP(3,J)
C Cartesian limits on surface
                PL(10)=RMIN(J)
                PL(11)=RMAX(J)
                PL(12)=AMIN(J)
                PL(13)=AMAX(J)
                PL(14)=NA
C Finally set parameters in common
                CALL SRT_SETF(0,33,14,PL,IDEF,0,NEXT,NEXT,ISTAT)
        ENDDO
        END
