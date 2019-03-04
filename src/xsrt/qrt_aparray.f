*+QRT_APARRAY Set up an array of rectangular apertures
        SUBROUTINE QRT_APARRAY(NA,AP,AN,AR,XHW,YHW)
*NA       input        number of apertures
*AP       input        positions of apertures
*AN       input        normals to aperture planes
*AR       input        reference axes in aperture planes
*XHW      input        aperture half widths in reference axes
*YHW      input        aperture half widths in remaining axes
Cf2py  intent(in) NA,AP,AN,AR,XHW,YHW
        IMPLICIT NONE
        INTEGER NA
        DOUBLE PRECISION AP(3,NA),AN(3,NA),AR(3,NA),XHW(NA),YHW(NA)
*-Author Dick Willingale 2019-Mar-02
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
                PL(10)=-XHW(J)
                PL(11)=-YHW(J)
                PL(12)=XHW(J)
                PL(13)=YHW(J)
                PL(14)=NA
C Finally set parameters in common
                CALL SRT_SETF(0,32,14,PL,IDEF,0,NEXT,NEXT,ISTAT)
        ENDDO
        END
