*+QRT_SQMPOARR Set up square pore MPO array
        SUBROUTINE QRT_SQMPOARR(PCEN,PNORM,RAXIS,RCUR,HWID,IDF,NAR,AR)
*PCEN       input        position of centre of array
*RNORM      input        normal at centre of array
*RAXIS      input        reference axi at centre of array
*RCUR       input        radius of curvature of array
*HWID       input        half width of array aperture
*IDF        input        deformation index for array
*NAR        input        number of additional parameters
*AP         input        array of additional parameters for each MPO in array
Cf2py  intent(in) PCEN,PNORM,RAXIS,RCUR,HWID,IDF,NAR,AR
        IMPLICIT NONE
        INTEGER IDF,NAR
        DOUBLE PRECISION PCEN(3),PNORM(3),RAXIS(3),RCUR,HWID
        DOUBLE PRECISION AR(NAR)
*-Author Dick Willingale 2017-Oct-23
        INTEGER IDEF(2),J,NTHIS,NP,MAXP,NMOD,IQ
        PARAMETER (MAXP=50)
        DOUBLE PRECISION OTHER(3),SCEN(3),PP(MAXP)
        INCLUDE 'QR_COM'
C
        IF(ISTAT.NE.0) RETURN
C Set surface deformation
        IDEF(1)=IDF
        IDEF(2)=2
C Set default surface quality 
        IQ=0
C normalise reference vectors
        CALL SRT_VNRM(PNORM,ISTAT)
        CALL SRT_VNRM(RAXIS,ISTAT)
C Find centre of sphere
        SCEN(1)=PCEN(1)-PNORM(1)*RCUR
        SCEN(2)=PCEN(2)-PNORM(2)*RCUR
        SCEN(3)=PCEN(3)-PNORM(3)*RCUR
C force reference axis to be at 90 degrees to normal
        CALL SRT_VCRS(PNORM,RAXIS,OTHER)
        CALL SRT_VNRM(OTHER,ISTAT)
        CALL SRT_VCRS(OTHER,PNORM,RAXIS)
        CALL SRT_VNRM(RAXIS,ISTAT)
C Pack parameters into single array
        DO J=1,3
                PP(J)=PNORM(J)
                PP(J+3)=RAXIS(J)
                PP(J+6)=SCEN(J)
        ENDDO
        PP(10)=RCUR
        PP(11)=HWID
C The following are void (spare) 
        PP(12)=0.0
        PP(13)=0.0
        PP(14)=0.0
        PP(15)=0.0
        PP(16)=0.0
C The next 2 are place markers for the x,y of the ray intersection in locals
        PP(17)=0.0
        PP(18)=0.0
C index of surface and surface quality
        CALL SRT_NSUR(NTHIS,ISTAT)
        PP(19)=NTHIS
        NP=19
C Additional MPO parameters
        IF(MOD(NAR,15).NE.0) THEN
               WRITE(*,*)
     +         'QRT_SQMPOARR - error - requires 15 parameters per MPO'
               ISTAT=1
               RETURN
        ENDIF
        NMOD=NAR/15
        CALL SRT_SETMPOARR(HWID,NMOD,
     +        AR(1),AR(NMOD+1),AR(NMOD*2+1),AR(NMOD*3+1),
     +        AR(NMOD*4+1),AR(NMOD*5+1),AR(NMOD*6+1),
     +        AR(NMOD*7+1),AR(NMOD*8+1),AR(NMOD*9+1),AR(NMOD*10+1),
     +        AR(NMOD*11+1),AR(NMOD*12+1),AR(NMOD*13+1),
     +        AR(NMOD*14+1),ISTAT)
C Set parameters in common for sphere with radial/cart+return limits, type=28
        CALL SRT_SETF(0,28,NP,PP,IDEF,0,0,-1,ISTAT)
C Set place marker for parallelogram aperture of pore
        CALL SRT_SETF(0,1,14,PP,0,0,0,-1,ISTAT)
C Set place markers for 4 sides of pore
        CALL SRT_SETF(0,5,13,PP,0,IQ,NTHIS+3,-1,ISTAT)
        CALL SRT_SETF(0,5,13,PP,0,IQ,NTHIS+4,-1,ISTAT)
        CALL SRT_SETF(0,5,13,PP,0,IQ,NTHIS+5,-1,ISTAT)
        CALL SRT_SETF(0,5,13,PP,0,IQ,NTHIS+2,-1,ISTAT)
        END
