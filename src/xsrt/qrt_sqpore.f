*+QRT_SQPORE Set up slumped square pore MCP
        SUBROUTINE QRT_SQPORE(PCEN,PNORM,RAXIS,RCUR,IPACK,RAP,PITCH,
     +  WALL,PLEN,IDF,IQ,PLMIN,PLMAX,FIBRE,NAR,AR)
*PCEN       input        position of centre of plate
*RNORM      input        normal at centre of plate
*RAXIS      input        reference axi at centre of plate
*RCUR       input        radius of curvature
*IPACK      input        pore packing
*               1 cart, 2 rad, 3 waff, 4 octag, 5 rand, 6 MIXS, 7 NFL
*RAP        input        half width of plate aperture
*PITCH      input        pitch of pores on a cartesian grid
*WALL       input        pore wall thickness
*PLEN       input        length of pores
*IDF        input        deformation index
*IQ         input        surface quality index
*PLMIN      input        minimum pore length
*PLMAX      input        maximum pore length
*FIBRE      input        size of fibre bundle in packing
*NAR        input        number of additional parameters
*AP         input        array of additional parameters
Cf2py  intent(in) PCEN,PNORM,RAXIS,RCUR,IPACK,RAP,PITCH
Cf2py  intent(in) WALL,PLEN,IDF,IQ,PLMIN,PLMAX,FIBRE,NAR,AR
        IMPLICIT NONE
        INTEGER IPACK,IDF,IQ,NAR
        DOUBLE PRECISION PCEN(3),PNORM(3),RAXIS(3),RCUR,RAP
        DOUBLE PRECISION PITCH,WALL,PLEN,PLMIN,PLMAX,AR(NAR)
        DOUBLE PRECISION FIBRE
*-Author Dick Willingale 2012-May-20
        INTEGER IDEF(2),J,NTHIS,NP,MAXP,NMOD,ITMOD
        PARAMETER (MAXP=50)
        DOUBLE PRECISION OTHER(3),SCEN(3),PP(MAXP)
        INCLUDE 'QR_COM'
C
        IF(ISTAT.NE.0) RETURN
C Set surface deformation
        IDEF(1)=IDF
        IDEF(2)=2
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
        PP(11)=RAP
        PP(12)=PITCH
        PP(13)=PITCH
        PP(14)=WALL
        PP(15)=WALL
        PP(16)=PLEN
C The next 2 are place markers for the x,y of the ray intersection in locals
        PP(17)=0.0
        PP(18)=0.0
C index of surface and surface quality
        CALL SRT_NSUR(NTHIS,ISTAT)
        PP(19)=NTHIS
        PP(20)=IQ
C Packing of pores
        PP(21)=IPACK
C Pore length range
        PP(22)=PLMIN
        PP(23)=PLMAX
C Size of FIBRE in packing
        PP(24)=FIBRE
        NP=24
        IF(IPACK.EQ.6.OR.IPACK.EQ.7) THEN
C Additional aperture parameters for MIXS and NFL implementation
              IF(MOD(NAR,6).NE.0) THEN
               WRITE(*,*) 'QRT_SQPORE - incorrect MIXS or NFL ap. defn'
               ISTAT=1
               RETURN
              ENDIF
              NMOD=NAR/6
              if(IPACK.EQ.6) THEN
                        ITMOD=2
              else
                        ITMOD=1
              endif
              CALL SRT_SETCONCOM(ITMOD,0.D0,RAP,NMOD,
     +        AR(1),AR(NMOD+1),AR(NMOD*2+1),AR(NMOD*3+1),
     +        AR(NMOD*4+1),AR(NMOD*5+1),AR(NMOD*5+1),
     +        AR(NMOD*5+1),ISTAT)
        ENDIF
C Set parameters in common for sphere with radial/cart+return limits, type=21
        CALL SRT_SETF(0,21,NP,PP,IDEF,0,0,-1,ISTAT)
C Set place marker for parallelogram aperture of pore
        CALL SRT_SETF(0,1,14,PP,0,0,0,-1,ISTAT)
C Set place markers for 4 sides of pore
        CALL SRT_SETF(0,5,13,PP,0,IQ,NTHIS+3,-1,ISTAT)
        CALL SRT_SETF(0,5,13,PP,0,IQ,NTHIS+4,-1,ISTAT)
        CALL SRT_SETF(0,5,13,PP,0,IQ,NTHIS+5,-1,ISTAT)
        CALL SRT_SETF(0,5,13,PP,0,IQ,NTHIS+2,-1,ISTAT)
        END
