*+SRT_SETSLECOM    Set values in constellation common block for SLE
        SUBROUTINE SRT_SETSLECOM(IPACK,FL,PLMIN,PLMAX,SLOT,
     +        RMIN,RMAX,CWIDTH,CHEIGHT,NMAX,RC,PC,TC,AC,NSET)
*IPACK        input        module packing 0 1 module, 1 sunflower, 2 cartesian
*FL        input        focal length
*PLMIN        input        minimum module length
*PLMAX        input        maximum module length
*SLOT        input        slot aperture width
*RMIN        input        minimum radius of aperture
*RMAX        input        maximum radius of aperture
*CSIZE        input        size of each module mm
*NMAX        input        maximum number of element positions returned
*RC,PC        output        radius and azimuth of each module
*TC        output        rotation angle of module
*AC        output        axial length of each module
*NSET        output        number of elements set
* IF IPACK=0 then single module with centre at x=y=(rmax+rmin)/2/SQRT(2)
*
              IMPLICIT NONE
        INTEGER IPACK,NMAX,NSET
        DOUBLE PRECISION FL,PLMIN,PLMAX,SLOT,CWIDTH,CHEIGHT
        DOUBLE PRECISION RMIN,RMAX,RC(NMAX),PC(NMAX),TC(NMAX)
        DOUBLE PRECISION AC(NMAX)
*-Derived from SRT_SETKBSCOM by Vladimir Tichy
* Last modification 24 Jan 2018
        INCLUDE 'SRT_COM'
        INTEGER J,NDO
        DOUBLE PRECISION THETAG
C
C#define VDEBUG
C
        inrad=RMIN
        outrad=RMAX
        NDO=MIN(maxlist,NMAX)
C there should be CSIZE instead of 0 when packing is used
        CALL srt_slemakeconstellation(IPACK,inrad,outrad,0,
     +          NDO,rlist,philist,thlist,nrofelts)
        iptype=1
        DO J=1,nrofelts
                wlist(J)=CWIDTH
                hlist(J)=CHEIGHT
                THETAG=ATAN2(rlist(J)/SQRT(2.0),ABS(FL))/2.0
                llist(J)=MIN(MAX(PLMIN,SLOT/THETAG),PLMAX)
                IF(FL.LT.0.0) THEN
C If 2nd stack then decrease radius to correct position for top of module
                        rlist(j)=rlist(j)*(FL*2.0)/(FL*2.0+PLMAX)
                ENDIF
        ENDDO
C 
        NSET=MIN(NMAX,nrofelts)
        DO J=1,NSET
                IF(J.LE.NMAX) THEN
                        RC(J)=rlist(J)
                        PC(J)=philist(J)
                        TC(J)=thlist(J)
                        AC(J)=llist(J)
                ENDIF
        ENDDO
c#ifdef VDEBUG
c        WRITE(*,*) 'DEBUG: MODULES SET:', nrofelts
c        DO J=1,nrofelts
c                WRITE(*,*) '>',j,': r=',rlist(j),' phi=',philist(j),
c     +          ' theta=',thlist(j),' al=',llist(j)
c                WRITE(*,*) 'wlist :',wlist(J)
c                WRITE(*,*) 'hlist :',hlist(J)
c        ENDDO
c#endif
        END
