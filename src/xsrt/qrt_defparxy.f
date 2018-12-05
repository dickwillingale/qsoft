*+QRT_DEFPARXY Set deformation parameters in x and y axes
        SUBROUTINE QRT_DEFPARXY(ID,IS,KNX,KNY,XPAR,YPAR)
*ID        input        deformation index
*IM        input        sub-surface index
*KNX       input        number of x parameters
*KNY       input        number of y parameters
*XSAM      input        x values
*YSAM      input        y values
Cf2py  intent(in) ID,IM,KNX,KNY,XPAR,YPAR
        IMPLICIT NONE
        INTEGER ID,IS
        INTEGER KNX,KNY
        DOUBLE PRECISION XPAR(KNX),YPAR(KNY)
*-Author Vladimir Tichy (2017)
        INCLUDE 'SRT_COM'
        INTEGER IX,IY
        INTEGER NS,NX,NY
        INCLUDE 'QR_COM'
C
c	WRITE (*,*) 'QRT_DEFPARXY DEBUG   KNX = ',KNX
c	WRITE (*,*) 'QRT_DEFPARXY DEBUG   KNY = ',KNY
        IF(ISTAT.NE.0) RETURN
C Check index in range
        IF(ID.LT.1.OR.ID.GT.MAXDF) THEN
         WRITE(*,*) 'QRT_DEFPARXY err, dfrm. index out of range',ID
         ISTAT=1
         RETURN
        ENDIF
C Check sub-matrix index in range
        IF(IS.LT.1.OR.IS.GT.IDFM(1,ID)) THEN
         WRITE(*,*) 'QRT_DEFPARXY err, subsurf. index out of range',ID
         ISTAT=1
         RETURN
        ENDIF
C Get indices from common
        NS=IDFM(1,ID)
        NX=IDFM(2,ID)
        NY=IDFM(3,ID)
        IX=IDFP(1,ID)
        IY=IDFP(2,ID)
C Check array sizes
        IF(KNX.GT.NX.OR.KNY.GT.NY) THEN
                WRITE(*,*) 'QRT_DEFPARXY error - dimensions incorrect'
                ISTAT=1
        ENDIF
C Copy parameters
        CALL SRT_DEFPARXY(NS,IS,KNX,KNY,XPAR,YPAR,
     +        DSAM(IX),DSAM(IY),ISTAT)
        END
