*+QRT_MATS Set deformation matrix
        SUBROUTINE QRT_MATS(ID,IM,KNX,KNY,XSAM,YSAM,ZDEF,NW,W1,W2)
*ID        input        deformation index
*IM        input        sub-matrix index
*KNX       input        number of x samples
*KNY       input        number of y samples
*XSAM      input        x values
*YSAM      input        y values
*ZDEF      input        deformation matrix
*NW        input        size of work arrays (at least MAX(KNX,KNY))
*W1        input        work array
*W2        input        work array
Cf2py  intent(in) ID,IM,KNX,KNY,XSAM,YSAM,ZDEF,NW,W1,W2
        IMPLICIT NONE
        INTEGER ID,IM
        INTEGER KNX,KNY,NW
        DOUBLE PRECISION XSAM(KNX),YSAM(KNY),ZDEF(KNX,KNY)
        DOUBLE PRECISION W1(NW),W2(NW)
*-Author Dick Willingale 2012-May-15
        INCLUDE 'SRT_COM'
        INTEGER IX,IY,IMZ,IMX,IMY
        INTEGER NM,NX,NY
        INCLUDE 'QR_COM'
C
        IF(ISTAT.NE.0) RETURN
C Check index in range
        IF(ID.LT.1.OR.ID.GT.MAXDF) THEN
         WRITE(*,*) 'QRT_MATS error - deformation index out of range',ID
         ISTAT=1
         RETURN
        ENDIF
C Check sub-matrix index in range
        IF(IM.LT.1.OR.IM.GT.IDFM(1,ID)) THEN
         WRITE(*,*) 'QRT_MATS error - sub-matrix index out of range',IM
         ISTAT=1        
         RETURN
        ENDIF
C Get indices from common
        NM=IDFM(1,ID)
        NX=IDFM(2,ID)
        NY=IDFM(3,ID)
        IX=IDFP(1,ID)
        IY=IDFP(2,ID)
        IMZ=IDFP(3,ID)
        IMX=IDFP(4,ID)
        IMY=IDFP(5,ID)
C Check array sizes
        IF(KNX.GT.NX.OR.KNY.GT.NY) THEN
                WRITE(*,*) 'QRT_MATS error - dimensions incorrect'
                ISTAT=1
        ENDIF
        IF(NW.LT.MAX(KNX,KNY)) THEN
                WRITE(*,*) 'QRT_MATS error - work arrays too small'
                ISTAT=1
        ENDIF
C Set matices
        CALL SRT_MATS(NM,IM,KNX,KNY,XSAM,YSAM,ZDEF,NW,W1,W2,
     +        DSAM(IX),DSAM(IY),DSAM(IMZ),DSAM(IMX),DSAM(IMY),ISTAT)
        END
