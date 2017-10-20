*+QRI_LECBEAM Analyse source above background in a lobster eye cross beam
        SUBROUTINE QRI_LECBEAM(NELS1,NELS2,ARRAY,S,H,BLEV,BVAR,NT,
     +        QUA,QUAN,NSAM,BFLUX,BSIGMA,FLUX,FSIGMA,PEAK,CEN,HEW,W90,
     +        AHEW,AW90,FPEAK)
        IMPLICIT NONE
        INTEGER NELS1,NELS2,NSAM,NT
        DOUBLE PRECISION QUA(NT,NT),QUAN(NT,NT)
        DOUBLE PRECISION ARRAY(NELS1,NELS2),S,H,BLEV,BVAR
        DOUBLE PRECISION BFLUX,BSIGMA,FLUX,FSIGMA,PEAK(2)
        DOUBLE PRECISION CEN(2),HEW,W90,AHEW,AW90,FPEAK
Cf2py  intent(in) NELS1,NELS2,ARRAY,S,H,BLEV,BVAR,NT
Cf2py  intent(out) QUA,QUAN,NSAM,BFLUX,BSIGMA,FLUX,FSIGMA,PEAK,CEN
Cf2py  intent(out) HEW,W90,AHEW,AW90,FPEAK
*NELS1       input        first dimension of array
*NELS2       input        second dimension of array
*ARRAY       input        image array
*S           input        size of square area in pixels
*H           input        height of cross-arm quadrant in pixels (=2d/L)
*BLEV        input        average background level per pixel (to be subtracted)
*BVAR        input        variance on BLEV (-ve for counting statistics)
*NT          input        dimension of output quadrant flux distribution
*QUA         output       quadrant surface brightness distribution
*QUAN        output       quadrant  pixel occupancy
*NSAM        output       number of pixels in beam
*BFLUX       output       background in beam (e.g. counts)
*BSIGMA      output       standard deviation of background
*FLUX        output       source flux above background in beam (e.g. counts)
*FSIGMA      output       standard deviation of source flux
*PEAK        output       source x,y peak position
*CEN         output       source x,y centroid position
*HEW         output       half energy width (pixels)
*W90         output       W90 (90% width) (pixels)
*AHEW        output       half energy area (sq pixels)
*AW90        output       W90 (90% width) area (sq pixels)
*FPEAK       output       flux in peak pixel
*Pixels assumed to run:
*        X 1 left to NELS1 right
*        Y 1 bottom to NELS2 top
*Coordinate system for XBEAM,YBEAM and other position is:
*        X runs from 0.0 on left to NELS1 on right
*        Y runs from 0.0 on bottom to NELS2 on top
*Centre of bottom left pixel is therefore 0.5,0.5
*Centre of top right pixel is NELS1-0.5,NELS2-0.5
*-Author Dick Willingale 2017-Sept-14 based on QRI_BEAM routine
        INCLUDE 'QRI_TRANSCOM'
        DOUBLE PRECISION XBEAM,YBEAM,XY(2),P(2),TOTAL
        DOUBLE PRECISION YP,XP,XR,YR,VAL
        DOUBLE PRECISION S2,BB2,BH,HTOT,WTOT
        INTEGER I,J,K,IXP,IYP,IRAD,NXL,NXH,NYL,NYH,NBUF
        PARAMETER (NBUF=1000)
        DOUBLE PRECISION BUF(NBUF),BUF1(NBUF),FRAC,RR
        INCLUDE 'QR_COM'
C
        IF(ISTAT.NE.0) RETURN
C Get position
        XY(1)=XYNOW(1)
        XY(2)=XYNOW(2)
        CALL QRI_LTOP(XY,P)
        XBEAM=P(1)
        YBEAM=P(2)
C Initialise numbers to be determined
        TOTAL=0.0
        BFLUX=0.0
        BSIGMA=0.0
        FLUX=0.0
        FSIGMA=0.0
        PEAK(1)=0.0
        PEAK(2)=0.0
        CEN(1)=0.0
        CEN(2)=0.0
        HEW=0.0
        W90=0.0
        NSAM=0
        DO J=1,NT
           DO K=1,NT
               QUA(K,J)=0.0
               QUAN(K,J)=0.0
           ENDDO
        ENDDO
C Find pixel ranges about centre
        IXP=INT(XBEAM+1.0)
        IYP=INT(YBEAM+1.0)
        IF(IXP.LT.1.OR.IXP.GT.NELS1.OR.IYP.LT.1.OR.IYP.GT.NELS2) THEN
                RETURN
        ENDIF
        S2=S*0.5
        IRAD=NINT(MIN(S2,H))
        NXL=MAX(MIN(IXP-IRAD,NELS1),1)
        NXH=MAX(MIN(IXP+IRAD,NELS1),1)
        NYL=MAX(MIN(IYP-IRAD,NELS2),1)
        NYH=MAX(MIN(IYP+IRAD,NELS2),1)
C Initialise internal stuff
        FPEAK=-1.E32
C Initialize buffers
        DO I=1,NBUF
                BUF(I)=0.0
                BUF1(I)=0.0
        ENDDO
C Loop round pixels
        DO J=NYL,NYH
                YP=(DBLE(J)-0.5)
                YR=ABS(YP-YBEAM)
                DO K=NXL,NXH
                        XP=(DBLE(K)-0.5)
                        XR=ABS(XP-XBEAM)
                        IF((XR.LT.S2).AND.(YR.LT.S2)) THEN
                            NSAM=NSAM+1
                            TOTAL=TOTAL+ARRAY(K,J)
                            VAL=ARRAY(K,J)-BLEV
                            FLUX=FLUX+VAL
                            IF(VAL.GT.FPEAK) THEN
                                FPEAK=VAL
                                PEAK(1)=XP
                                PEAK(2)=YP
                            ENDIF
                            CEN(1)=CEN(1)+VAL*XP
                            CEN(2)=CEN(2)+VAL*YP
                            DO I=1,NBUF
                                BB2=I
                                BH=BB2/H
                                IF((YR.LT.(BB2-BH*XR)).OR.
     +                          (XR.LT.(BB2-BH*YR))) THEN
                                    BUF(I)=BUF(I)+VAL
                                    BUF1(I)=BUF1(I)+1.0
                                ENDIF
                            ENDDO
                            IXP=INT(XR)+1
                            IYP=INT(YR)+1
                            IF(IXP.LE.NT.AND.IYP.LE.NT) THEN
                                QUA(IXP,IYP)=QUA(IXP,IYP)+VAL
                                QUAN(IXP,IYP)=QUAN(IXP,IYP)+1.0
                            ENDIF
                        ENDIF
                ENDDO
        ENDDO
C Normalise background to beam
        BFLUX=DBLE(NSAM)*BLEV
        IF(BVAR.LT.0) THEN
                BSIGMA=SQRT(BFLUX)
        ELSEIF(BVAR.EQ.0) THEN
                BSIGMA=0.0
        ELSE
                BSIGMA=SQRT(DBLE(NSAM)*BVAR)
        ENDIF
        IF(FLUX.LE.0.0) THEN
                RETURN
        ENDIF
C Calculate error on flux
        IF(BVAR.LT.0) THEN
                FSIGMA=SQRT(TOTAL+BFLUX)
        ELSEIF(BVAR.EQ.0) THEN
                FSIGMA=0.0
        ELSE
                FSIGMA=BSIGMA
        ENDIF
C Calculate centroid position
        CEN(1)=CEN(1)/FLUX
        CEN(2)=CEN(2)/FLUX
C Scan buffer to find HEW and W90
        HTOT=FLUX*0.5
        WTOT=FLUX*0.9
        DO I=2,NBUF
            IF(BUF(I).GT.HTOT.AND.HEW.EQ.0.0) THEN
                RR=BUF(I)-BUF(I-1)
                IF(RR.NE.0.0) THEN
                   FRAC=(HTOT-BUF(I-1))/RR
                ELSE
                   FRAC=0.0
                ENDIF
                HEW=(DBLE(I-1)+FRAC)
                AHEW=BUF1(I-1)+FRAC*(BUF1(I)-BUF1(I-1))
            ENDIF
            IF(BUF(I).GT.WTOT.AND.W90.EQ.0.0) THEN
                RR=BUF(I)-BUF(I-1)
                IF(RR.NE.0.0) THEN
                    FRAC=(WTOT-BUF(I-1))/RR
                ELSE
                    FRAC=0.0
                ENDIF
                W90=(DBLE(I-1)+FRAC)
                AW90=BUF1(I-1)+FRAC*(BUF1(I)-BUF1(I-1))
            ENDIF
        ENDDO
C Correct HEW and W09 to side of central square
        HEW=2.0*(HEW-HEW**2/H)/(1.0-HEW**2/H**2)
        W90=2.0*(W90-W90**2/H)/(1.0-W90**2/H**2)
C Normalise quadrant to surface brightness
        DO J=1,NT
           DO K=1,NT
               IF(QUAN(K,J).GT.0.0) THEN
                   QUA(K,J)=QUA(K,J)/QUAN(K,J)
               ENDIF
           ENDDO
        ENDDO
        END
