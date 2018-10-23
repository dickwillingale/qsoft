*+QRI_SQBEAM Analyse source above background in a square beam
        SUBROUTINE QRI_SQBEAM(NELS1,NELS2,ARRAY,HBEAM,BLEV,BVAR,
     +   NP,XPI,XPR,YPI,YPR,BUF,
     +   NX,NY,BFLUX,BSIGMA,FLUX,FSIGMA,PEAK,CEN,RMSX,RMSY,
     +   PI5,PI25,MED,PI75,PI95,
     +   HEWX,HEWY,W90X,W90Y)
        IMPLICIT NONE
        INTEGER NELS1,NELS2,NP,NX,NY
        DOUBLE PRECISION XPI(NP),XPR(NP),YPI(NP),YPR(NP),BUF(NP)
        DOUBLE PRECISION ARRAY(NELS1,NELS2),HBEAM,BLEV,BVAR
        DOUBLE PRECISION BFLUX,BSIGMA,FLUX,FSIGMA,PEAK(2)
        DOUBLE PRECISION CEN(2),RMSX,RMSY
        DOUBLE PRECISION PI5(2),PI25(2),MED(2),PI75(2),PI95(2)
        DOUBLE PRECISION HEWX,HEWY,W90X,W90Y
Cf2py  intent(in) NELS1,NELS2,ARRAY,HBEAM,BLEV,BVAR,NP
Cf2py  intent(out) XPI,XPR,YPI,YPR,BUF
Cf2py  intent(out) NX,NY,BFLUX,BSIGMA,FLUX,FSIGMA,PEAK,CEN
Cf2py  intent(out) RMSX,RMSY,PI5,PI25,MED,PI75,PI95,HEWX,HEWY,W90X,W90Y
*NELS1       input        first dimension of array
*NELS2       input        second dimension of array
*ARRAY       input        image array
*HBEAM       input        half width of beam in pixels
*BLEV        input        average background level per pixel (to be subtracted)
*BVAR        input        variance on BLEV (-ve for counting statistics)
*NP          input        dimension of output arrays
*XPI,XPR     output       x output arrays, pixel position and flux
*YPI,YPR     output       y output arrays, pixel position and flux
*BUF         output       buffer array used to bin up cummulative distr.
*NX,NY       output       dimension of beam pixels
*BFLUX       output       background in beam (e.g. counts)
*BSIGMA      output       standard deviation of background
*FLUX        output       source flux above background in beam (e.g. counts)
*FSIGMA      output       standard deviation of source flux
*PEAK        output       source x,y peak position
*CEN         output       source x,y centroid position
*RMSX        output       rms width in x (pixels) about centroid
*RMSY        output       rms width in y (pixels) about centroid
*PI5         output       5% x,y position
*PI25        output       25% x,y position
*MED         output       median (50%) x,y position
*PI75        output       75% x,y position
*PI95        output       95% x,y position
*HEWX        output       half energy width x (pixels)
*HEWY        output       half energy width y (pixels)
*W90X        output       W90 (90% width) x (pixels)
*W90Y        output       W90 (90% width) y (pixels)
*Pixels assumed to run:
*        X 1 left to NELS1 right
*        Y 1 bottom to NELS2 top
*Coordinate system for XBEAM,YBEAM and other position is:
*        X runs from 0.0 on left to NELS1 on right
*        Y runs from 0.0 on bottom to NELS2 on top
*Centre of bottom left pixel is therefore 0.5,0.5
*Centre of top right pixel is NELS1-0.5,NELS2-0.5
*-Author Dick Willingale 2018-Oct-02
        INCLUDE 'QRI_TRANSCOM'
        DOUBLE PRECISION XBEAM,YBEAM,XY(2),P(2),TOTAL
        DOUBLE PRECISION PIBY2,PVAL,XSQ,YSQ,VAL
        DOUBLE PRECISION XMSQ,YMSQ
        INTEGER J,K,IXP,IYP,IRAD,NXL,NXH,NYL,NYH
        DOUBLE PRECISION P5,P25,P50,P75,P95
        DOUBLE PRECISION W5,W25,W50,W75,W95,RR,FRAC
        INCLUDE 'QR_COM'
C
        IF(ISTAT.NE.0) RETURN
C Get position
        XY(1)=XYNOW(1)
        XY(2)=XYNOW(2)
        CALL QRI_LTOP(XY,P)
        XBEAM=P(1)
        YBEAM=P(2)
C
         PIBY2=ASIN(1.D0)
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
        RMSX=0.0
        RMSY=0.0
        MED(1)=0.0
        MED(2)=0.0
        HEWX=0.0
        HEWY=0.0
        W90X=0.0
        W90Y=0.0
C Find pixel ranges about centre
        IXP=INT(XBEAM+1.0)
        IYP=INT(YBEAM+1.0)
        IF(IXP.LT.1.OR.IXP.GT.NELS1.OR.IYP.LT.1.OR.IYP.GT.NELS2) THEN
                RETURN
        ENDIF
        IRAD=INT(HBEAM)
        NXL=MAX(MIN(IXP-IRAD,NELS1),1)
        NXH=MAX(MIN(IXP+IRAD,NELS1),1)
        NYL=MAX(MIN(IYP-IRAD,NELS2),1)
        NYH=MAX(MIN(IYP+IRAD,NELS2),1)
        NX=NXH-NXL+1
        NY=NYH-NYL+1
        IF(NX+1.GT.NP.OR.NY+1.GT.NP) THEN
                WRITE(*,*) 'NX NY NP',NX,NY,NP
                WRITE(*,*) 'QRI_SQBEAM error - output arrays too small'
                RETURN
        ENDIF        
C Initialise internal stuff
        PVAL=-1.E32
        XSQ=0.0
        YSQ=0.0
C Loop round pixels
        DO J=NYL,NYH
                IYP=J-NYL+1
                YPI(IYP)=(DBLE(J)-0.5)
                DO K=NXL,NXH
                        IXP=K-NXL+1
                        XPI(IXP)=(DBLE(K)-0.5)
                        TOTAL=TOTAL+ARRAY(K,J)
                        VAL=ARRAY(K,J)-BLEV
                        FLUX=FLUX+VAL
                        IF(VAL.GT.PVAL) THEN
                                PVAL=VAL
                                PEAK(1)=XPI(IXP)
                                PEAK(2)=YPI(IYP)
                        ENDIF
                        XSQ=XSQ+VAL*(XPI(IXP)**2)
                        YSQ=YSQ+VAL*(YPI(IYP)**2)
                        CEN(1)=CEN(1)+VAL*XPI(IXP)
                        CEN(2)=CEN(2)+VAL*YPI(IYP)
                        XPR(IXP)=XPR(IXP)+VAL
                        YPR(IYP)=YPR(IYP)+VAL
                ENDDO
        ENDDO
C Normalise background to beam
        BFLUX=DBLE(NX*NY)*BLEV
        IF(BVAR.LT.0) THEN
                BSIGMA=SQRT(BFLUX)
        ELSEIF(BVAR.EQ.0) THEN
                BSIGMA=0.0
        ELSE
                BSIGMA=SQRT(DBLE(NX*NY)*BVAR)
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
C Calculate rms widths
        XSQ=XSQ/FLUX
        YSQ=YSQ/FLUX
        XMSQ=CEN(1)**2
        YMSQ=CEN(2)**2
        IF(XSQ.GT.XMSQ) THEN
                RMSX=SQRT(XSQ-XMSQ)
        ELSE
                RMSX=0.0
        ENDIF
        IF(YSQ.GT.YMSQ) THEN
                RMSY=SQRT(YSQ-YMSQ)
        ELSE
                RMSY=0.0
        ENDIF
C Bin up cummulative distribution in x
        BUF(1)=0.0
        BUF(2)=XPR(1)
        DO J=3,NX+1
                BUF(J)=BUF(J-1)+XPR(J-1)
        ENDDO
        W5=BUF(NX+1)*0.05
        P5=0.0
        W25=BUF(NX+1)*0.25
        P25=0.0
        W50=BUF(NX+1)*0.5
        P50=0.0
        W75=BUF(NX+1)*0.75
        P75=0.0
        W95=BUF(NX+1)*0.95
        P95=0.0
C Scan to find Median position, HEW and W90
        DO J=2,NX+1
                IF(BUF(J).GT.W5.AND.P5.EQ.0.0) THEN
                        RR=BUF(J)-BUF(J-1)
                        IF(RR.NE.0.0) THEN
                                FRAC=(W5-BUF(J-1))/RR
                        ELSE
                                FRAC=0.0
                        ENDIF
                        P5=(XPI(J-1)-0.5+FRAC)
                ENDIF
                IF(BUF(J).GT.W25.AND.P25.EQ.0.0) THEN
                        RR=BUF(J)-BUF(J-1)
                        IF(RR.NE.0.0) THEN
                                FRAC=(W25-BUF(J-1))/RR
                        ELSE
                                FRAC=0.0
                        ENDIF
                        P25=(XPI(J-1)-0.5+FRAC)
                ENDIF
                IF(BUF(J).GT.W50.AND.P50.EQ.0.0) THEN
                        RR=BUF(J)-BUF(J-1)
                        IF(RR.NE.0.0) THEN
                                FRAC=(W50-BUF(J-1))/RR
                        ELSE
                                FRAC=0.0
                        ENDIF
                        P50=(XPI(J-1)-0.5+FRAC)
                ENDIF
                IF(BUF(J).GT.W75.AND.P75.EQ.0.0) THEN
                        RR=BUF(J)-BUF(J-1)
                        IF(RR.NE.0.0) THEN
                                FRAC=(W75-BUF(J-1))/RR
                        ELSE
                                FRAC=0.0
                        ENDIF
                        P75=(XPI(J-1)-0.5+FRAC)
                ENDIF
                IF(BUF(J).GT.W95.AND.P95.EQ.0.0) THEN
                        RR=BUF(J)-BUF(J-1)
                        IF(RR.NE.0.0) THEN
                                FRAC=(W95-BUF(J-1))/RR
                        ELSE
                                FRAC=0.0
                        ENDIF
                        P95=(XPI(J-1)-0.5+FRAC)
                ENDIF
        ENDDO
        PI5(1)=P5
        PI25(1)=P25
        MED(1)=P50
        PI75(1)=P75
        PI95(1)=P95
        HEWX=P75-P25
        W90X=P95-P5
C Bin up cummulative distribution in y
        BUF(1)=0.0
        BUF(2)=YPR(1)
        DO J=3,NY+1
                BUF(J)=BUF(J-1)+YPR(J-1)
        ENDDO
        W5=BUF(NY+1)*0.05
        P5=0.0
        W25=BUF(NY+1)*0.25
        P25=0.0
        W50=BUF(NY+1)*0.5
        P50=0.0
        W75=BUF(NY+1)*0.75
        P75=0.0
        W95=BUF(NY+1)*0.95
        P95=0.0
C Scan to find Median position, HEW and W90
        DO J=2,NY+1
                IF(BUF(J).GT.W5.AND.P5.EQ.0.0) THEN
                        RR=BUF(J)-BUF(J-1)
                        IF(RR.NE.0.0) THEN
                                FRAC=(W5-BUF(J-1))/RR
                        ELSE
                                FRAC=0.0
                        ENDIF
                        P5=(YPI(J-1)-0.5+FRAC)
                ENDIF
                IF(BUF(J).GT.W25.AND.P25.EQ.0.0) THEN
                        RR=BUF(J)-BUF(J-1)
                        IF(RR.NE.0.0) THEN
                                FRAC=(W25-BUF(J-1))/RR
                        ELSE
                                FRAC=0.0
                        ENDIF
                        P25=(YPI(J-1)-0.5+FRAC)
                ENDIF
                IF(BUF(J).GT.W50.AND.P50.EQ.0.0) THEN
                        RR=BUF(J)-BUF(J-1)
                        IF(RR.NE.0.0) THEN
                                FRAC=(W50-BUF(J-1))/RR
                        ELSE
                                FRAC=0.0
                        ENDIF
                        P50=(YPI(J-1)-0.5+FRAC)
                ENDIF
                IF(BUF(J).GT.W75.AND.P75.EQ.0.0) THEN
                        RR=BUF(J)-BUF(J-1)
                        IF(RR.NE.0.0) THEN
                                FRAC=(W75-BUF(J-1))/RR
                        ELSE
                                FRAC=0.0
                        ENDIF
                        P75=(YPI(J-1)-0.5+FRAC)
                ENDIF
                IF(BUF(J).GT.W95.AND.P95.EQ.0.0) THEN
                        RR=BUF(J)-BUF(J-1)
                        IF(RR.NE.0.0) THEN
                                FRAC=(W95-BUF(J-1))/RR
                        ELSE
                                FRAC=0.0
                        ENDIF
                        P95=(YPI(J-1)-0.5+FRAC)
                ENDIF
        ENDDO
        PI5(2)=P5
        PI25(2)=P25
        MED(2)=P50
        PI75(2)=P75
        PI95(2)=P95
        HEWY=P75-P25
        W90Y=P95-P5
        END
