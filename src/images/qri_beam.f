*+QRI_BEAM Analyse source above background in a circular beam
        SUBROUTINE QRI_BEAM(NELS1,NELS2,ARRAY,RBEAM,BLEV,BVAR,
     +        NSAM,BFLUX,BSIGMA,FLUX,FSIGMA,PEAK,CEN,THA,RMSA,RMSB,
     +        FWHM,HEW,W90,FWHMP,HEWP,W90P,FWHMC,HEWC,W90C)
        IMPLICIT NONE
        INTEGER NELS1,NELS2,NSAM
        DOUBLE PRECISION ARRAY(NELS1,NELS2),RBEAM,BLEV,BVAR
        DOUBLE PRECISION BFLUX,BSIGMA,FLUX,FSIGMA,PEAK(2)
        DOUBLE PRECISION CEN(2),THA,RMSA,RMSB,FWHM,HEW,W90
        DOUBLE PRECISION FWHMP,HEWP,W90P,FWHMC,HEWC,W90C
Cf2py  intent(in) NELS1,NELS2,ARRAY,RBEAM,BLEV,BVAR
Cf2py  intent(out) NSAM,BFLUX,BSIGMA,FLUX,FSIGMA,PEAK,CEN,THA
Cf2py  intent(out) RMSA,RMSB,FWHM,HEW,W90,FWHMP,HEWP,W90P
Cf2py  intent(out) FWHMC,HEWC,W90C
*NELS1       input        first dimension of array
*NELS2       input        second dimension of array
*ARRAY       input        image array
*RBEAM       input        radius of beam in pixels
*BLEV        input        average background level per pixel (to be subtracted)
*BVAR        input        variance on BLEV (-ve for counting statistics)
*NSAM        output       number of pixels in beam
*BFLUX       output       background in beam (e.g. counts)
*BSIGMA      output       standard deviation of background
*FLUX        output       source flux above background in beam (e.g. counts)
*FSIGMA      output       standard deviation of source flux
*PEAK        output       source x,y peak position
*CEN         output       source x,y centroid position
*THA         output       angle (degrees) of major axis wrt x (x to y +ve)
*RMSA        output       source max rms width (major axis) (pixels)
*RMSB        output       source min rms width (minor axis) (pixels)
*FWHM        output       full width half maximum (pixels) about beam centre
*HEW         output       half energy width (pixels) about beam centre
*W90         output       W90 (90% width) (pixels) about beam centre
*FWHMP       output       full width half maximum (pixels) about peak
*HEWP        output       half energy width (pixels) about peak
*W90P        output       W90 (90% width) (pixels) about peak
*FWHMC       output       full width half maximum (pixels) about centroid
*HEWC        output       half energy width (pixels) about centroid
*W90C        output       W90 (90% width) (pixels) about centroid
*Pixels assumed to run:
*        X 1 left to NELS1 right
*        Y 1 bottom to NELS2 top
*Coordinate system for XBEAM,YBEAM and other position is:
*        X runs from 0.0 on left to NELS1 on right
*        Y runs from 0.0 on bottom to NELS2 on top
*Centre of bottom left pixel is therefore 0.5,0.5
*Centre of top right pixel is NELS1-0.5,NELS2-0.5
*-Author Dick Willingale 2012-May-14
        INCLUDE 'QRI_TRANSCOM'
        DOUBLE PRECISION XBEAM,YBEAM,XY(2),P(2),TOTAL
        DOUBLE PRECISION PIBY2,PVAL,XSQ,YSQ,XYS,YP,XP,RAD,VAL
        DOUBLE PRECISION XMSQ,YMSQ,XXX,TEMP,RMSX,RMSY
        INTEGER J,K,IXP,IYP,IRAD,NXL,NXH,NYL,NYH
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
        THA=0.0
        RMSA=0.0
        RMSB=0.0
        FWHM=0.0
        HEW=0.0
        W90=0.0
        NSAM=0
C Find pixel ranges about centre
        IXP=INT(XBEAM+1.0)
        IYP=INT(YBEAM+1.0)
        IF(IXP.LT.1.OR.IXP.GT.NELS1.OR.IYP.LT.1.OR.IYP.GT.NELS2) THEN
                RETURN
        ENDIF
        IRAD=INT(RBEAM+2.0)
        NXL=MAX(MIN(IXP-IRAD,NELS1),1)
        NXH=MAX(MIN(IXP+IRAD,NELS1),1)
        NYL=MAX(MIN(IYP-IRAD,NELS2),1)
        NYH=MAX(MIN(IYP+IRAD,NELS2),1)
C Initialise internal stuff
        PVAL=-1.E32
        XSQ=0.0
        YSQ=0.0
        XYS=0.0
C Loop round pixels
        DO J=NYL,NYH
                YP=(DBLE(J)-0.5)
                DO K=NXL,NXH
                        XP=(DBLE(K)-0.5)
                        RAD=SQRT((XP-XBEAM)**2+(YP-YBEAM)**2)
                        IF(RAD.LE.RBEAM) THEN
                                NSAM=NSAM+1
                                TOTAL=TOTAL+ARRAY(K,J)
                                VAL=ARRAY(K,J)-BLEV
                                FLUX=FLUX+VAL
                                IF(VAL.GT.PVAL) THEN
                                        PVAL=VAL
                                        PEAK(1)=XP
                                        PEAK(2)=YP
                                ENDIF
                                XSQ=XSQ+VAL*(XP**2)
                                YSQ=YSQ+VAL*(YP**2)
                                XYS=XYS+VAL*YP*XP
                                CEN(1)=CEN(1)+VAL*XP
                                CEN(2)=CEN(2)+VAL*YP
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
C Calculate rms widths
        XSQ=XSQ/FLUX
        YSQ=YSQ/FLUX
        XYS=XYS/FLUX
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
        XYS=XYS-CEN(1)*CEN(2)
C Calculate angle of stationary value (maximum width) wrt angle
        XXX=RMSX**2-RMSY**2
        IF(XYS.NE.0.) THEN
                THA=0.5*ATAN2(2.*XYS,XXX)
        ELSE
                THA=0.
        ENDIF
C Calculate major and minor widths
        RMSA=(RMSX*COS(THA))**2+(RMSY*SIN(THA))**2+
     +        2.*XYS*COS(THA)*SIN(THA)
        IF(RMSA.GT.0) THEN
                RMSA=SQRT(RMSA)
        ELSE
                RMSA=0.0
        ENDIF
        RMSB=(RMSX*SIN(THA))**2+(RMSY*COS(THA))**2-
     +        2.*XYS*COS(THA)*SIN(THA)
        IF(RMSB.GT.0.0) THEN
                RMSB=SQRT(RMSB)
        ELSE
                RMSB=0.0
        ENDIF
        IF(RMSA.LT.RMSB) THEN
C If found minor axis then swap
                TEMP=RMSA
                RMSA=RMSB
                RMSB=TEMP
                THA=THA+PIBY2
        ENDIF
C Force angular range of major axis to be 0 - pi
        IF(THA.LT.0.) THEN
                THA=THA+2.*PIBY2
        ELSEIF(THA.GT.2.*PIBY2) THEN
                THA=THA-2.*PIBY2
        ENDIF
C Convert angle of major axis to degrees
        THA=THA*90./PIBY2
C Find FWHM, HEW and W90 about beam centre
        CALL QRI_WBEAM(NELS1,NELS2,ARRAY,NXL,NXH,NYL,NYH,
     +        BLEV,XBEAM,YBEAM,RBEAM,HEW,FWHM,W90)
C Find FWHM and HEW and W90 about peak
        CALL QRI_WBEAM(NELS1,NELS2,ARRAY,NXL,NXH,NYL,NYH,
     +        BLEV,PEAK(1),PEAK(2),RBEAM,HEWP,FWHMP,W90P)
C Find FWHM and HEW and W90 about centroid
        CALL QRI_WBEAM(NELS1,NELS2,ARRAY,NXL,NXH,NYL,NYH,
     +        BLEV,CEN(1),CEN(2),RBEAM,HEWC,FWHMC,W90C)
        END
*+QRI_WBEAM        Calculate HEW and FWHM of source in beam about given centre
        SUBROUTINE QRI_WBEAM(NELS1,NELS2,ARRAY,NXL,NXH,NYL,NYH,BLEV,
     +        XCEN,YCEN,RBEAM,HEW,FWHM,W90)
        IMPLICIT NONE
        INTEGER NELS1,NELS2,NXL,NXH,NYL,NYH
        DOUBLE PRECISION  ARRAY(NELS1,NELS2),BLEV,XCEN,YCEN,RBEAM
        DOUBLE PRECISION  HEW,FWHM,W90
*NELS1       input        dimension of array
*NELS2       input        dimension of array
*ARRAY       input        data array
*NXL,NXH     input        x pixel range (column range) of box
*NYL,NYH     input        y pixel range (row range) of box
*BLEV        input        average background level (to be subtracted)
*XCEN        input        x centre
*YCEN        input        y centre
*RBEAM       input        radius of beam
*HEW         output        half energy width
*FWHM        output        full width half maximum
*W90         output        90% width
*-Author Dick Willingale 2012-May-14
        INTEGER NBUF,IRAD,J,K,NRAD
        PARAMETER (NBUF=1000)
        DOUBLE PRECISION BUF(NBUF),BUF1(NBUF)
        DOUBLE PRECISION XP,YP,DMAX,FRAC,HMAX,RAD,WTOT,RR,HTOT
C Initialize buffers
        DO J=1,NBUF
                BUF(J)=0.0
                BUF1(J)=0.0
        ENDDO
C Bin up about centre
        NRAD=0
        DO J=NYL,NYH
                YP=(DBLE(J)-0.5)-YCEN
                DO K=NXL,NXH
                        XP=(DBLE(K)-0.5)-XCEN
                        RAD=SQRT(XP**2+YP**2)
                        IF(RAD.LE.RBEAM) THEN
                                IRAD=INT(RAD+1.0)
                                IRAD=MIN(IRAD,NBUF)
                                IRAD=MAX(IRAD,1)
                                NRAD=MAX(IRAD,NRAD)
                                BUF(IRAD)=BUF(IRAD)+ARRAY(K,J)-BLEV
                                BUF1(IRAD)=BUF1(IRAD)+1.0
                        ENDIF
                ENDDO
        ENDDO
C Normalize to surface brightness and find maximum
        DMAX=0.0
        DO J=1,NRAD
                IF(BUF1(J).GT.0) THEN
                        BUF1(J)=BUF(J)/BUF1(J)
                        DMAX=MAX(DMAX,BUF1(J))
                ENDIF
        ENDDO
C Find FWHM
        HMAX=DMAX*0.5
        DO J=NRAD,1,-1
                IF(BUF1(J).GT.HMAX) THEN
                        RR=BUF1(J)-BUF1(J+1)
                        IF(RR.NE.0.0) THEN
                                FRAC=(HMAX-BUF1(J+1))/RR
                        ELSE
                                FRAC=0.0
                        ENDIF
                        FWHM=((DBLE(J+1)-0.5)-FRAC)*2.0
                        GOTO 100
                ENDIF
        ENDDO
100        CONTINUE
C Bin up ready to find HEW
        BUF1(1)=0.0
        BUF1(2)=BUF(1)
        DO J=3,NRAD+1
                BUF1(J)=BUF1(J-1)+BUF(J-1)
        ENDDO
        HTOT=BUF1(NRAD+1)*0.5
        WTOT=BUF1(NRAD+1)*0.9
C Scan to find HEW and W90
        DO J=2,NRAD+1
                IF(BUF1(J).GT.HTOT.AND.HEW.EQ.0.0) THEN
                        RR=BUF1(J)-BUF1(J-1)
                        IF(RR.NE.0.0) THEN
                                FRAC=(HTOT-BUF1(J-1))/RR
                        ELSE
                                FRAC=0.0
                        ENDIF
                        HEW=(DBLE(J-2)+FRAC)*2.0
                ENDIF
                IF(BUF1(J).GT.WTOT.AND.W90.EQ.0.0) THEN
                        RR=BUF1(J)-BUF1(J-1)
                        IF(RR.NE.0.0) THEN
                                FRAC=(WTOT-BUF1(J-1))/RR
                        ELSE
                                FRAC=0.0
                        ENDIF
                        W90=(DBLE(J-2)+FRAC)*2.0
                        RETURN
                ENDIF
        ENDDO
        END
