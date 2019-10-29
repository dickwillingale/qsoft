*+QRA_GRBMSB        Minimum significance binning of GRB lightcurves
        SUBROUTINE QRA_GRBMSB(NC,NP,TIME,ARRAY,ERROR,SMIN,NMIN,NMAX,
     +        WRK1,WRK2,IWRK,
     +        NRAT,FRATE,EFRATE,RATES,ERATES,RTIME,ETIME,TLO,THI)
             IMPLICIT NONE
              INTEGER NC,NP
        DOUBLE PRECISION TIME(NP)
        DOUBLE PRECISION ARRAY(NC,NP),ERROR(NC,NP),SMIN
        DOUBLE PRECISION WRK1(NP),WRK2(NP)
        INTEGER IWRK(NP),NRAT,NMIN,NMAX
        DOUBLE PRECISION FRATE(NP),EFRATE(NP)
        DOUBLE PRECISION RATES(NC,NP),ERATES(NC,NP),RTIME(NP),ETIME(NP)
        DOUBLE PRECISION TLO,THI
*NC        input        number of energy channels
*NP        input        number of samples in light curve
*TIME        input        array of light curve times
*ARRAY        input        light curve samples
*ERROR        input        light curve errors
*SMIN        input        minimum significance for detection of flux
*NMIN        input        minimum number of input bins per output bin
*NMAX        input        maximum number of input bins per output bin
*WRK1        in/out        work array
*WRK2        in/out        work array
*IWRK        in/out        work array
*NRAT        output        number of rate samples
*FRATE        output        rate summed over channels
*EFRATE        output        errors on rate summed over channels
*RATES        output        array of count rates in each channel
*ERATES        output        array of errors on count rates in each channel
*RTIME        output        array of rate times
*ETIME        output        array of time errors (half time bin width)
*TLO        output        Start of significant flux (wrt trigger)
*THI        output        End of significant flux (wrt trigger)
*Adapted from original do_grbmsb.f routine
*-Author Dick Willingale 2013-Jul-05
        include 'QR_COM'
        INTEGER I,J,JJ,NN,J1
        DOUBLE PRECISION TBIN
C Check status
        IF(ISTAT.NE.0) RETURN
C
        TBIN=TIME(2)-TIME(1)
C Construct full lightcurve summing across all energy channels
        DO J=1,NP
                WRK1(J)=0.0
                WRK2(J)=0.0
                DO I=1,NC
                        WRK1(J)=WRK1(J)+ARRAY(I,J)
                        WRK2(J)=WRK2(J)+ERROR(I,J)**2
                ENDDO
        ENDDO
C Now do minimum significance binning
        CALL DO_MSB(NP,TIME,WRK1,WRK2,SMIN,NMIN,NMAX,WRK1,
     +        WRK2,NRAT,FRATE,
     +        EFRATE,RTIME,ETIME,IWRK,ISTAT)
C Find significant time span and clear output arrays
        TLO=1.e6
        THI=0.0
        DO J=1,NP
                IF(IWRK(J).GT.0.AND.TLO.EQ.1.e6) THEN
                        TLO=TIME(J)-TBIN
                ENDIF
                IF(IWRK(J).GT.0) THEN
                        THI=TIME(J)+TBIN
                ENDIF
                DO I=1,NC
                        RATES(I,J)=0.0
                        ERATES(I,J)=0.0
                ENDDO
                WRK1(J)=0.0
        ENDDO
C Bin up individual channels using binning
        J1=IWRK(1)
        JJ=1
        DO J=1,NP
                IF(IWRK(J).NE.J1) THEN
                        J1=IWRK(J)
                        JJ=JJ+1
                ENDIF
                IF(J1.NE.0) THEN
                        WRK1(JJ)=WRK1(JJ)+1.0
                        DO I=1,NC
                                RATES(I,JJ)=RATES(I,JJ)+ARRAY(I,J)
                                ERATES(I,JJ)=ERATES(I,JJ)+ERROR(I,J)**2
                        ENDDO
                ENDIF
        ENDDO
        DO J=1,NRAT
            IF(WRK1(J).GT.0.0) THEN
                DO I=1,NC
                        RATES(I,J)=RATES(I,J)/WRK1(J)
                        ERATES(I,J)=SQRT(ERATES(I,J))/WRK1(J)
                ENDDO
            ENDIF
        ENDDO
        END
