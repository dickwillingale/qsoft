*+QRB_BIOTEST        Map peaks intensities and calculate marker rates
        SUBROUTINE QRB_BIOTEST(NP,NS,ARRAY,FBOT,FMIN,FDIF,FPLU,FMAX,
     +  PDF,NMARK,
     +  IMARK,PHI,RATE,RMSDEV)
        INTEGER NP,NS,NMARK,IMARK(NMARK)
        DOUBLE PRECISION ARRAY(NP,NS),FMIN(NP),FDIF(NP),FPLU(NP)
        DOUBLE PRECISION FMAX(NP)
        DOUBLE PRECISION FBOT(NP),PHI(NP,NS),RATE(NS),PDF(NP)
        DOUBLE PRECISION RMSDEV(NS)
Cf2py    intent(in) NP,NS,ARRAY,FBOT,FMIN,FDIF,FPLU,FMAX,PDF,NMARK,IMARK
Cf2py    intent(inout) WRK1,WRK2
Cf2py    intent(out) PHI,RATE,RMSDEV
*NP        input        number of peaks
*NS        input        number of individuals
*ARRAY     input        flux values of peaks
*FBOT      input        lowest flux for peak
*FMIN      input        flux at lower disciminator
*FDIF      input        flux at maximum fractional difference
*FPLU      input        flux at upper disciminator
*FMAX      input        maximum flux for peak
*PDF       input        maximum cummulative difference (sign used on marker)
*NMARK     input        number of markers
*IMARK     input        index of markers
*PHI       output       mapped marker flux values
*RATE      output       marker rate for individuals
*RMSDEV    output       rms of marker rate for individuals
* QR version RW 2014-Feb-16
* qsoft version RW 2017-Jun-17
*-Author Dick Willingale 2005-Mar-05
        INCLUDE 'QR_COM'
        INTEGER J,K
        DOUBLE PRECISION WTS
C Check status
        IF(ISTAT.NE.0) RETURN
C Loop for all individuals
        DO J=1,NS
C Loop for all peaks
                DO K=1,NP
                        IF(ARRAY(K,J).LE.FMIN(K)) THEN
                                IF(FMIN(K).GT.FBOT(K)) THEN
                                     PHI(K,J)=(ARRAY(K,J)-FBOT(K))/
     +                               (FMIN(K)-FBOT(K))
                                     PHI(K,J)=0.5*PHI(K,J)-1.0
                                ELSE
                                        PHI(K,J)=-1.0
                                ENDIF
                        ELSEIF(ARRAY(K,J).LE.FDIF(K)) THEN
                                     PHI(K,J)=(ARRAY(K,J)-FMIN(K))/
     +                               (FDIF(K)-FMIN(K))
                                     PHI(K,J)=PHI(K,J)*0.5-0.5
                        ELSEIF(ARRAY(K,J).LE.FPLU(K)) THEN
                                     PHI(K,J)=(ARRAY(K,J)-FDIF(K))/
     +                               (FPLU(K)-FDIF(K))
                                     PHI(K,J)=PHI(K,J)*0.5
                        ELSE
                                     PHI(K,J)=(ARRAY(K,J)-FPLU(K))/
     +                               (FMAX(K)-FPLU(K))
                                     PHI(K,J)=PHI(K,J)*0.5+0.5
                        ENDIF
                        IF(PDF(K).EQ.0.0) THEN
                                PHI(K,J)=MAX(0.0,PHI(K,J))
                        ENDIF
                ENDDO
                RATE(J)=0.0
                RMSDEV(J)=0.0
C loop for marker peaks
                DO K=1,NMARK
                        KK=IMARK(K)
                        WTS=PHI(KK,J)*SIGN(1.0D0,PDF(KK))
C                        WTS=PHI(KK,J)*PDF(KK)
                        RATE(J)=RATE(J)+WTS
                        RMSDEV(J)=RMSDEV(J)+WTS**2
                ENDDO
        ENDDO
C Calculate rms deviation of rates
        DO J=1,NS
                RMSDEV(J)=RMSDEV(J)/NMARK-(RATE(J)/NMARK)**2
                IF(RMSDEV(J).GT.0.0) THEN
                        RMSDEV(J)=SQRT(RMSDEV(J))
                ENDIF
        ENDDO
        END
