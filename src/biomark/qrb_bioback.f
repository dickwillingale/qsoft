*+QRB_BIOBACK        Continuum subtraction from mass spectra
        SUBROUTINE QRB_BIOBACK(NP,NS,ARRAY,NBACK,WK1,WK2,NFRAC,
     +        PEAKS,BAK,VAR)
        INTEGER NP,NS,NBACK,NFRAC
        DOUBLE PRECISION ARRAY(NP,NS),WK1(NBACK),WK2(NBACK)
        DOUBLE PRECISION PEAKS(NP,NS),BAK(NP,NS),VAR(NP,NS)
Cf2py   intent(in) NP,NS,ARRAY,NBACK,WK1,WK2,NFRAC
Cf2py   intent(out) PEAKS,BAK,VAR
*NP        input        number of peaks
*NS        input        number of individuals
*ARRAY     input        weights values of peaks range -1 to +1
*NBACK     input        size of background window
*WK1       input        work array
*WK2       input        work array
*NFRAC     input        number of samples in window below background level
*PEAKS     output       peaks above background
*BAK       output       continuum
*VAR       output       variance
* QR version RW 2014-Feb-16
* qsoft version RW 2017-Jun-26
*-Author Dick Willingale 2005-Mar-05
        INCLUDE 'QR_COM'
        INTEGER J,K,I,N1,N2,NH,NW,IT,ITT
C Check status
        IF(ISTAT.NE.0) RETURN
C set half window width
        NH=MAX(NBACK/2,1)
C Loop for all individual spectra
        DO J=1,NS
C Loop for all bins in spectrum
                DO K=1,NP
C set range of window at this position
                        N1=MAX(1,K-NH+1)
                        N2=MIN(NP,K+NH-1)
                        NW=N2-N1+1
                        DO I=N1,N2
                                WK1(I-N1+1)=ARRAY(I,J)
                                WK2(I-N1+1)=I
                        ENDDO
C sort window samples into ascending order
                        CALL PDA_DSORT(WK1,WK2,NW,2,ISTAT)
C calculate background from low tail
                        IT=MIN(NFRAC,NW)
                        BAK(K,J)=WK1(IT)
                        PEAKS(K,J)=ARRAY(K,J)-BAK(K,J)
C Estimate variance about background level ignoring top fraction
                        ITT=NW-NFRAC
                        VAR(K,J)=0.0
                        DO I=1,ITT
                                VAR(K,J)=VAR(K,J)+(WK1(I)-WK1(IT))**2
                        ENDDO
                        VAR(K,J)=VAR(K,J)/ITT
                ENDDO
        ENDDO
        END
