*+QRB_BIOSCALE        Scale peaks from all samples
        SUBROUTINE QRB_BIOSCALE(NPN,INORM,NS,NP,CFLX,PFLX,
     +  PVAR,PBAK,
     +  RNORM,IMAX,PCTS,PSIG,PND)
        IMPLICIT NONE
        INTEGER NS,NP,IMAX(NP),NPN,INORM(NPN)
        DOUBLE PRECISION CFLX(NP)
        DOUBLE PRECISION PFLX(NP,NS),PVAR(NP,NS),PBAK(NP,NS)
        DOUBLE PRECISION RNORM(NS)
        DOUBLE PRECISION PCTS(NP),PSIG(NP),PND(NP,NS)
Cf2py    intent(in) NPN,INORM,NS,NP,CFLX
Cf2py    intent(inout) PFLX,PVAR,PBAK
Cf2py    intent(out) RNORM,IMAX,PCTS,PSIV,PND
*NPN       input        number of normalisation peaks
*INORM     input        index of normalisation peaks
*NS        input        number of individuals
*NP        input        number of peaks
*CFLX      input        combined mean flux for each peak
*PFLX      in/out       peak fluxes
*PFLX      in/out       peak variances
*PBAK      in/out       peak background
*RNORM     output       normalisation used for each instance
*IMAX      output       instance index of maximum peak value
*PCTS      output       mean flux in peak after normalisation
*PSIG      output       fractional rms variation in peaks
*PND       output       difference between normalised peaks and average spectrum
* QR version RW 2014-Feb-16
* qsoft version RW 2017-Jun-26
*-Author Dick Willingale 2005-Apr-21
        INCLUDE 'QR_COM'
C Check status
        IF(ISTAT.NE.0) RETURN
C Normalise
        CALL QRB_BIONORM(NS,NP,CFLX,NPN,INORM,RNORM,PFLX,
     +  PVAR,PBAK)
C calculate mean and scatter of peaks
        CALL QRB_BIOSCAT(NS,NP,CFLX,PFLX,PCTS,PSIG,IMAX,PND)
        END
