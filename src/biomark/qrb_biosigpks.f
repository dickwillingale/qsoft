*+QRB_BIOSIGPKS        look for significant peaks (local maxima)
        SUBROUTINE QRB_BIOSIGPKS(SIG,NS,SI,N,MZ,FL,VA,WK,IP,PMZ,PFL,
     +  PSI,NP)
        IMPLICIT NONE
        INTEGER NS,N,NP
        INTEGER IP(N)
        DOUBLE PRECISION SIG,SI(NS)
        DOUBLE PRECISION MZ(N),FL(N),VA(N),WK(N),PMZ(N),PFL(N),PSI(N)
Cf2py    intent(in) SIG,NS,SI,N,MZ,FL,VA
Cf2py    intent(inout) WK
Cf2py    intent(out) IP,PMZ,PFL,PSI,NP
*SIG  input     minimum significance for peaks
*NS   input     number of gaussian sigma values 1 or N
*SI   input     array of gaussian sigma values (in units of samples)
*N    input     number of data points
*MZ   input     array of data mass values
*FL   input     array of data flux values
*VA   input     array of data variance values
*WK   input     work array (size N)
*IP   output    indices of peaks found
*PMZ  output    mass values for peak positions
*PFL  output    flux values for peaks 
*PSI  output    significance values for peaks
*NP   output    number of peaks found
* qsoft version RW 2017-Jun-26
*-Author Dick Willingale 2014-Mar-31
        INTEGER K,K1,K2
C Smooth fluxes
        CALL QRB_BIOSPREAD(N,FL,WK,NS,SI,PFL)
C Calculate and smooth weighted fluxes
        DO K=1,N
                PSI(K)=FL(K)*MZ(K)
        ENDDO
        CALL QRB_BIOSPREAD(N,PSI,WK,NS,SI,PMZ)
C Calculate centroids of potential peaks
        DO K=1,N
                PMZ(K)=PMZ(K)/PFL(K)
        ENDDO
C Smooth variances
        CALL QRB_BIOSPREAD(N,VA,WK,NS,SI,PSI)
C Calculate significance and set mask
        DO K=1,N
                PSI(K)=PFL(K)/SQRT(PSI(K))
                IP(K)=0
        ENDDO
C Hunt for significant local maxima and collect together
        NP=0
        DO K=2,N-1
                K1=K-1
                K2=K+1
                IF((PFL(K).GT.PFL(K1)).AND.(PFL(K).GT.PFL(K2))
     +                .AND.(PSI(K).GT.SIG)) THEN
                        NP=NP+1
                        PFL(NP)=PFL(K)
                        PMZ(NP)=PMZ(K)
                        PSI(NP)=PSI(K)
                        IP(NP)=K
                ENDIF
        ENDDO
        END
