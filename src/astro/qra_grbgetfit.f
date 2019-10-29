*+QRA_GRBGETFIT get parameters and model light curves of fit
        SUBROUTINE QRA_GRBGETFIT(NR,NT,NP,TIME,RATS,PAR,NCL)
        IMPLICIT NONE
        INTEGER NR,NT,NP,NCL
        DOUBLE PRECISION TIME(NT),RATS(NR,NT),PAR(NP)
*NR        input        number of rates (8 channels)
*NT        input        number of time samples
*NP        input        total number of fit parameters 
*TIME        input        array of times
*RATS        output        array of rates
*PAR        outout        array of parameter values
*NCL        output        number of calls to qra_grbchisq routine
*-Author Dick Willingale 2013-Jan-24
        include 'QR_COM'
        include 'GRBS_COM'
        integer J,I,ipul
        double precision rates(8),flux,elo,ehi
C
        ipul=-1
        elo=0.3
        ehi=350.0
C
        IF(ISTAT.NE.0) RETURN
C
        IF(NP.ne.npars) THEN
                write(*,*) 'QRA_GRBGETFIT - wrong number of parameters'
                ISTAT=1
                RETURN
        ENDIF
        IF(NR.ne.8) THEN
                write(*,*) 'QRA_GRBGETFIT - wrong number of channels'
                ISTAT=1
                RETURN
        ENDIF
C Get rates
        DO J=1,NT
                CALL qra_grbrates(TIME(j),npars,pars,ipul,
     +          elo,ehi,8,rates,flux)
                DO I=1,NR
                        RATS(I,J)=rates(I)
                ENDDO
        ENDDO
C Get parameters values from common
        DO J=1,NP
                PAR(J)=pars(J)
        ENDDO
C
        NCL=ncall
        END
