*+QRA_GRBLOADPARS Load GRB light curve model parameters into common
        SUBROUTINE QRA_GRBLOADPARS(NP,PAR,NG,GRBN)
        IMPLICIT NONE
        INTEGER NP,NG
        DOUBLE PRECISION PAR(NP)
        CHARACTER*(NG) GRBN
*NP        input        number of parameters
*PAR        input        model parameters (see routine qra_grbrates for details)
*NG     input   number of characters in tabulation file name
*GRBN   input   GRB name
*-Author Dick Willingale 2013-Jan-24
        include 'QR_COM'
        include 'GRBS_COM'
        integer J
C
        IF(ISTAT.NE.0) RETURN
C
        IF(NP.GT.maxp) THEN
                write(*,*) 'QRA_GRBLOADPARS too many fit parameters'
                ISTAT=1
                RETURN
        ENDIF
C Set GRB name in common
        grbname=GRBN(1:NG)
C Reset count in common so that qra_specrates() re-loads look up table
        na=0
C Set fit parameters in common
        npars=NP
        DO J=1,NP
                pars(J)=PAR(J)
        ENDDO
        END
