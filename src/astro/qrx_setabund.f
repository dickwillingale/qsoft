*-QRX_SETABUND        Set abundances for XSPEC routines
        SUBROUTINE QRX_SETABUND(ABUN,AMET)
        implicit none
        CHARACTER*(30) ABUN
        REAL AMET
Cf2py  intent(in) abun,amet
*ABUN       input        XSPEC abundance table
*AMET       input        metalicity
*-Author Dick Willingale 2011-May-26
        include 'SPX_COM'
        INTEGER ILEN,LEN_TRIM,i
        EXTERNAL LEN_TRIM
C Get parameters
        ILEN=LEN_TRIM(ABUN)
        metalicity=AMET
        iabu=0
        do i=1,nabu
             if(abundance(i) .eq. ABUN(1:ILEN)) then
                     iabu=i
             endif
        enddo
        if(iabu.eq.0) then
             write(*,*) 'qrx_setabund error: unkown abundance ',abun
        endif
        END
