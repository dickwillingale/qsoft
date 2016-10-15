*+QRX_IISMTAU        X-ray optical depth of ionized ISM gas
      subroutine qrx_iismtau(nh,z,tk,pl,ist,ne,ear,tau,etau)
      implicit none
      INTEGER ne
      REAL nh,z,tk,pl,ist,ear(0:ne),tau(ne),etau(ne)
Cf2py  intent(in) nh,z,tk,pl,ist,ne,ear
Cf2py  intent(out) tau,etau
*nh    input        column density 10**21 cm-3 at redshift z
*z     input        redshift of source
*tk    input        temperature Kelvin
*pl    input        powerlaw index
*ist   input        ionization state
*ne    input        number of energy bins
*ear   input        array of energy bin boundaries keV in observer frame
*tau   output       optical depth for each energy
*-Dick Willingale 2011-May-30
        real param(6)
C set parameters for Xspec routine
        param(1)=pl
        param(2)=nh*0.1
        param(3)=tk
        param(4)=ist
        param(5)=z
        param(6)=1.0
        call xsabsori_tau(ear,ne,param,1,tau,etau)
        end 
