*+QRX_ISMTAU        X-ray optical depth of cold ISM gas
      subroutine qrx_ismtau(nh,z,ne,ear,tau,etau)
      implicit none
      INTEGER ne
      REAL nh,z,ear(0:ne),tau(ne),etau(ne)
Cf2py  intent(in) nh,z,ne,ear
Cf2py  intent(out) tau,etau
*nh    input        column density 10**21 cm-3 at redshift z
*z     input        redshift of source
*ne    input        number of energy bins
*ear   input        array of energy bin boundaries keV in observer frame
*tau   output       optical depth for each energy
*-Dick Willingale 2011-May-30
        real param(2)
        integer j
C set parameters for Xspec routine
        param(1)=nh*0.1
        param(2)=z
C Call Xpsec routine
        call xszphb(ear,ne,param,1,tau,etau)
C convert to optical depth
        do j=1,ne
                tau(j)=-log(tau(j))
        enddo
        end 
