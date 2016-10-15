*+QRX_IIGMTAU  X-ray optical depth of ionized IGM gas
        subroutine qrx_iigmtau(n0,pl,tk,ist,zmax,h0,omegam,omegal,
     +  ne,ear,xsec,tau)
        implicit none
        INTEGER ne
        REAL n0,pl,tk,ist,zmax,h0,omegam,omegal
        REAL ear(0:ne),xsec(ne),tau(ne)
Cf2py  intent(in) n0,pl,tk,ist,zmax,h0,omegam,omegal,ne,ear
Cf2py  intent(out) xsec,tau
*n0      input        number density cm-3 at z=0
*pl      input        powerlaw index
*tk      input        temperature Kelvin
*ist     input        ionization state
*zmax    input        redshift of source
*h0      input        cosmological H0 km s-1 Mpc-1
*omegam  input        cosmological omegam
*omegal  input        cosmological omegal
*ne      input        number of energy bins
*ear     input        array of energy bin boundaries keV in observer frame
*xsec    in/out       X-ray cross section array
*tau     output       optical depth for each energy
*-Dick Willingale 2011-May-24
        real cmpermpc,c
        parameter (cmpermpc=3.08568e24,c=2.99792458e5)
        real param(6),exsec
        INTEGER nz,j,i
        real zsam,zz,z1,zf
C initialise optical depth
        do j=1,ne
                tau(j)=0.0
        enddo
C integrate over redshift using delz approx 0.02
        nz=nint(zmax/0.02)
        zsam=zmax/nz
        zz=zsam*0.5
        do i=1,nz
                z1=zz+1.0
                zf=(z1**2/sqrt(omegam*z1**3+omegal))
                zf=zf*zsam*c*n0*cmpermpc/h0
C Call Xpsec routine to calculate cross-section vs. energy
                param(1)=pl
                param(2)=1.e-4
                param(3)=tk
                param(4)=ist
                param(5)=zz
                param(6)=1.0
                call xsabsori_tau(ear,ne,param,1,xsec,exsec)
C sum optical depth
                do j=1,ne
                        xsec(j)=xsec(j)*1.e-22/param(2)
                        tau(j)=tau(j)+xsec(j)*zf
                enddo
                zz=zz+zsam
        ENDDO
        END
