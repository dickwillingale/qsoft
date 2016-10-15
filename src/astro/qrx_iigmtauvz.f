*+QRX_IIGMTAUVZ    X-ray optical depth of ionized IGM gas vs. redshift
        subroutine qrx_iigmtauvz(n0,dind,m0,mind,pl,tk,ist,
     +  h0,omegam,omegal,ekev,nz,z,tau)
        implicit none
        INTEGER nz
        REAL n0,dind,m0,mind,pl,tk,ist
        REAL h0,omegam,omegal,ekev,z(nz),tau(nz)
Cf2py  intent(in) n0,dind,m0,mind,pl,tk,ist,h0,omegam,omegal,ekev,nz,z
Cf2py  intent(out) tau
*n0     input        number density cm-3 at z=0
*dind   input        density index wrt z+1
*m0     input        metalicity log[X/H] at z=0
*mind   input        metalicity log[X/H] index wrt z
*pl     input        powerlaw index
*tk     input        temperature Kelvin
*ist    input        ionization state
*h0     input        cosmological H0 km s-1 Mpc-1
*omegam input        cosmological omegam
*omegal input        cosmological omegal
*ekev   input        energy keV
*nz     input        number of redshift values
*z      input        array of redshift values
*tau    output       optical depth for each energy
*-Dick Willingale 2015-May-09
        include 'SPX_COM'
        real cmpermpc,c
        parameter (cmpermpc=3.08568e24,c=2.99792458e5)
        real param(6),exsec,xsec,smet
        INTEGER ne,i
        real dz,zz,z1,zf,ear(0:1)
C initialise energy bin
        ear(0)=ekev-0.001
        ear(1)=ekev+0.001
        ne=1
C Save current metalicity
        smet=metalicity
C set up parameters
        param(1)=pl
        param(2)=1.e-4
        param(3)=tk
        param(4)=ist
        param(5)=0.0
        param(6)=1.0
C integrate over redshift
        do i=1,nz
                if(i.eq.1) then
                        zz=z(1)*0.5
                        dz=z(1)
                else
                        zz=(z(i)+z(i-1))*0.5
                        dz=z(i)-z(i-1)
                endif
                z1=zz+1.0
                zf=(z1**2/sqrt(omegam*z1**3+omegal))
                zf=zf*dz*c*n0*cmpermpc/h0
C density evolution
                zf=zf*(z1**dind)
C metalicity evolution
                metalicity=10.0**(m0+zz*mind)
C Call Xpsec routine to calculate cross-section vs. energy
                param(5)=zz
                call xsabsori_tau(ear,ne,param,1,xsec,exsec)
C convert absorption factor to cross sections
                xsec=xsec*1.e-22/param(2)
C Sum optical depth
                if(i.eq.1) then
                        tau(i)=xsec*zf
                else
                        tau(i)=tau(i-1)+xsec*zf
                endif
      ENDDO
C reset metalicity
      metalicity=smet
      END
