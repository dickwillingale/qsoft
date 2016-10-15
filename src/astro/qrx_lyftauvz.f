*+QRX_LYFTAUVZ    X-ray optical depth of Lyman Forest vs. z
        subroutine qrx_lyftauvz(n0,dind,m0,mind,h0,omegam,omegal,
     +  ekev,nz,z,tau)
        implicit none
        INTEGER nz
        REAL n0,dind,m0,mind,h0,omegam,omegal,ekev
        REAL z(nz),tau(nz)
Cf2py  intent(in) n0,dind,m0,mind,h0,omegam,omegal,ekev,nz,z
Cf2py  intent(out) tau
*n0        input        number density cm-3 at z=0
*dind      input        density index wrt z+1
*m0        input        metalicity log[X/H] at z=0
*mind      input        metalicity log[X/H] index wrt z
*h0        input        cosmological H0 km s-1 Mpc-1
*omegam    input        cosmological omegam
*omegal    input        cosmological omegal
*ekev      input        photon energy keV
*nz        input        number of redshift bins
*z         input        array of redshifts
*tau       output       optical depth for each redshift
*-Dick Willingale 2015-May-10
        include 'SPX_COM'
        real cmpermpc,c
        parameter (cmpermpc=3.08568e24,c=2.99792458e5)
        real param(2),exsec,xsec
        INTEGER ne,i
        real dz,zz,z1,zf,smet,ear(0:1)
C initialise energy bin
        ear(0)=ekev-0.001
        ear(1)=ekev+0.001
        ne=1
C Save current metalicity
        smet=metalicity
C set up parameters
        param(1)=1.e-4
        param(2)=1.0
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
C call Xpsec routine to calculate cross-section vs. energy
                param(2)=zz
                call xszphb(ear,ne,param,1,xsec,exsec)
C convert absorption factor to cross sections
                xsec=-log(xsec)*1.e-18
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
