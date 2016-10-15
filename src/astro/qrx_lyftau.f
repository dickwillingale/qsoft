*+QRX_LYFTAU        X-ray optical depth of Lyman Forest
        subroutine qrx_lyftau(n0,dind,m0,mind,zmax,h0,omegam,omegal,
     +  ne,ear,xsec,tau)
        implicit none
        INTEGER ne
        REAL n0,dind,m0,mind,zmax,h0,omegam,omegal
        REAL ear(0:ne),xsec(ne),tau(ne)
Cf2py  intent(in) n0,dind,m0,mind,zmax,h0,omegam,omegal,ne,ear
Cf2py  intent(out) xsec,tau
*n0        input        number density cm-3 at z=0
*dind      input        density index wrt z+1
*M0        input        metalicity log[X/H] at z=0
*mind      input        metalicity log[X/H] index wrt z
*zmax      input        redshift of source
*h0        input        cosmological H0 km s-1 Mpc-1
*omegam    input        cosmological omegam
*omegal    input        cosmological omegal
*ne        input        number of energy bins
*ear       input        array of energy bin boundaries keV in observer frame
*xsec      in/out       X-ray cross section array
*tau       output       optical depth for each energy
*-Dick Willingale 2013-Feb-06
        include 'SPX_COM'
        real cmpermpc,c
        parameter (cmpermpc=3.08568e24,c=2.99792458e5)
        real param(2),exsec
        INTEGER nz,j,i
        real zsam,zz,z1,zf,smet
C Save current metalicity
        smet=metalicity
C initialise optical depth
        do j=1,ne
                tau(j)=0.0
        enddo
C set up parameters
        param(1)=1.e-4
        param(2)=1.0
C integrate over redshift using delz approx 0.02
        nz=nint(zmax/0.02)
        zsam=zmax/nz
        zz=zsam*0.5
        do i=1,nz
                z1=zz+1.0
                zf=(z1**2/sqrt(omegam*z1**3+omegal))
                zf=zf*zsam*c*n0*cmpermpc/h0
C density evolution
                zf=zf*(z1**dind)
C metalicity evolution
                metalicity=10.0**(m0+zz*mind)
C call Xpsec routine to calculate cross-section vs. energy
                param(2)=zz
                call xszphb(ear,ne,param,1,xsec,exsec)
C convert absorption factor to cross sections and sum optical depth
                do j=1,ne
                        xsec(j)=-log(xsec(j))*1.e-18
                        tau(j)=tau(j)+xsec(j)*zf
                enddo
                zz=zz+zsam
       ENDDO
C reset metalicity
       metalicity=smet
       END
