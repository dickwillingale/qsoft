*+QRA_COSMO Cosmology calculations
        subroutine qra_cosmo(h0,omegam,omegal,dz,maxnum,z,ez,dc,dm,
     +         da,dl,distm,dvc,tlbak,vc,thsec,thgyr,dhmpc)
        implicit none
        integer maxnum
        double precision h0,omegam,omegal,dz,thsec,thgyr,dhmpc
        double precision z(maxnum),ez(maxnum)
        double precision dc(maxnum),dm(maxnum),da(maxnum),dl(maxnum)
        double precision distm(maxnum),dvc(maxnum),tlbak(maxnum)
        double precision vc(maxnum)
Cf2py  intent(in) h0,omegam,omegal,dz
Cf2py  intent(in) maxnum
Cf2py  depend(maxnum) z,ez,dc,dm,da,dl,distm,dvc,tlbak,vc,thsec,thgyr,dhmpc
Cf2py  intent(out) z,ez,dc,dm,da,dl,distm,dvc,tlbak,vc,thsec,thgyr,dhmpc
*h0      input   Hubble constant km s-1 Mpc-1
*omegam  input   Omega_m
*omegal  input   Omega_lamdba
*dz      input   redshift increment for quadrature (integrals start at z=0)
*maxnum  input   number of integration steps - number of redshift values
*z       output  redshift values
*ez      output  scaling function
*dc      output  comoving line-of-sight distance
*dm      output  transverse comoving distance
*da      output  angular diameter distance
*dl      output  luminosity distance
*distm   output  distance modulus
*dvc     output  comoving volume element
*tlbak   output  look back time
*vc      output  integrated comoving volume over whole sky
*thsec   output  Hubble time in seconds
*thgyr   output  Hubble time in Giga years
*dhmpc   output  Hubble length in Mpc
*Equations from David Hogg astro-ph/9905116
*-Dick Willingale 2011-Feb-25
        double precision ez0,tmp1,tmp2,tmp3,omegak
        double precision escale,arcsinh
        integer i
        double precision pi,fourpi,ckms,spgy,mpctokm
        parameter (pi=3.1415926535898)
        parameter (fourpi=pi*4.d0)
        parameter (ckms=2.99792458d5)
        parameter (spgy=3.15576e16)
        parameter (mpctokm=3.08568d19)
c        write(*,*) 'start cosmo maxnum',maxnum
c
        omegak = 1.0d0 - omegam - omegal
c Hubble distance in Mpc
        dhmpc = ckms / h0
c Hubble time in sec
        thsec = mpctokm / h0
c Hubble time in Gyrs
        thgyr = thsec / spgy
c set up redshift values
        do i=1,maxnum
              z(i)=dz*i
        enddo
C        write(*,*) 'set z'
C dc = dh * integral(0,z, dz'/E(z'))
        ez0 = escale(0.0d0,omegam,omegal,omegak)
        ez(1) = escale(z(1),omegam,omegal,omegak)
        dc(1) = dhmpc * dz * (1./ez(1) + 1./ez0) / 2.
        tlbak(1) = thgyr * dz * ( 1./(1.+z(1))/ez(1) + 1./ez0) / 2.
        do i=2,maxnum
                ez(i) = escale(z(i),omegam,omegal,omegak)
                dc(i) = dc(i-1) + dhmpc * dz * (1./ez(i)+1./ez(i-1))/2.
                tlbak(i) = tlbak(i-1) + thgyr * dz *
     +                (1./(1.+z(i))/ez(i) + 1./(1.+z(i-1))/ez(i-1))/2.
        enddo
C        write(*,*) 'set ez'
        do i=1,maxnum
C                write(*,*) i,dc(i),dm(i)
c  transverse comoving distance
                if (abs(omegak) .lt. 1.e-4) then
                        dm(i) = dc(i)
                elseif (omegak .gt. 0.) then
                        dm(i) = dhmpc / sqrt(omegak) * 
     +                        sinh(sqrt(omegak)*dc(i)/dhmpc)
                else
                        dm(i) = dhmpc / sqrt(abs(omegak)) * 
     +                        sin(sqrt(abs(omegak))*dc(i)/dhmpc)
                endif
C                write(*,*) 'set dm'
c  angular diameter distance
                da(i) = dm(i) / (1.0d0+z(i))
C                write(*,*) 'set da'
c  luminosity distance
                dl(i) = dm(i) * (1.0d0+z(i))
C                write(*,*) 'set dl'
c  distance modulus.  dl(i) is in Mpc and magnitudes are referred to 10 pc.
                distm(i) = 5.0d0 * log10(dl(i)*1.0d5)
C                write(*,*) 'set distm'
c  comoving volume element, one needs to multiply this by
c  dOmega dz (i.e. solid angle * dz) to get dV_c
                dvc(i) = dhmpc*(1.d0+z(i))**2*da(i)**2/ez(i)
C                write(*,*) 'set dvc'
c actual integrated comoving volume, out to redshift z, over
c the whole sky (so multiply by dOmega/4pi for a given solid angle)
                if (abs(omegak) .lt. 1.d-4) then
                      vc(i) = fourpi/3.0d0 * dm(i)**3
                else 
                      tmp1=fourpi*dhmpc**3 / 2.0d0 / omegak
                      tmp2=dm(i)/dhmpc*sqrt(1.+omegak*(dm(i)/dhmpc)**2)
                      tmp3=sqrt(abs(omegak)) * dm(i) / dhmpc
                      if (omegak .gt. 0.0d0) then
                              vc(i) = tmp1*(tmp2 - 1./sqrt(abs(omegak))
     +                              *arcsinh(tmp3))
                      else
                              vc(i) = tmp1*(tmp2 - 1./sqrt(abs(omegak))
     +                              *asin(tmp3))
                      endif
                endif
C                write(*,*) 'set vc'
        enddo
        end
*+ESCALE scaling function for cosmological distance calculation
        double precision function escale(z,omegam,omegal,omegak)
        implicit none
        double precision z,omegam,omegal,omegak
*-
        escale = omegam*(1.+z)**3 + omegal + omegak*(1.+z)**2
        escale=sqrt(escale)
        end
*+ARCSINH        for cosmological distance integrals
        function arcsinh(x)
        implicit none
        double precision arcsinh,x
        arcsinh = sign(1.d0,x)*log(abs(x) + sqrt(1.0d0 + x**2))
        end
