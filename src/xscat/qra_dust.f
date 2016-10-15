*+QRA_DUSTRINGS Fitting to dust X-ray scattering halo rings
        subroutine qra_dustrings(nrings,ndata,dat,der,dsou,ne,ekv,srate,
     +  dts,td,zd,amin,amax,qa,sig1,angs,ndust,edust,model,
     +  chisq,ndof)
        implicit none
        integer nrings,ndata,ne
        double precision dat(nrings,ndata),der(nrings,ndata)
        double precision dsou,ekv(ne),srate(ne),dts,td(ndata)
        double precision zd(nrings),amin,amax,qa,sig1
        double precision angs(nrings,ndata),ndust(nrings),edust(nrings)
        double precision model(nrings,ndata),chisq
        integer ndof
Cf2py  intent(in) nrings,ndata,dat,der,dsou,ne,ekv,srate,dts
Cf2py  intent(in) td,zd,amin,amax,qa,sig1
Cf2py  intent(out) angs,ndust,edust,model,chisq,ndof
C nrings     input number of expanding rings
C ndata      input number of observation times per ring
C dat        input array of rings surface brightness data cts/s/str
C der        input array of errors on data
C dsou       input distance to source PC
C ne         input nummber of energy samples
C ekv        input array of energies keV (equally spaced across observed band)
C srate      input source spectrum cts/s/keV
C dts        input source burst duration
C td         input delay time of observations secs
C zd         input fraction of source distance to rings
C amin,amax  input grain size radius range microns
C qa         input grain size distribution index N(a)=A.a^-qa
C sig1       input differential cross-section of 1 grain cm2, 1 keV, 0.1 microns
C angs       output angs
C ndust      output N dust columns cm-2
C edust      output errors on N dust columns cm-2
C model      output model cts/s/str in rings
C chisq      output Chi-Squared between data and model
C ndof       output ndof
*-Dick Willingale 2015-Aug-05
        double precision esam,pi,ao
        double precision sigma,sumd,summ,sume,diffsec,fcos,asc,zz
        double precision cms,pctom
        parameter (cms=2.99792458d8)
        parameter (pctom=3.08568d16)
        parameter (pi=3.141592653589793)
        external diffsec
        integer i,k,j,nn
        esam=ekv(2)-ekv(1)
C Loop for rings
        do i=1,nrings
C Loop for data at delay times
            sumd=0.0
            summ=0.0
            do k=1,ndata
C get scatter angle from delay time
                call qra_dustthetascat(td(k),dsou,zd(i),asc)
C observed ring angle
                zz=1.0-zd(i)
                ao=asc*zz
                angs(i,k)=ao
C Calculate scaling factor for solid angle and integration wrt td and
                fcos=esam*dts*cos(asc-ao)*zd(i)/zz/td(k)
                model(i,k)=0.0
C Loop over energy
                do j=1,ne
C differential cross-section for delay time
                    sigma=diffsec(ekv(j),amin,amax,qa,sig1,asc)
                    model(i,k)=model(i,k)+sigma*srate(j)*fcos
                enddo
                sumd=sumd+dat(i,k)
                summ=summ+model(i,k)
            enddo
            ndust(i)=sumd/summ
            sume=0.0
            do k=1,ndata
                model(i,k)=model(i,k)*ndust(i)
                sume=sume+(dat(i,k)-model(i,k))**2
            enddo
            edust(i)=ndust(i)*sqrt(sume)/sumd
        enddo
C Calculate Chi-squared
        chisq=0.0
        nn=0
        do i=1,nrings
            do k=1,ndata
                if(der(i,k).gt.0.0) then
                    chisq=chisq+(dat(i,k)-model(i,k))**2/der(i,k)**2
                    nn=nn+1
                endif
            enddo
        enddo
        ndof=nn-(nrings+2)
        end
*+QRA_THETASCAT Scattering angle radians
        subroutine qra_dustthetascat(td,ds,zd,arad)
        implicit none
        double precision td,ds,zd,arad
Cf2py  intent(in) td,ds,zd
Cf2py  intent(out) arad
* td delay time seconds after source flare (assumed delta function)
* ds distance to source parsecs (convert to m - 3.086e16/parsec)
* zd fractional distance to dust 
*-Dick Willingale 2015-Aug-05
        double precision cms,pctom
        parameter (cms=2.99792458d8)
        parameter (pctom=3.08568d16)
        arad=sqrt(td*2.0*cms/(1.0-zd)/(zd*ds*pctom))
        end
C Internal routines
C Angular dependence of the differential scattering cross-section
        double precision function scatfun(agr,ts,lambda)
        implicit none
        double precision agr,ts,lambda
C agr                grain size microns
C ts                scattering angle radians
C lambda        wavelength microns
*-Dick Willingale 2016-May-12
        double precision x,pi,aol
        parameter (pi=3.141592653589793)
        aol=4.0*pi*sin(ts/2.0)/lambda
        x=aol*agr
C Rayleigh-Gans approximation for angular dependence
C (1-cos(theta)^2)*(j1(x)/x)**2 the Rayleigh-Gans approximation for angular dependence
C of X-ray scattering from dust grains
        scatfun=(((sin(x)/x**2-cos(x)/x)/x)**2)*(1.0+cos(ts)**2)*0.5
        end
C Differential cross-section Rayleigh-Gans approximation
        double precision function diffsec(ekv,amin,amax,qa,sig1,ts)
        implicit none
        double precision ekv,amin,amax,qa,sig1,ts
C differential scattering cross section of dust grains
C Based on Mauche and Gorenstein 1986 and Smith and Dwek 1998
C ekv photon energy keV
C amin, amax grain size range microns
C qa grain size distribution index
C sig1 diff cross section at 1 keV for grain size 0.1 microns
C ts scattering angle radians
        integer na,k
        double precision dela,q1,ain,anorm,pi,sig,agr
        double precision scatfun,micpkev,lambda
        parameter (micpkev=12.398425d-4)
        parameter (pi=3.141592653589793)
        external scatfun
C Loop to integrate over grain size distribution
        na=20
        dela=(amax-amin)/na
C Normalisation for the integral of N(a)=A.a**qa
        q1=1.0-qa
        ain=dela*q1/(amax**q1-amin**q1)
C Peak value of (j1(x)/x)**2 distribution is 1/9
C Factor 1.0e6 because sig1 for a=0.1 microns
        anorm=9.0e6*sig1*ain
C Note convert keV to wavelength in microns lam=12.4e-4/ekv
        lambda=micpkev/ekv
        sig=0.0
        do k=0,na-1
            agr=amin+dela*(k+0.5)
            sig=sig+scatfun(agr,ts,lambda)*agr**(6.0-qa)
        enddo
        diffsec=sig*anorm
        end
