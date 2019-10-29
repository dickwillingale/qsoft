*+QRA_GRBDUST	include dust scattering flux component in afterglow
	SUBROUTINE QRA_GRBDUST(tdata,np,pr,ipul,elo,ehi,nc,rates,flux)
	integer np,nc,ipul
	double precision tdata,pr(np),elo,ehi,rates(nc),flux
*tdata	input	data times (wrt BAT trigger)
*np	input	number of parameters
*pr	input	parameters
*ipul	input	0 afterglow, >0 individual pulse, -1 all, -2 all prompt
*elo	input	low energy of flux band keV
*ehi	input	high energy of flux band keV
*nc	input	number of energy channels (8, 1-4 XRT, 5-8 BAT)
*rates	output	predicted count rates
*flux	output	predicted flux keV cm-2 s-1 in band elo-ehi keV
* The parameters are:
*1       Ns      number of pulses (fixed in fitting)
*2       b1      low spectral index (with sign)
*3       amb	 b1-b2 change in spectral index (>0.0) (with sign)
*4       d       temporal index of nu_p  (-1 in initial description)
*5       a       temporal index of L(nu_p) (1 in initial description)
*6	 fr      rise time expressed as fraction of Tf
*		 time of 1st photon Tz=Tf*(1-fr)
*7       Tf      time of last photon
*8       Epk     the break or characteristic energy at Tpk
*9       Tpk     peak of pulse wrt trigger (data) time
*10    	 Spk     Flux kev cm-2 s-1 0.3 to 350 keV at time Tpk
* Repeat of pulse parameters if Ns>1
* Possible header parameter for afterglow
*11	 na	 type of afterglow model
*		 100	with evolution of spectral index
*	 vv	 spectral index evolution index
*		 200	a dust scattering model
*	 nd	 dust column density, cm-2
*	 dd      dust size distribution index
*     	 rl      minimum distance to dust, pc
*     	 rh      maximum distance to dust, pc
* A normal afterglow component can be included when na=200
* In original version na was missing or na=100
*11      ba	 spectral index of afterglow at Ta
*12	 ga      Break energy (if 0 then 1000 keV)
*13	 Tr	 rise time of after glow
*14	 Ta	 transition time from exponential decay to powerlaw in afterglow
*15	 Fm	 maximum flux kev cm-2 s-1 0.3 to 10.0 keV 
*16	 aa	 temporal index of final afterglow powerlaw decay
*17	 Tb	 time of late break
*18	 ab	 increment in powerlaw decay at Tb
*19	 Tb1	 time of further break
*20	 ab1	 decay index after this break
*21	 Tb2	 time of yet another break
*22	 ab2	 decay index after this break
* The injection time Tin=Tpk-Tf
* Tf timescale associated with pulse measured wrt Tin as above
* Tr and Ta are timescales for afterglow measured wrt Tpk of 1st pulse
* In this qR version we assume the wide fitted band of 0.3-350 keV
*-Author Dick Willingale 2013-Mar-14
	integer i,istat,k,kk,Ns,na
	double precision alpha,amb,ec,obsi,Epk,Spk,vv
	double precision b1,b2,d,a,Tz,Tf,Tin,Fpk,Apk,ba,ga,Tr
	double precision Ta,Fm,aa,lfa,Fa,Tpk,fr,Tb,ab,a2,Tb1,ab1,Tb2,ab2,dif
	double precision rata(10),t,rt
C Initialise rates array
	do i=1,nc
		rates(i)=0.0
	enddo
	flux=0.0
C Get number of prompt pulses
	kk=1
	Ns=pr(kk)
C Loop for prompt pulses
	do k=1,Ns
	    kk=kk+1
	    b1=pr(kk)
	    kk=kk+1
	    amb=pr(kk)
	    kk=kk+2
	    a=pr(kk)
	    d=-a
	    kk=kk+1
	    fr=pr(kk)
	    kk=kk+1
	    Tf=pr(kk)
	    kk=kk+1
	    Epk=pr(kk)
	    kk=kk+1
	    Tpk=pr(kk)
	    kk=kk+1
	    Spk=pr(kk)
C Trap Tf too small
      	    Tf=max(Tf,1.e-2)
C Trap rise fraction out of usable range
      	    fr=max(min(fr,0.99),0.01)
C Trap Spk < 0
      	    Spk=max(Spk,0.0)
C Set Tzero (time of first photon)
	    Tz=Tf*(1.0-fr)
C Set time wrt injection time
	    Tin=Tpk-Tf
	    t=tdata-Tin
	    if( (ipul.lt.0.or.ipul.eq.k) .and. (t.gt.Tz) ) then
C set low photon index for prompt
		alpha=b1-1.0
C Integrate band function to convert flux to normalisation at 1 keV
		istat=0
C Note here we use the wide energy band 0.3 - 350.0 keV
		call qra_bandint(0.3D0,350.0D0,alpha,alpha-amb,Epk,obsi)
		Fpk=Spk/obsi
C time evolution of the characteristic energy of Band function of prompt
		ec=Epk*(t/Tf)**d
C Apk is normalisation at 1 keV at time t
		Apk=Fpk*(Epk/ec)**alpha
C Get spectrum rates
		CALL qra_specrates(alpha,amb,ec,nc,rata,istat)
C Calculate flux
		call qra_bandint(elo,ehi,alpha,alpha-amb,ec,obsi)
		if(istat.eq.0) then
C Set overall normalisation depending on time
			rt=((t/Tf)**(-1.0))*Apk
			a2=a+2.0
			rt=rt*(min(t,Tf)**a2-Tz**a2)
			rt=rt/(Tf**a2-Tz**a2)
C Multiply rates by normalisation
			do i=1,nc
				rates(i)=rates(i)+rata(i)*rt
			enddo
C Add in flux
			flux=flux+obsi*rt
		endif
	    endif
	enddo
C return if no afterglow parameters set or afterglow not required
	if(np.eq.kk.or.ipul.gt.0.or.ipul.eq.-2) return
C Select afterglow type if set
	if(pr(kk+1).ge.100) then
			na=pr(kk+1)
			if(na.eq.100) then
C Afterglow includes evolution of spectral index
				kk=kk+2
				vv=pr(kk)
			elseif(na.eq.200) then
C Dust scattering included
				kk=kk+5
				vv=0.0
			endif
	else
		vv=0.0
	endif
C Add in afterglow with reference time peak of 1st prompt pulse
	t=tdata-pr(9)
	ba=pr(kk+1)
	ga=pr(kk+2)
	Tr=pr(kk+3)
	Ta=pr(kk+4)
	Fm=pr(kk+5)
	aa=pr(kk+6)
	Tb=0.0
	ab=0.0
	Tb1=0.0
	ab1=0.0
	Tb2=0.0
	ab2=0.0
	if(np.gt.kk+6) then
		Tb=pr(kk+7)
		ab=pr(kk+8)
		if(np.gt.kk+8) then
			Tb1=pr(kk+9)
			ab1=pr(kk+10)
			if(np.gt.kk+10) then
				tb2=pr(kk+11)
				ab2=pr(kk+12)
			endif
		endif
	endif
	if(t.gt.0.0) then
C set low photon index at time Ta
		alpha=ba-1.0
		amb=10.0
		if(ga.eq.0) then
			ec=1000.
		else
			ec=ga
		endif
C Evolve cut-off energy with time
C		ec=ec*(t/Ta)**d
C Evolve spectral index
		if(t.gt.ta) then
			if(t.gt.tb.and.ab.gt.0.0) then
				alpha=alpha*(tb/Ta)**vv
			else
				alpha=alpha*(t/Ta)**vv
			endif
		endif
C Integrate band function to convert flux to normalisation at 1 keV
		istat=0
C For afterglow use the XRT energy band only
		call qra_bandint(0.3D0,10.0D0,alpha,alpha-amb,ec,obsi)
		Fm=Fm/obsi
C Get spectrum rates for afterglow
		call qra_specrates(alpha,amb,ec,nc,rata,istat)
C Calculate flux
		call qra_bandint(elo,ehi,alpha,alpha-amb,ec,obsi)
C Find Fa from Fm
		if(Tr.le.aa*Ta) then
			Fa=Fm*exp(-aa+2.0*sqrt(Tr*aa/Ta))
		else
			Fa=Fm*exp(aa)*(Tr/(aa*Ta))**aa
		endif
C Set normalisation for afterglow
		if(t.lt.Ta) then
			lfa=log10(Fa)+(aa-Tr/t-t*aa/Ta)*log10(exp(1.0))
		else
			lfa=log10(Fa)-aa*log10(t/Ta)-(Tr/t)*log10(exp(1.0))
		endif
		if(ab.ne.0.0.and.t.gt.Tb) then
			lfa=lfa-ab*log10(t/Tb)
		endif
		if(ab1.ne.0.0.and.t.gt.Tb1) then
			dif=ab1-(aa+ab)
			lfa=lfa-dif*log10(t/Tb1)
		endif
		if(ab2.ne.0.0.and.t.gt.Tb2) then
			dif=ab2-ab1
			lfa=lfa-dif*log10(t/Tb2)
		endif
		if(istat.eq.0) then
			rt=10.0**(max(lfa,-20.))
C add in afterglow
			do i=1,nc
				rates(i)=rates(i)+rata(i)*rt
		
			enddo
C add flux
			flux=flux+obsi*rt
		endif
	endif
	END
