*+XSABSORI_TAU fortran mutation of XSPEC C++ routine
      SUBROUTINE xsabsori_tau(ear, ne, param, ifl, tau, etau)
      implicit none
      INTEGER ne, ifl
      REAL ear(0:ne), param(6), tau(ne), etau(ne)
*- Dick Willingale 2014-Apr-10
c Arguments:
c     ear       r        i: the energy ranges on which to calculate the model
c     ne        i        i: the number of energy ranges
c     param     r        i: see below
c     ifl       i        i: the dataset
c     tau       r        r: opacity
C Parameters:
C       1       gamma - photon index of continuum power spectrum
C       2       nh - Hydrogen column of the absorper in units of 10^22
C       3       temp  - absorber temperature in K
C       4       xi    - absorber ionization state = L/nR^2
C       5       z     - redshift
C       6       ab0   - iron abundance relative to the solar iron abundance
C algorithm:
C               a(x) = a(x)*exp(-Nh*sigma(x))
	include 'SPX_COM'
	real Gamma, nH, Temp, Xi, zshift, FeAbund
	integer nq,istat,id,i,j,k,iatom,nion,nionsave
	character qsoft*100,dfile*100
	real amass,alpha,lemin,lemax,dele,denom,alphai
	real xil,t4,tfact,mult,ratsum,intgral
	real arec,z2,y,ssum,enev
	real emin,emax,re,rp,rp1,rp2
	integer ie
C Copy parameters into named varables
	Gamma=param(1)
	nH=param(2)
	Temp=param(3)
	Xi=param(4)
	zshift=1.0+param(5)
	FeAbund=param(6)
C Get warm absorber data
	if(iwarmset.eq.0) then
C Open ion data file
		istat=0
		CALL SYS_GETENV('QSOFT',qsoft,nq,istat)
		dfile=qsoft(1:nq)//'/data/iond_mm.dat'
		CALL SYS_GETLUN(id,istat)
		OPEN(UNIT=id,FILE=dfile,STATUS='OLD')
C Read ion data
		nion=0
		do i=1,maxatom
			read(id,*) AtomicNumber(i),amass
			do j=1,AtomicNumber(i)
				nion=nion+1
				read(id,*) (ion(nion,k),k=1,maxcoef)
C change units of coefficients
		      		ion(nion,2) = ion(nion,2)*1.0E+10;
				ion(nion,4) = ion(nion,4)*1.0E+04;
				ion(nion,5) = ion(nion,5)*1.0E-04;
				ion(nion,7) = ion(nion,7)*1.0E-04;
			enddo
		enddo
		close(id)
C Open photoionization data from Reilman and Manson ApJ supp 40 1979
		dfile=qsoft(1:nq)//'/../data/mansig.dat'
		CALL SYS_GETLUN(id,istat)
		OPEN(UNIT=id,FILE=dfile,STATUS='OLD')
C Read photoionization data
		nion=0
		do i=1,maxatom
		 	read(id,*) iatom 
			do j=1,AtomicNumber(i)
C now go over all ion states - values come from Reilman
C ApJ supp 1979 except hydrogenic which is analytic
				nion=nion+1
				if(j.gt.1) then
				    read(id,*) iatom 
				endif
				do k=1,maxenerg
					read(id,*) Energy(k),sigma(nion,k)
					sigma(nion,k)=sigma(nion,k)/6.6e-27
				enddo
			enddo
		enddo
		close(id)
		iwarmset=1
	endif
C Calculate power law photon spectrum on cross-section energy grid
C spec() is N(E)*delE
	alpha=2.0-Gamma
	emin=Energy(1)
	emax=Energy(maxenerg)
	if(abs(alpha).lt.1.e-10) then
		denom=log(emax)-log(emin)
	else
		denom=(emax**alpha-emin**alpha)/alpha
	endif
	do k=1,maxenerg
		if(k.eq.1) then
			dele=(Energy(2)-Energy(1))/2.0
		elseif(k.eq.maxenerg) then
			dele=(Energy(maxenerg)-Energy(maxenerg-1))/2.0
		else
			dele=(Energy(k+1)-Energy(k-1))/2.0
		endif
		spec(k)=(dele/denom)*Energy(k)**(-Gamma)
	enddo
C Set handy variables
	if(Xi.le.0.0) then
		xil=-100.0
	else
		xil=log(Xi)
	endif
	t4=Temp*.0001
	tfact=1.033e-3/sqrt(t4)
C Loop over elements	
	nion=0
	do i=1,maxatom
		mult=0.0
		ratsum=0.0
		nionsave=nion
		do j=1,AtomicNumber(i)
			nion=nion+1
C Do integral on cross-section energy grid
		 	intgral=0.0
		 	do k=1,maxenerg
		 		intgral=intgral+sigma(nion,k)*spec(k)
		 	enddo
C calculate recombination coefficients - radiative and dielectric
C from Aldrovandi+Pequignot 1973 Astro, astrophys 25 137
		       if (j.lt.AtomicNumber(i)-1) then
  		             arec=ion(nion,2)*t4**(-ion(nion,3))+
     +			     ion(nion,4)*t4**(-1.5)*exp(-ion(nion,5)/t4)*
     +			     (1.0+ion(nion,6)*exp(-ion(nion,7)/t4))
			else
C  do radiative recomb to hydrogenic ion NB assumed gaunt factor=1
C  from Gould+Thakur 1970 ann phys 61 351.
			     z2=AtomicNumber(i)*AtomicNumber(i)
			     y=15.8*z2/t4
		             arec=tfact*z2*(1.735+log(y)+1.0/(6.*y))
			endif
			ratio(j)=log(3.2749e-6*intgral/arec)
			ratsum=ratsum+ratio(j)
C calculate population numbers
C  set sum=1 as the sum=1+A21+A32A21+....
C  complicated because floating range limit
		        mul(j) = ratsum + j*xil
		        if ( mul(j) .gt. mult ) then
			 	mult = mul(j)
			endif
		enddo
		ssum=0.0
		do j=1,AtomicNumber(i)
			mul(j)=mul(j)-mult
			ssum=ssum+exp(mul(j))
		enddo
		ssum=ssum+exp(-mult)
		nion=nionsave+1
		num(nion)=-mult-log(ssum)
		do j=2,AtomicNumber(i)
			nion=nion+1
			num(nion)=num(nion-1)+ratio(j-1)+xil
		enddo
		nion=nionsave
		do j=1,AtomicNumber(i)
			nion=nion+1
			num(nion)=exp(num(nion))
		enddo
	enddo
C Set abundances
	do i=1,maxatom
		call fgabn(ElementNames(i),AtomAbunds(i))
		if(AtomicNumber(i).eq.26) then
			AtomAbunds(i)=AtomAbunds(i)*FeAbund
		endif
	enddo
	k=1
	do ie=1,ne
		tau(ie)=0.0
C Redshift energy of bin centre in eV
		enev=(ear(ie)+ear(ie-1))*1.e3*zshift/2.0
		if(enev.ge.emax) then
			re=(enev/emax)**(-3.0)
			nion=0
			do i=1,maxatom
				rp=AtomAbunds(i)*re
				do j=1,AtomicNumber(i)
					nion=nion+1
					tau(ie)=tau(ie)+num(nion)*
     +					sigma(nion,721)*rp
				enddo
			enddo
		elseif(enev.lt.emin) then
			tau(ie)=1.e20
		else
C pick closest energy bin
			do while(Energy(k).lt.enev)
				k=k+1
			enddo
			re=(enev-Energy(k-1))/(Energy(k)-Energy(k-1))
			nion=0
			do i=1,maxatom
				rp1=AtomAbunds(i)*re
				rp2=AtomAbunds(i)*(1.0-re)
				do j=1,AtomicNumber(i)
					nion=nion+1
					tau(ie)=tau(ie)+num(nion)*
     + 					(rp1*sigma(nion,k)+
     +					rp2*sigma(nion,k-1))
				enddo
			enddo
		endif
	enddo
	do ie=1,ne
		tau(ie)=nH*tau(ie)*6.6e-5
	enddo
	end
