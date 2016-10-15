#!/usr/bin/env Rscript
# Calculate reflectivity of Ir with B4C overcoat
	emin<- 0.1
	emax<- 15.0
	esam<- 0.05
	#lemin<- log10(emin)
	#lemax<- log10(emax)
	#lesam<- (lemax-lemin)/(ne-1)
	#lekev<- ((1:ne)-1)*lesam+lemin
	#ekev<- 10^lekev
	ekev<- seq(from=emin,to=emax,by=esam)
	ne<- length(ekev)
#
	irdensity<- 22.65
	fir<- qrt_xopt("Ir",irdensity,ekev,1)
	nri<- -fir$alpha/2.0+1.0
	nii<- fir$gamma/2.0
#
	b4cdensity<- 2.52
	fbc<- qrt_xopt("B4 C",b4cdensity,ekev,1)
	nrc<- -fbc$alpha/2.0+1.0
	nic<- fbc$gamma/2.0
#
	nrlay<- double(length=3)
	nilay<- double(length=3)
	dlay<- double(length=3)
	nrlay[1]<- 1
	nilay[1]<- 0.0
	dlay[1]<- 0.0
	dlay[2]<- 80.0
	dlay[3]<- 0.0
#
	na<- 75
	amin<- 87.0
	amax<- 89.9
	asam<- (amax-amin)/(na-1)
	iangle<- ((1:na)-1)*asam+amin
	gangle<- -iangle+90.0
#
	refs<- double(length=na*ne)
	dim(refs)=c(na,ne)
#
	for(k in 1:ne) {
		nrlay[2]=nrc[k]
		nilay[2]=nic[k]
		nrlay[3]=nri[k]
		nilay[3]=nii[k]
		mdat<- qrt_mlayer(iangle,ekev[k],nrlay,nilay,dlay,1)
		refs[1:na,k:k]<- (mdat$rsig+mdat$rpi)*0.5
	}
#
	X11()
#
	plot(ekev,refs[1:1,1:ne],type="l",log="x",xlab="keV",ylab="Ref")
	for(k in 2:na) {
		lines(ekev,refs[k:k,1:ne],type="l")
	}
#
	go_on<- locator(1)
# Save result
	ir_b4cover<-c()
	ir_b4cover$ekev<- ekev
	ir_b4cover$angs<- iangle
	ir_b4cover$refs<- refs
	save(ir_b4cover,file="ir_b4cover.RData")
