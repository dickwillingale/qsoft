# QR routines for astrophysics
if(!exists("qrx_init")) {
qrx_init<-function() {
	.Fortran("qrx_init")
invisible()
}
qrx_init()
# Utilities
qr_quaderr<-function(x0,y0,x1,y1,y) {
# Quadratic estimator for confidence limit
# x0,y0         parameter and statistic at minimum
# x1,y1         parameter and statistic near minimum
# y             required statistic value
# returns       estimate of parameter corresponding to y
        a<- (y1-y0)/(x1-x0)^2
        b<- -2*a*x0
        c<- y0-a*x0^2-b*x0
        ba<- b^2-4*a*(c-y)
        if((ba>=0)&(a!=0)) {
                return(sqrt(ba)/(2*a))
        } else {
                return(0)
        }
}
qr_srchmin<-function(pars,pl,ph,stat,delstat,derr) {
# Search for minimum statistic and return best fit parameters and
# confidence limits of the parameters
# pars		initial parameter values
# pl		hard lower limit of parameter values
# ph		hard upper limit of parameter values
# stat		the statistic function to be minimised
# delstat	the change in statistic for confidence limits
# derr		logical array to specify parameters for confidence estimates
# returns	list from optim() plus confidence limits parlo and parhi
	np<-length(pars)
# make local copies of hard limits
	prl<-pl
	prh<-ph
# find statistic minimum and estimate hessian matrix
	ft<-optim(pars,stat,gr=NULL,method="L-BFGS-B",
	lower=prl,upper=prh,hessian=T)
# estimate errors on parameters
	ft$delstat<-delstat
	dig<- abs(diag(ft$hessian))
	delpar<- double(length=np)
# fix hard limits of parameters for which we don't want error estimate
	for(k in 1:np) {
# Trap zero on hessian diagonal
		if(dig[k]==0) {
			delpar[k]<- 0
			if(!derr[k]) {
				prl[k]<-ft$par[k]-ft$par[k]/200
				prh[k]<-ft$par[k]+ft$par[k]/200
			}
		} else {
			delpar[k]<-sqrt(2*delstat/dig[k])
			if(!derr[k]) {
				prl[k]<-ft$par[k]-delpar[k]/200
				prh[k]<-ft$par[k]+delpar[k]/200
			}
		}
	}
	cat("Min statistic",ft$value,"\n")
	cat("best fit parameters",ft$par,"\n")
	cat("diag",dig,"\n")
	cat("delpar",delpar,"\n")
	nfr<- sum(derr)
	cat("Find limits for ",nfr," parameters\n")
	parhi<-double(length=np)
	parlo<-double(length=np)
	for(k in 1:np) {
		if(derr[k]&&delpar[k]>0) {
			cat("searching for error range",k,ft$par[k],delpar[k],
			prl[k],prh[k],"\n")
# upper limit
			epar<- 0
			while(epar==0) {
				pval<-ft$par[k]+delpar[k]*nfr
				part<-ft$par
				partl<-prl
				parth<-prh
				part[k]<-pval
				partl[k]<-pval-delpar[k]/200
				parth[k]<-pval+delpar[k]/200
				ftp<-optim(part,stat,gr=NULL,method="L-BFGS-B",
				lower=partl,upper=parth)
				if(ftp$value<ft$value) {
					ft$par<- ftp$par
					ft$value<- ftp$value
					cat("new min found",ft$value,"\n")
					cat("new fit parameters",ftp$par,"\n")
					if(ftp$par[k]>prh[k]) {
						epar<- 999
					}
				} else {
					if(ftp$par[k]>prh[k]) {
					 epar<- 999
					} else {
					 epar<-qr_quaderr(ft$par[k],ft$value,
					 pval,ftp$value,ft$value+delstat)
					}
				}
			}
			cat("upper",ft$par[k],ft$value,pval,ftp$value,
			ft$value+delstat,epar,"\n")
			parhi[k]=min(ft$par[k]+epar,prh[k])
# lower limit
			epar<- 0
			while(epar==0) {
				pval<-ft$par[k]-delpar[k]*nfr
				part<-ft$par
				partl<-prl
				parth<-prh
				part[k]<-pval
				partl[k]<-pval-delpar[k]/200
				parth[k]<-pval+delpar[k]/200
				ftp<-optim(part,stat,gr=NULL,method="L-BFGS-B",
				lower=partl,upper=parth)
				if(ftp$value<ft$value) {
					ft$par<- ftp$par
					ft$value<- ftp$value
					cat("new min found",ft$value,"\n")
					if(ftp$par[k]<prl[k]) {
						epar<- 999
					}
					cat("new fit parameters",ftp$par,
					prl[k],epar,"\n")
				} else {
					if(ftp$par[k]<prl[k]) {
					 epar<- 999
					} else {
					 epar<-qr_quaderr(ft$par[k],ft$value,
					 pval,ftp$value,ft$value+delstat)
					}
				}
			}
			cat("lower",ft$par[k],ft$value,pval,ftp$value,
			ft$value+delstat,epar,"\n")
			parlo[k]=max(ft$par[k]-epar,prl[k])
		} else {
    			parhi[k]=0.0
    			parlo[k]=0.0
    		}
	}
	ft$parlo<-parlo
	ft$parhi<-parhi
	return(ft)
}
# Astrophysical routines - prefix qra_
qra_cosmo<-function(h0,omegam,omegal,zmax) {
	dz<-0.01
	maxnum<-trunc(zmax/dz)+1
	.Fortran("qra_cosmo",
	h0=as.double(h0),
	omegam=as.double(omegam),
	omegal=as.double(omegal),
	as.double(dz),
	as.integer(maxnum),
	z=double(length=maxnum),
	ez=double(length=maxnum),
	dc=double(length=maxnum),
	dm=double(length=maxnum),
	da=double(length=maxnum),
	dl=double(length=maxnum),
	distm=double(length=maxnum),
	dvc=double(length=maxnum),
	tlbak=double(length=maxnum),
	vc=double(length=maxnum),
	thsecc=double(length=1),
	thgyr=double(length=1),
	dhmpc=double(length=1))
}
qra_habs<- function(cd,ekev) {
	ne<- length(ekev)
	a<- .Fortran("qra_habs",
	as.single(cd),
	as.integer(ne),
	as.single(ekev),
	fact=single(length=ne))
	return(a$fact)
}
qra_nhtot<-function(equinox,radeg,decdeg,ssize,disio) {
	np<-length(radeg)
	a<-.Fortran("qra_nhtot",
	equinox=as.double(equinox),
	as.integer(np),
	radeg=as.double(radeg),
	decdeg=as.double(decdeg),
	ssize=as.double(ssize),
	disio=as.double(disio),
	nha=double(length=np),
	nhw=double(length=np),
	ebva=double(length=np),
	ebvw=double(length=np),
	nhma=double(length=np),
	nhmw=double(length=np),
	nhtota=double(length=np),
	nhtotw=double(length=np))
}
qra_kcorrb<-function(e1src,e2src,e1obs,e2obs,z,gamma1,gamma2,ecobs) {
	np<-length(e1obs)
	.Fortran("qra_kcorrb",
	eslo=as.double(e1src),
	eshi=as.double(e2src),
	as.integer(np),
	eolo=as.double(e1obs),
	eohi=as.double(e2obs),
	zs=as.double(z),
	gam1=as.double(gamma1),
	gam2=as.double(gamma2),
	ecob=as.double(ecobs),
	kcorr=double(length=np),
	obint=double(length=np),
	boint=double(length=np))
}
qra_grbmsb<-function(time,arr,err,smin,nmin,nmax) {
	np<- length(time)
	na<- length(arr)
	nc<- na/np
	a<-.Fortran("qra_grbmsb",
	as.integer(nc),
	as.integer(np),
	time=as.double(time),
	arr=as.double(arr),
	err=as.double(err),
	smin=as.double(smin),
	nmin=as.integer(nmin),
	nmax=as.integer(nmax),
	double(length=np),
	double(length=np),
	integer(length=np),
	nrat=integer(length=1),
	frate=double(length=np),
	efrate=double(length=np),
	rates=double(length=na),
	erates=double(length=na),
	rtime=double(length=np),
	etime=double(length=np),
	tlo=double(length=1),
	thi=double(length=1))
	dim(a$rates)<-c(nc,np)
	dim(a$erates)<-c(nc,np)
	nrat<- a$nrat
	a$rates<- a$rates[1:nc,1:nrat]
	a$erates<- a$erates[1:nc,1:nrat]
	a$rtime<- a$rtime[1:nrat]
	a$etime<- a$etime[1:nrat]
	return(a)
}
qra_grbloadpars<-function(par,grbn) {
	np<-length(par)
	ng<-nchar(grbn)
	.Fortran("qra_grbloadpars",
	as.integer(np),
	as.double(par),
	as.integer(ng),
	as.character(grbn))
}
qra_grbmodlc<-function(time,elo,ehi,ipul) {
	nt<-length(time)
	nc<-8
	nr<-nc*nt
	a<-.Fortran("qra_grbmodlc",
	as.integer(nc),
	as.integer(nt),
	time=as.double(time),
	as.double(elo),
	as.double(ehi),
	as.integer(ipul),
	rates=double(length=nr),
	flux=double(length=nt))
	dim(a$rates)<-c(nc,nt)
	return(a)
}
qra_grbloadlc<-function(tbat,rbat,ebat,txrt,rxrt,exrt,mxrt) {
	ntbat<-length(tbat)
	ntxrt<-length(txrt)
	ncbat<-length(rbat)/ntbat
	ncxrt<-length(rxrt)/ntxrt
	.Fortran("qra_grbloadlc",
		as.integer(ncbat),
		as.integer(ntbat),
		as.double(tbat),
		as.double(rbat),
		as.double(ebat),
		as.integer(ncxrt),
		as.integer(ntxrt),
		as.double(txrt),
		as.double(rxrt),
		as.double(exrt),
		as.integer(mxrt))
}
qra_grbsetfit<-function(itofit,fbat,fxrt) {
	npf<-length(itofit)
	ntbat<-length(fbat)
	ntxrt<-length(fxrt)
	.Fortran("qra_grbsetfit",
	as.integer(npf),
	as.integer(itofit),
	as.integer(ntbat),
	as.integer(fbat),
	as.integer(ntxrt),
	as.integer(fxrt),
	spars=double(length=npf),
	lpars=double(length=npf),
	upars=double(length=npf),
	bgood=integer(length=1),
	xgood=integer(length=1),
	ndof=integer(length=1))
}
qra_grbchisq<-function(fpars) {
	npf<-length(fpars)
	.Fortran("qra_grbchisq",
	as.integer(npf),
	as.double(fpars),
	chisq=double(length=1))$chisq
}
qra_grbgetfit<-function(times) {
	nr<-8
	nt<-length(times)
	nn<-nt*nr
	np<-.Fortran("qra_grbgetnpars",
	np=integer(length=1))$np
	.Fortran("qra_grbgetfit",
	as.integer(nr),
	as.integer(nt),
	as.integer(np),
	as.double(times),
	rates=double(length=nn),
	pars=double(length=np),
	ncall=integer(length=1))
}
qra_grblistmod<- function(grbfit) {
	ns<-grbfit$NS
	np<-grbfit$NP
#
	cat("** ",grbfit$GRBN,",",ns,"pulses, ",np," parameters **\n")
# Loop for all pulses
	for(i in 1:ns) {
		i1<-(i-1)*9+2
		i2<- i1+8
		cat("** Pulse",i," **\n")
   		for(k in i1:i2) {
			cat(k,grbfit$PNAME[k],grbfit$PROMF[k]," range ",
			grbfit$PROML[k],grbfit$PROMH[k],"\n")
		}   
	}
	if(np>i2) {
# Do afterglow
		i1<- ns*9+2
		i2<- length(grbfit$PROMF)
		cat("** Afterglow **\n")
		for(k in i1:i2) {
			cat(k,grbfit$PNAME[k],grbfit$PROMF[k]," range ",
			grbfit$PROML[k],grbfit$PROMH[k],"\n")
		}
	}
}
qra_grbplotsum<- function(grbfit,xrtbat,tlo,thi) {
	tblo<- grbfit$TBLO
	tbhi<- grbfit$TBHI
	txlo<- grbfit$TXLO
	txhi<- grbfit$TXHI
	txrt<- grbfit$TXRT
	etxrt<- grbfit$ETXRT
	rxrt<- grbfit$RXRT
	exrt<- grbfit$EXRT
	mxrt<- as.logical(grbfit$MXRT)
	tbat<- grbfit$TBAT
	etbat<- grbfit$ETBAT
	rbat<- grbfit$RBAT
	ebat<- grbfit$EBAT
	mbat<- as.logical(grbfit$MBAT)
	promf<- grbfit$PROMF
	grbn<- grbfit$GRBN
	pp<- qra_grbloadpars(promf,grbn)
	if(xrtbat==1) {
# plot summed BAT data
		nt<- length(tbat)
		tbat<- subset(tbat,as.logical(mbat))
		etbat<- subset(etbat,as.logical(mbat))
		rbat1<- subset(rbat[1:1,1:nt],as.logical(mbat))
		ebat1<- subset(ebat[1:1,1:nt],as.logical(mbat))
		rbat2<- subset(rbat[2:2,1:nt],as.logical(mbat))
		ebat2<- subset(ebat[2:2,1:nt],as.logical(mbat))
		rbat3<- subset(rbat[3:3,1:nt],as.logical(mbat))
		ebat3<- subset(ebat[3:3,1:nt],as.logical(mbat))
		rbat4<- subset(rbat[4:4,1:nt],as.logical(mbat))
		ebat4<- subset(ebat[4:4,1:nt],as.logical(mbat))
		xlab<- "seconds since trigger"
		ylab<- "cts/s"
		rat<- rbat1+rbat2+rbat3+rbat4
		erat<- sqrt(ebat1^2+ebat2^2+ebat3^2+ebat4^2)
		rmin<- min(rat-erat*2)
		rmax<- max(rat+erat*2)
		if(tlo==0&&thi==0) {
			xlim<-c(tblo,tbhi)
		} else {
			xlim<-c(tlo,thi)
		}
		ylim<-c(rmin,rmax)
		title<- paste("GRB",grbn," BAT  15-350 keV",sep="")
		plot(tbat,rat,type="p",cex=0,
		xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=title)
		segments(tbat-etbat,rat, tbat+etbat,rat)
		segments(tbat,rat-erat,tbat,rat+erat)
# Generate model
		nm<- 1000
	        tsm<-(tbhi-tblo)/nm
	        tbm<-((1:nm)-1)*tsm+tblo
	        pb<-qra_grbmodlc(tbm,10.0,350.0,-1)
	        mbt<-pb$time
	        mbr<-pb$rates[5:5,1:nm]
	        mbr<-mbr+pb$rates[6:6,1:nm]
	        mbr<-mbr+pb$rates[7:7,1:nm]
	        mbr<-mbr+pb$rates[8:8,1:nm]
# plot model
		lines(mbt,mbr,type="l")
	} else {
# Plot the summed XRT data
		ntx<-length(txrt)
		txrt<- subset(txrt,as.logical(mxrt))
		etxrt<- subset(etxrt,as.logical(mxrt))
		rxrt1<- subset(rxrt[1:1,1:ntx],as.logical(mxrt))
		exrt1<- subset(exrt[1:1,1:ntx],as.logical(mxrt))
		rxrt2<- subset(rxrt[2:2,1:ntx],as.logical(mxrt))
		exrt2<- subset(exrt[2:2,1:ntx],as.logical(mxrt))
		xlab<-"seconds since trigger"
		ylab<-"cts/s"
		title<- paste("GRB",grbn," XRT  0.3-10 keV",sep="")
		rat<- rxrt1+rxrt2
		erat<- sqrt(exrt1^2+exrt2^2)
		rmin<- min(rat-erat)/2
		rmax<- max(rat+erat)*2
		if(tlo==0&&thi==0) {
			xlim<- c(txlo,txhi)
			plot(txrt,rat,log="xy",type="p",cex=0,
			xlim=xlim,ylim=c(rmin,rmax),
			xlab=xlab,ylab=ylab,main=title)
		} else {
			xlim<- c(tlo,thi)
			plot(txrt,rat,log="y",type="p",cex=0,
			xlim=xlim,ylim=c(rmin,rmax),
			xlab=xlab,ylab=ylab,main=title)
		}
		segments(txrt-etxrt,rat,txrt+etxrt,rat)
		segments(txrt,rat-erat,txrt,rat+erat)
# generate model
		nm<- 3000
	        ltsm<- (log10(txhi)-log10(txlo))/nm
		ltxm<- ((1:nm)-1)*ltsm+log10(txlo)
		txm<- 10^ltxm
		px<- qra_grbmodlc(txm,0.3,10.0,-1)
		mxt<- px$time
		mxrwt<- px$rates[1:1,1:nm]+px$rates[3:3,1:nm]
		mxrpc<- px$rates[2:2,1:nm]+px$rates[4:4,1:nm]
# plot model
		lines(mxt,mxrwt,type="l")
		lines(mxt,mxrpc,type="l")
	}
}
qra_grbaddpul<- function(grbfit,batxrt,tlo,thi) {
	if(batxrt==1) {
		cts2flux<- 1300
	} else {
		cts2flux<- 0.044
	}
#
	cat("Use cursor to mark start and peak of pulse.\n")
	X11()
	qra_grbplotsum(grbfit,batxrt,tlo,thi)
	st<- locator(1)
	pk<- locator(1)
	dev.off()
#
	ppar<- double(length=9)
	pnam<- character(length=9)
	ppar[1]<- -0.5
	pnam[1]<- "b1"
	ppar[2]<- 10
	pnam[2]<- "amb"
	ppar[3]<- -1
	pnam[3]<- "d"
	ppar[4]<- 1
	pnam[4]<- "a"
	ppar[5]<- 0.2
	pnam[5]<- "fr"
	ppar[6]<- (pk$x-st$x)/(1-ppar[5])
	pnam[6]<- "Tf"
#
	tpk<-pk$x
	zgrb<- grbfit$ZGRB
	tpkz<- tpk/(zgrb+1)
	if(tpkz>100) {
		#epk<- max(500*(tpkz/100)^(-1),10)/(zgrb+1)
		epk<- 500/(zgrb+1)
	} else {
		epk<- 500/(zgrb+1)
	}
	ppar[7]<- epk
	pnam[7]<- "Epk"
	ppar[8]<- tpk
	pnam[8]<- "Tpk"
	ppar[9]<- pk$y*cts2flux
	pnam[9]<- "Spk"
# search for place to add
	ns<- grbfit$NS
	np<- grbfit$NP
	nss<- ns+1
	npp<- np+9
	ip<-0
	if(ns>0) {
	    for(i in 1:ns) {
		ii<-(i-1)*9+2
		if(ip==0&(grbfit$PROMF[ii+7]>ppar[8])) {
			ip<- ii
		}
    	    }
	}
	if(ip==0) ip<- ns*9+2
	ii<- ip-1
	ipp<- ip+8
	ippp<-ip+9
# create new arrays
	promf<- double(length=npp)
	proml<- double(length=npp)
	promh<- double(length=npp)
	pname<- character(length=npp)
#
	promf[1:ii]<- grbfit$PROMF[1:ii]
	proml[1:ii]<- grbfit$PROML[1:ii]
	promh[1:ii]<- grbfit$PROMH[1:ii]
	pname[1:ii]<- grbfit$PNAME[1:ii]
	promf[1:1]<- nss
#
	promf[ip:ipp]<- ppar
	pname[ip:ipp]<- pnam
	proml[ip:ipp]<- 0
	promh[ip:ipp]<- 0
#
	if(npp>ipp) {
		promf[ippp:npp]<- grbfit$PROMF[ip:np]
		proml[ippp:npp]<- grbfit$PROML[ip:np]
		promh[ippp:npp]<- grbfit$PROMH[ip:np]
		pname[ippp:npp]<- grbfit$PNAME[ip:np]
	}
#
	grbfit$NS<- nss
	grbfit$NP<- npp
	grbfit$PROMF<- promf
	grbfit$PROML<- proml
	grbfit$PROMH<- promh
	grbfit$PNAME<- pname
#
	return(grbfit)
}
qra_grbdelpul<- function(grbfit,ipul) {
	ns<- grbfit$NS
	if(ipul<1|ipul>ns) {
		cat("qra_delpul error - pulse ",ipul," not set\n")
		return(grbfit)
	}
	np<- grbfit$NP
	nss<- ns-1
	npp<- np-9
	ip<-(ipul-1)*9+2
	ii<- ip-1
	ipp<- ip+9
# create new arrays
	promf<- double(length=npp)
	proml<- double(length=npp)
	promh<- double(length=npp)
	pname<- character(length=npp)
#
	promf[1:ii]<- grbfit$PROMF[1:ii]
	proml[1:ii]<- grbfit$PROML[1:ii]
	promh[1:ii]<- grbfit$PROMH[1:ii]
	pname[1:ii]<- grbfit$PNAME[1:ii]
	promf[1:1]<- nss
#
	promf[ip:npp]<- grbfit$PROMF[ipp:np]
	proml[ip:npp]<- grbfit$PROML[ipp:np]
	promh[ip:npp]<- grbfit$PROMH[ipp:np]
	pname[ip:npp]<- grbfit$PNAME[ipp:np]
#
	grbfit$NS<- nss
	grbfit$NP<- npp
	grbfit$PROMF<- promf
	grbfit$PROML<- proml
	grbfit$PROMH<- promh
	grbfit$PNAME<- pname
#
	return(grbfit)
}
qra_grbaddaft<- function(grbfit) {
	batxrt<- 2
	cts2flux<- 0.044
#
	cat("Use cursor to mark afterglow maximum and end of plateau.\n")
	X11()
	qra_grbplotsum(grbfit,batxrt,0,0)
	pk<- locator(1)
	en<- locator(1)
	dev.off()
#
	apar<- double(length=8)
	pnam<- character(length=8)
	apar[1]<- 100
	pnam[1]<- "na"
	apar[2]<- 0
	pnam[2]<- "vv"
	apar[3]<- -1.0
	pnam[3]<- "ba"
	apar[4]<- 0 
	pnam[4]<- "ga"
	apar[5]<- 100 
	pnam[5]<- "Tr"
	apar[6]<- en$x
	pnam[6]<- "Ta"
	apar[7]<- pk$y*cts2flux
	pnam[7]<- "Fm"
	apar[8]<- 1.5
	pnam[8]<- "aa"
# set place to add
	ns<- grbfit$NS
	np<- grbfit$NP
	nss<- ns
	npp<- ns*9+9
	ip<- ns*9+2
	ii<- ip-1
# create new arrays
	promf<- double(length=npp)
	proml<- double(length=npp)
	promh<- double(length=npp)
	pname<- character(length=npp)
#
	promf[1:ii]<- grbfit$PROMF[1:ii]
	proml[1:ii]<- grbfit$PROML[1:ii]
	promh[1:ii]<- grbfit$PROMH[1:ii]
	pname[1:ii]<- grbfit$PNAME[1:ii]
	promf[1:1]<- nss
#
	promf[ip:npp]<- apar
	pname[ip:npp]<- pnam
	proml[ip:npp]<- 0
	promh[ip:npp]<- 0
#
	grbfit$NS<- nss
	grbfit$NP<- npp
	grbfit$PROMF<- promf
	grbfit$PROML<- proml
	grbfit$PROMH<- promh
	grbfit$PNAME<- pname
#
	return(grbfit)
}
qra_grbadddust<- function(grbfit) {
	apar<- double(length=15)
	pnam<- character(length=15)
	apar[1]<- 200
	pnam[1]<- "na"
	apar[2]<- 1
	pnam[2]<- "od"
	apar[3]<- 0.01
	pnam[3]<- "am"
	apar[4]<- 0.25
	pnam[4]<- "ap"
	apar[5]<- 3.5
	pnam[5]<- "dd"
	apar[6]<- 0.347
	pnam[6]<- "zd"
	apar[7]<- 100
	pnam[7]<- "rl"
	apar[8]<- 10
	pnam[8]<- "rd"
	apar[9]<- 1
	pnam[9]<- "nr"
	apar[10]<- -1.0
	pnam[10]<- "ba"
	apar[11]<- 0
	pnam[11]<- "ga"
	apar[12]<- 100 
	pnam[12]<- "Tr"
	apar[13]<- 1000
	pnam[13]<- "Ta"
	apar[14]<- 0.0
	pnam[14]<- "Fm"
	apar[15]<- 1.5
	pnam[15]<- "aa"
# set place to add
	ns<- grbfit$NS
	np<- grbfit$NP
	nss<- ns
	npp<- ns*9+16
	ip<- ns*9+2
	ii<- ip-1
# create new arrays
	promf<- double(length=npp)
	proml<- double(length=npp)
	promh<- double(length=npp)
	pname<- character(length=npp)
#
	promf[1:ii]<- grbfit$PROMF[1:ii]
	proml[1:ii]<- grbfit$PROML[1:ii]
	promh[1:ii]<- grbfit$PROMH[1:ii]
	pname[1:ii]<- grbfit$PNAME[1:ii]
	promf[1:1]<- nss
#
	promf[ip:npp]<- apar
	pname[ip:npp]<- pnam
	proml[ip:npp]<- 0
	promh[ip:npp]<- 0
#
	grbfit$NS<- nss
	grbfit$NP<- npp
	grbfit$PROMF<- promf
	grbfit$PROML<- proml
	grbfit$PROMH<- promh
	grbfit$PNAME<- pname
#
	return(grbfit)
}
qra_grbaddbrk<- function(grbfit) {
	batxrt<- 2
# set place to add
	ns<- grbfit$NS
	np<- grbfit$NP
	nss<- ns
	npa<- ns*9+7
	npp<- np+2
	ip<- np+1
	ii<- ip-1
	if(np<npa) {
		cat("qra_grbaddbrk error - afterglow not set\n")
		return(grbfit)
	}
	if(npa==np) {
		ab<- 0.5
	} else {
		ab<- grbfit$PROMF[npa]+0.5
	}
#
	cat("Use cursor to mark time of break\n")
	X11()
	qra_grbplotsum(grbfit,batxrt,0,0)
	br<- locator(1)
	dev.off()
# create new arrays
	promf<- double(length=npp)
	proml<- double(length=npp)
	promh<- double(length=npp)
	pname<- character(length=npp)
#
	promf[1:ii]<- grbfit$PROMF[1:ii]
	proml[1:ii]<- grbfit$PROML[1:ii]
	promh[1:ii]<- grbfit$PROMH[1:ii]
	pname[1:ii]<- grbfit$PNAME[1:ii]
	promf[1:1]<- nss
#
	promf[ip:ip]<- br$x
	pname[ip:ip]<- "Tb"
	proml[ip:ip]<- 0
	promh[ip:ip]<- 0
	promf[npp:npp]<- ab
	pname[npp:npp]<- "ab"
	proml[npp:npp]<- 0
	promh[npp:npp]<- 0
#
	grbfit$NS<- nss
	grbfit$NP<- npp
	grbfit$PROMF<- promf
	grbfit$PROML<- proml
	grbfit$PROMH<- promh
	grbfit$PNAME<- pname
#
	return(grbfit)
}
qra_grbdelbrk<- function(grbfit) {
	ns<- grbfit$NS
	np<- grbfit$NP
	npa<- ns*9+7
	if(np<npa) {
		cat("qra_grdelbrk error - afterglow not set\n")
		return(grbfit)
	}
	if(np==npa) {
		cat("qra_grdelbrk error - break not set\n")
		return(grbfit)
	}
	npp<- np-2
#
	grbfit$NP<- npp
	grbfit$PROMF<- grbfit$PROMF[1:npp]
	grbfit$PROML<- grbfit$PROML[1:npp]
	grbfit$PROMH<- grbfit$PROMH[1:npp]
	grbfit$PNAME<- grbfit$PNAME[1:npp]
#
	return(grbfit)
}
qra_grbdelaft<- function(grbfit) {
	ns<- grbfit$NS
	np<- grbfit$NP
	npa<- ns*9+7
	if(np<npa) {
		cat("qra_grdelaft error - afterglow not set\n")
		return(grbfit)
	}
	npp<- ns*9+1
#
	grbfit$NP<- npp
	grbfit$PROMF<- grbfit$PROMF[1:npp]
	grbfit$PROML<- grbfit$PROML[1:npp]
	grbfit$PROMH<- grbfit$PROMH[1:npp]
	grbfit$PNAME<- grbfit$PNAME[1:npp]
#
	return(grbfit)
}
qra_grbfitaft<- function(grbfit) {
# Set up data masks
	bgood<- grbfit$MBAT&(grbfit$TBAT>grbfit$TBFLO)&(grbfit$TBAT<grbfit$TBFHI)
	xgood<- grbfit$MXRT&(grbfit$TXRT>grbfit$TXFLO)&(grbfit$TXRT<grbfit$TXFHI)
# Get parameters
	np<-grbfit$NP
	ns<-grbfit$NS
	zgrb<- grbfit$ZGRB
	pname<-grbfit$PNAME
#
	ii<-ns*9+1
	if(np == ii) {
		cat("No afterglow for ",grbfit$GRBN,"\n")
		return(grbfit)
	}
	cat("Fit afterglow of ",grbfit$GRBN,"\n")
# Select afterglow parameters
	nfl<-0
	ifl<-vector(mode="integer",length=nfl)
	derr<-vector(mode="logical",length=nfl)
	na<- grbfit$PROMF[ii+1]
	if(na==100) {
	        ii=ii+2
		#nfl<-nfl+1
		#ifl[nfl]<-ii
		#derr[nfl]<- T
	}
	if(na==200) {
	        ii=ii+2
# od dust column optical depth
		nfl<-nfl+1
		ifl[nfl]<-ii
		derr[nfl]<- T
# am dust minimum size
		ii=ii+1
# ap dust maximum size
		ii=ii+1
# dd dust size distribution index
	        ii=ii+1
		nfl<-nfl+1
		ifl[nfl]<-ii
		derr[nfl]<- T
# zd dust redshift
	        ii=ii+1
# rl minimum distance to dust
	        ii=ii+1
		nfl<-nfl+1
		ifl[nfl]<-ii
		derr[nfl]<- T
# rd distance range of dust layers
	        ii=ii+1
# nd number of dust layers
	        ii=ii+1
	}
# ba
	nfl<-nfl+1
	ii=ii+1
	ifl[nfl]<-ii
	derr[nfl]<- T
# ga
#	nfl<-nfl+1
#	ii=ii+1
#	ifl[nfl]<-ii
#	derr[nfl]<- T
# Tr
#	nfl<-nfl+1
#	ii=ii+1
#	ifl[nfl]<-ii
#	derr[nfl]<- T
# ta
	nfl<-nfl+1
	ii<-ii+3
	ifl[nfl]<-ii
	derr[nfl]<- T
# fa
	nfl<-nfl+1
	ii=ii+1
	ifl[nfl]<-ii
	derr[nfl]<- T
	if(grbfit$PROMF[ii]>0) {
# aa
		nfl<-nfl+1
		ii<-ii+1
		ifl[nfl]<-ii
		derr[nfl]<- T
		if(np>ii) {
			ibr<-(np-ii)/2
# Cycle through all late breaks set
			for(i in 1:ibr) {
# Tb late break
				nfl<-nfl+1
				ii<-ii+1
				ifl[nfl]<-ii
				derr[nfl]<- F
# ab
				nfl<-nfl+1
				ii<-ii+1
				ifl[nfl]<-ii
				derr[nfl]<- F
			}
		}
	}
# load parameters in common
	s<- qra_grbloadpars(grbfit$PROMF,grbfit$GRBN)
# Do fit
	b<-qra_grbsetfit(ifl,bgood,xgood)
	ndata<- b$bgood*4+b$xgood*2
	cat("Initial afterglow parameters\n")
	for(k in 1:nfl) {
		kk<- ifl[k]
		cat(k,kk,grbfit$PNAME[kk],
		b$spars[k]," range ",b$lpars[k],b$upars[k],"\n")
	}
	nfr<- sum(derr)
	delstat<-qchisq(0.9,nfr)

	f<- qr_srchmin(b$spars,b$lpars,b$upars,qra_grbchisq,delstat,derr)
	g<-qra_grbgetfit(1)
	cat("Fit afterglow with",g$ncall,"calls to qra_grbchisq\n")
	for(k in 1:nfl) {
		kk<- ifl[k]
		grbfit$PROMF[kk]<- f$par[k]
		grbfit$PROML[kk]<- f$parlo[k]
		grbfit$PROMH[kk]<- f$parhi[k]
		cat(k,kk,grbfit$PNAME[kk],grbfit$PROMF[kk]," range ",
		grbfit$PROML[kk],grbfit$PROMH[kk],"\n")
	}
	chisq<-f$value
	cat("Chi-squared",chisq,"with",ndata,"data points and",nfr,
	"free parameters\n")
	cat(b$ndof,"DoF, Reduced Chi-squared",chisq/b$ndof,"\n")
#
	grbfit$CHIMIN<- chisq
	grbfit$NFREE<- nfr
	grbfit$NDATA<- ndata
	grbfit$NDOF<- b$ndof
	return(grbfit)
}
qra_grbfitpul<- function(grbfit,ipul,efit) {
# Get parameters
	np<-grbfit$NP
	ns<-grbfit$NS
	zgrb<- grbfit$ZGRB
	pname<-grbfit$PNAME
#
	cat("Fit pulse ",ipul,"of ",grbfit$GRBN,"\n")
	if(efit) {
		cat("Epk floating, fr fixed\n")
	} else {
		cat("fr floating, Epk fixed\n")
	}
	nfl<-0
	ifl<-vector(mode="integer",length=nfl)
	derr<-vector(mode="logical",length=nfl)
	ii<-(ipul-1)*9+1
# b1
	nfl<-nfl+1
	ii<-ii+1
	ifl[nfl]<-ii
	derr[nfl]<- T
	b1<-grbfit$PROMF[ii]
# amb
	ii<-ii+1
# fr
	ii<-ii+3
	fr<- grbfit$PROMF[ii]
	if(!efit) {
		nfl<-nfl+1
		ifl[nfl]<-ii
		derr[nfl]<- F
	}
# tf
	nfl<-nfl+1
	ii<-ii+1
	tf<- grbfit$PROMF[ii]
	ifl[nfl]<-ii
	derr[nfl]<- T
# epk 
	ii<-ii+1
	epk<-grbfit$PROMF[ii]
	if(efit) {
		nfl<-nfl+1
		ifl[nfl]<-ii
		derr[nfl]<- F
	}
# Tpk
	ii<-ii+1
	tpk<-grbfit$PROMF[ii]
# Change Epk
        if(epk==200) {
                tpkz<- tpk/(zgrb+1)
		if(tpkz>100) {
                	#epkn<- max(500*(tpkz/100)^(-1),10)/(zgrb+1)
			epkn<- 500/(zgrb+1)
		} else {
			epkn<- 500/(zgrb+1)
		}
                grbfit$PROMF[ii-1]<- epkn
		cat(ipul,epk,epkn,"\n")
        }
# fpk 
	nfl<-nfl+1
	ii<-ii+1
	ifl[nfl]<-ii
	derr[nfl]<- T
# load parameters in common
	s<- qra_grbloadpars(grbfit$PROMF,grbfit$GRBN)
# Set up data masks
	tll<- tpk-tf*fr*2
	cat("tpk tll",tpk,tll,"\n")
	bgood<- grbfit$MBAT&(grbfit$TBAT>tll)&(grbfit$TBAT<grbfit$TBFHI)
	xgood<- grbfit$MXRT&(grbfit$TXRT>tll)&(grbfit$TXRT<grbfit$TXFHI)
# Do fit
	b<-qra_grbsetfit(ifl,bgood,xgood)
	cat("Initial pulse",ipul,"parameters\n")
	for(k in 1:nfl) {
		kk<- ifl[k]
		cat(k,kk,grbfit$PNAME[kk],
		b$spars[k]," range ",b$lpars[k],b$upars[k],"\n")
	}
	nfr<- sum(derr)
	delstat<-qchisq(0.9,nfr)
	f<- qr_srchmin(b$spars,b$lpars,b$upars,qra_grbchisq,delstat,derr)
	g<-qra_grbgetfit(1)
	cat("** Fit pulse",ipul,"with",g$ncall,"calls to qra_grbchisq\n")
	for(k in 1:nfl) {
		kk<- ifl[k]
		grbfit$PROMF[kk]<- f$par[k]
		grbfit$PROML[kk]<- f$parlo[k]
		grbfit$PROMH[kk]<- f$parhi[k]
		cat(k,kk,grbfit$PNAME[kk],grbfit$PROMF[kk]," range ",
		grbfit$PROML[kk],grbfit$PROMH[kk],"\n")
	}
	chisq<-f$value
	ndata<- b$bgood*4+b$xgood*2
	cat("Chi-squared",chisq,"with",ndata,"data points and",
	nfr,"free parameters\n")
	cat(b$ndof,"DoF, Reduced Chi-squared",chisq/b$ndof,"\n")
#
	grbfit$CHIMIN<- chisq
	grbfit$NFREE<- nfr
	grbfit$NDATA<- ndata
	grbfit$NDOF<-b$ndof
	return(grbfit)
}
qra_grbfitgrp<- function(grbfit,ilo,ihi) {
# Get parameters
	np<-grbfit$NP
	ns<-grbfit$NS
	zgrb<- grbfit$ZGRB
	pname<-grbfit$PNAME
#
	cat("Fit group ",ilo," to ",ihi,"of ",grbfit$GRBN,"\n")
	nfl<-0
	ifl<-vector(mode="integer",length=nfl)
	derr<-vector(mode="logical",length=nfl)
#
	for(ipul in ilo:ihi) {
		ii<-(ipul-1)*9+1
# b1
		ii<-ii+1
	#nfl<-nfl+1
	#ifl[nfl]<-ii
	#derr[nfl]<- F
	#b1<-grbfit$PROMF[ii]
# amb
		ii<-ii+1
# fr
		ii<-ii+3
		fr<- grbfit$PROMF[ii]
		nfl<-nfl+1
		ifl[nfl]<-ii
		derr[nfl]<- F
# tf
		nfl<-nfl+1
		ii<-ii+1
		tf<- grbfit$PROMF[ii]
		ifl[nfl]<-ii
		derr[nfl]<- F
# epk 
		ii<-ii+1
		epk<-grbfit$PROMF[ii]
# Tpk
		ii<-ii+1
		tpk<-grbfit$PROMF[ii]
		if(ipul==ilo) {
			tss<- tpk
		}
# fpk 
		nfl<-nfl+1
		ii<-ii+1
		ifl[nfl]<-ii
		derr[nfl]<- F
	}
# load parameters in common
	s<- qra_grbloadpars(grbfit$PROMF,grbfit$GRBN)
# Set up data masks
	tll<- tss-tf*fr*2
	cat("tss tll",tpk,tll,"\n")
	bgood<- grbfit$MBAT&(grbfit$TBAT>tll)&(grbfit$TBAT<grbfit$TBFHI)
	xgood<- grbfit$MXRT&(grbfit$TXRT>tll)&(grbfit$TXRT<grbfit$TXFHI)
# Do fit
	b<-qra_grbsetfit(ifl,bgood,xgood)
	cat("Initial group",ilo," to ",ihi,"parameters\n")
	for(k in 1:nfl) {
		kk<- ifl[k]
		cat(k,kk,grbfit$PNAME[kk],
		b$spars[k]," range ",b$lpars[k],b$upars[k],"\n")
	}
	nfr<- sum(derr)
	delstat<-qchisq(0.9,nfr)
	f<- qr_srchmin(b$spars,b$lpars,b$upars,qra_grbchisq,delstat,derr)
	g<-qra_grbgetfit(1)
	cat("** Fit group",ilo," to ",ihi,"with",g$ncall,
	"calls to qra_grbchisq\n")
	for(k in 1:nfl) {
		kk<- ifl[k]
		grbfit$PROMF[kk]<- f$par[k]
		grbfit$PROML[kk]<- f$parlo[k]
		grbfit$PROMH[kk]<- f$parhi[k]
		cat(k,kk,grbfit$PNAME[kk],grbfit$PROMF[kk]," range ",
		grbfit$PROML[kk],grbfit$PROMH[kk],"\n")
	}
	chisq<-f$value
	ndata<- b$bgood*4+b$xgood*2
	cat("Chi-squared",chisq,"with",ndata,"data points and",
	nfr,"free parameters\n")
	cat(b$ndof,"DoF, Reduced Chi-squared",chisq/b$ndof,"\n")
#
	grbfit$CHIMIN<- chisq
	grbfit$NFREE<- nfr
	grbfit$NDATA<- ndata
	grbfit$NDOF<-b$ndof
	return(grbfit)
}
qra_grbfitpars<- function(grbfit,ifl) {
# Set up data masks
	bgood<- grbfit$MBAT&(grbfit$TBAT>grbfit$TBFLO)&(grbfit$TBAT<grbfit$TBFHI)
	xgood<- grbfit$MXRT&(grbfit$TXRT>grbfit$TXFLO)&(grbfit$TXRT<grbfit$TXFHI)
# Get parameters
	np<- grbfit$NP
	pname<-grbfit$PNAME
#
	cat("Fit parameters of ",grbfit$GRBN,"\n")
	nfl<- length(ifl)
	derr<-vector(mode="logical",length=nfl)
# Check parameter indices
	for(k in 1:nfl) {
		if(ifl[k]<2 | ifl[k]>np) {
			cat("qra_grbfitpars - error parameter index",ifl[k],
			"out of range\n")
			return(grbfit)
		}
		derr[k]<- FALSE
	}
# load parameters in common
	s<- qra_grbloadpars(grbfit$PROMF,grbfit$GRBN)
# Do fit
	b<-qra_grbsetfit(ifl,bgood,xgood)
	cat("Initial parameter values and allowed ranges\n")
	for(k in 1:nfl) {
		kk<- ifl[k]
		cat(k,kk,grbfit$PNAME[kk],
		b$spars[k]," range ",b$lpars[k],b$upars[k],"\n")
	}
	nfr<- nfl
	delstat<-qchisq(0.9,nfr)
	f<- qr_srchmin(b$spars,b$lpars,b$upars,qra_grbchisq,delstat,derr)
	g<-qra_grbgetfit(1)
	cat("** Parameters fitted with",g$ncall,"calls to qra_grbchisq\n")
	for(k in 1:nfl) {
		kk<- ifl[k]
		grbfit$PROMF[kk]<- f$par[k]
		grbfit$PROML[kk]<- f$parlo[k]
		grbfit$PROMH[kk]<- f$parhi[k]
		cat(k,kk,grbfit$PNAME[kk],grbfit$PROMF[kk],"\n")
	}
	chisq<-f$value
	ndata<- b$bgood*4+b$xgood*2
	cat("Chi-squared",chisq,"with",ndata,"data points and",
	nfr,"free parameters\n")
	cat(b$ndof,"DoF, Reduced Chi-squared",chisq/b$ndof,"\n")
#
	grbfit$CHIMIN<- chisq
	grbfit$NFREE<- nfr
	grbfit$NDATA<- ndata
	grbfit$NDOF<-b$ndof
	return(grbfit)
}
# X-ray astrophysics routines - prefix qrx_
qrx_ismtau<-function(nh,z,ear) {
	ie<-length(ear)
	ibin<-ie-1
	ebin<-(ear[1:ibin]+ear[2:ie])/2
	a<-.Fortran("qrx_ismtau", DUP=T,
	as.single(nh),
	as.single(z),
	as.integer(ibin),
	as.single(ear),
	tau=single(length=ibin),
	etau=single(length=ibin))
	return(list(nh=nh,z=z,ear=ear,ebin=ebin,tau=a$tau))
}
qrx_iismtau<-function(nh,z,tk,pl,ist,ear) {
	ie<-length(ear)
	ibin<-ie-1
	ebin<-(ear[1:ibin]+ear[2:ie])/2
	a<-.Fortran("qrx_iismtau", DUP=T,
	as.single(nh),
	as.single(z),
	as.single(tk),
	as.single(pl),
	as.single(ist),
	as.integer(ibin),
	as.single(ear),
	tau=single(length=ibin),
	etau=single(length=ibin))
	return(list(nh=nh,z=z,tk=tk,pl=pl,ist=ist,ear=ear,ebin=ebin,tau=a$tau))
}
qrx_brems<- function(ekev,t) {
	ne<- length(ekev)
	a<- .Fortran("qrx_brems",
	as.integer(ne),
	as.double(ekev),
	as.double(t),
	ph=double(length=ne))
	return(a$ph)
}
}
