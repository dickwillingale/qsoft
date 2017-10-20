# qsoft biomark routines
if(!exists("qrb_init")) {
qrb_init<-function() {
	.Fortran("qrb_init")
	invisible()
}
qrb_init()
# Bin up raw mass spec data
qrb_biobin<- function(x,w,xleft,xright,nx) {
        n<- length(x)
        a<-.Fortran("qrb_biobin",
        as.integer(n),
        as.double(x),
        as.double(w),
        as.double(xleft),
        as.double(xright),
        as.integer(nx),
        flux=double(length=nx),
        mz=double(length=nx))
        return(list(mz=a$mz,flux=a$flux))
}
# Perform background substraction of binned mass spectrum
qrb_bioback<- function(w,nback,nfrac) {
        np<- nrow(w)
        ns<- ncol(w)
        nn<- np*ns
        a<-.Fortran("qrb_bioback",
        as.integer(np),
        as.integer(ns),
        as.double(w),
        as.integer(nback),
        double(length=nback),
        double(length=nback),
        as.integer(nfrac),
        pk=double(length=nn),
        bk=double(length=nn),
        va=double(length=nn))
        return(list(peaks=a$pk,back=a$bk,var=a$va))
}
# Perform smoothing of binned mass spectrum
qrb_biospread<- function(w,gsig) {
        np<- length(w)
        ns<- length(gsig)
        a<- .Fortran("qrb_biospread",
        as.integer(np),
        as.double(w),
        double(length=np),
        as.integer(ns),
        as.double(gsig),
        sp=double(length=np))
        return(a$sp)
}
# Find significant peaks in the summed spectrum
qrb_biosigpks<- function(minsig,gsig,mz,fl,va) {
	n<- length(mz)
	ng<- length(gsig)
	a<- .Fortran("qrb_biosigpks",
	as.double(minsig),
	as.integer(ng),
	as.double(gsig),
	as.integer(n),
	as.double(mz),
	as.double(fl),
	as.double(va),
	double(length=n),
	ip=integer(length=n),
	pmz=double(length=n),
	pfl=double(length=n),
	psi=double(length=n),
	np=integer(length=1))
        return(list(np=a$np,ip=a$ip[1:a$np],pmz=a$pmz[1:a$np],
	pfl=a$pfl[1:a$np],psi=a$psi[1:a$np]))
}
# Extract peaks from all samples
qrb_biopeaks<- function(ipeak,flx,var,bak,mz) {
	nb<- nrow(flx)
	ns<- ncol(flx)
	np<- length(ipeak)
	ihwid<- integer(length=nb)
	ihwid[1:nb]<- 1
#
	a<- .Fortran("qrb_biopeaks",
	as.integer(nb),
	as.integer(ns),
	as.integer(np),
	as.integer(ipeak),
	as.integer(ihwid),
	as.double(flx),
	as.double(var),
	as.double(bak),
	as.double(mz),
	cflx=double(length=np),
	ppos=double(length=np),
	pflx=double(length=np*ns),
	pvar=double(length=np*ns),
	pbak=double(length=np*ns))
	dim(a$pflx)<- c(np,ns)
	dim(a$pvar)<- c(np,ns)
	dim(a$pbak)<- c(np,ns)
	return(list(npeak=np,nsample=ns,mz=mz,cflx=a$cflx,ppos=a$ppos,
               pflx=a$pflx,pvar=a$pvar,pbak=a$pbak))
}
# Scale peaks from all samples
qrb_bioscale<- function(inorm,cflx,pflx,pvar,pbak) {
	npn<- length(inorm)
	ns<- ncol(pflx)
	np<- nrow(pflx)
#
	a<- .Fortran("qrb_bioscale",
	as.integer(npn),
	as.integer(inorm),
	as.integer(ns),
	as.integer(np),
	as.double(cflx),
	sflx=as.double(pflx),
	svar=as.double(pvar),
	sbak=as.double(pbak),
	rnorm=double(length=ns),
	imax=integer(length=np),
	pcts=double(length=np),
	psig=double(length=np),
	pnd=double(length=np*ns))
	dim(a$pnd)=c(np,ns)
	dim(a$sflx)=c(np,ns)
	dim(a$svar)=c(np,ns)
	dim(a$sbak)=c(np,ns)
	return(list(rnorm=a$rnorm,sflx=a$sflx,svar=a$svar,sbak=a$sbak,
		imax=a$imax,pcts=a$pcts,psig=a$psig,pnd=a$pnd))
}
}
