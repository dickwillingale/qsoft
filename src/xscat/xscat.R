# qR xscat routines
if(!exists("qrs_init")) {
qrs_init<-function() {
	.Fortran("qrs_init")
	invisible()
}
qrs_init()
qrt_xopt<-function(mspec,rho,ekev,itype) {
	.Fortran("qrt_xopt",
	as.integer(nchar(mspec)),
	as.character(mspec),
	as.double(rho),
	as.integer(length(ekev)),
	as.double(ekev),
	as.integer(itype),
	alpha=double(length(ekev)),
	gamma=double(length(ekev)),
	absl=double(length(ekev)),
	f1=double(length(ekev)),
	f2=double(length(ekev)))
}
qrt_xfresnel<-function(alpha,gamma,angs) {
	.Fortran("qrt_xfresnel",
	as.double(alpha),
	as.double(gamma),
	as.integer(length(angs)),
	as.double(angs),
	rs=double(length(angs)),
	rp=double(length(angs)),
	runp=double(length(angs)))
}
qrt_mlayer<-function(angs,ekev,nr,ni,d,nper) {
	.Fortran("qrt_mlayer",
	as.integer(length(angs)),
	as.double(angs),
	as.double(ekev),
	as.integer(length(nr)),
	as.double(nr),
	as.double(ni),
	as.double(d),
	as.integer(nper),
	rsig=double(length(angs)),
	rpi=double(length(angs)),
	tsig=double(length(angs)),
	tpi=double(length(angs)))
}
}
