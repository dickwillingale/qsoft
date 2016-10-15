#!/usr/bin/env Rscript
# Athena ray tracing - on-axis PSF
args <- commandArgs(trailingOnly = TRUE)
if(length(args)!=1) {
	        cat("Innvocation: athena_onaxis_psf.R kev\n")
        quit()
}
kev<- as.numeric(args[1])
# useful vectors
sn<-c(1,0,0)
nn<-c(0,0,0)
rx<-c(0,1,0)
# Mirror module parameters
load("athena_ladapt_modules.RData")
# reflecting surface definition
is<-1
# Use look up table of reflectivities of Ir + BC4 overcoat
load("ir_b4cover.RData")
angs<-ir_b4cover$angs
ekev<-ir_b4cover$ekev
refs<-ir_b4cover$refs
nag<-length(angs)
iekev=which(ekev==kev)
if(length(iekev)==0) {
        cat("Energy",kev,"keV not available\n")
        quit()
}
ist<-2
cat("index energy",iekev,ekev[iekev],"\n")
# Use optical constants
xop<-qrt_xopt(mspec,rho,ekev,1)
alpha<-xop$alpha
gamma<-xop$gamma
#ist<-1
cat("surface roughness",nrough,rind,brk,"\n")
# Set position of aperture centre
a2j<-max(lm)
pcen<-c(flen+a2j,0,0)
# set source to cover aperture
ssize<-rmax*1.05
smin<-rmin*0.95
slim<-c(smin,ssize)
spos<-c(flen+300,0,0)
nray=500000
# Detector
drad=ssize
dpos<-c(0,0,0)
dlim<-c(0,drad)
radet<-0.0
# set deformation index
id<-1
# set image parameters - pixel 0.25 arc secs
dpix<- (0.25/3600)*pi/180*flen
irad<-4
nx<- trunc(irad/dpix*2)
ny<- nx
irad<- nx*dpix/2
cat("image parameters",dpix,nx,irad,"\n")
# Loop for offaxis angles
offarcmin<-0
na<-length(offarcmin)
area<- 1:na
hew<- 1:na
fwhm<- 1:na
dshft<- 1:na
w90<- 1:na
#
cat("offarcmin area  hew w90 fwhm dshift\n")
#
#proc.time()
for(i in 1:na) {
 offrad<-(offarcmin[i]/60)*pi/180.
 ccy<-sin(offrad)
 di<-c(-sqrt(1.0-ccy^2),ccy,0.0)
 qrt_reset()
 qr_rseed(55)
 qrt_surface(is,ist,ekev[iekev],nrough,brk,rind,alpha,gamma,angs,
 refs[1:nag,iekev:iekev],0,0,0)
 qrt_deform(id,1,5,nm,1)
 qrt_defmat(id,1,xd,yd,dx)
 qrt_defmat(id,2,xd,yd,dy)
 qrt_defmat(id,3,xd,yd,dz)
 qrt_defmat(id,4,xd,yd,db)
 qrt_defmat(id,5,xd,yd,dc)
 qrt_sipore(pcen,sn,rx,flen,rpitch,apitch,wall,rm,pm,te,wm,hm,lm,cm,gm,wfr,
 a2j,id,is)
 qrt_source(1,di,nn,spos,sn,rx,slim,0,nray,0)
 qrt_detector(1,dpos,sn,rx,dlim,radet)
 results<-qrt_trace(0,drad,-2)
 proc.time()
 dshft[i]=results$dshft
 detpos<-read.table("detected.dat",header=TRUE)
 rad<- flen*max(gr)*offrad/8
 xcen<- offrad*flen+rad
 ycen<- 0
 xmin<- xcen-irad
 xmax<- xcen+irad
 ymin<- ycen-irad
 ymax<- ycen+irad
 a<-qri_binxy(detpos$YDET,detpos$ZDET,0,detpos$AREA,xmin,xmax,nx,ymin,ymax,ny)
 bb<-qri_beam(a$data_array,nx,0,0)
 #print(bb)
 area[i]<-bb$flux/100.
 hew[i]<-(bb$hew*a$xsam/flen)*(180/pi)*3600
 w90[i]<-(bb$w90*a$xsam/flen)*(180/pi)*3600
 fwhm[i]<-(bb$fwhm*a$xsam/flen)*(180/pi)*3600
 cat(offarcmin[i],area[i],hew[i],w90[i],fwhm[i],dshft[i],"\n")
}
#proc.time()
# create radial surface brightness distribution
radius<-sqrt(detpos$YDET^2+detpos$ZDET^2)
theta<-atan2(detpos$ZDET,detpos$YDET)
nr<- 100
rasbin<- 1
rbin<- rasbin/3600*pi/180*flen
rasmax<- nr*rasbin
rmax<- rasmax/3600*pi/180*flen
p<-qri_binxy(radius,theta,0,detpos$AREA,0,rmax,nr,-pi,pi,1)
ras<- seq(from=rasbin, to=rasmax, by=rasbin)
rm<- ras/3600*pi/180*flen
sb<-approx(p$xp,p$data_array/(2*pi*p$xp*rbin),rm,rule=2)
# Put surface brightness profile on file
#tt<- list(ras=ras,sb=sb$y/sb$y[1])
#write.table(tt,file="athena_ladapt_psf_15defocus.dat",col.names=T,row.names=F)
# Plot results
pdfile=paste("ladapt_PSF_",kev,"keV.pdf",sep=" ")
pdf(file=pdfile)
#par(mfrow=c(2,1))
a$xlab<- "mm"
a$ylab<- "mm"
qri_displayimage(a)
#go_on<- locator(1)
plot(ras,sb$y,log="y",type="l",
main=paste("ladapt PSF",kev,"keV",sep=" "),
xlab="radius arc secs",ylab="surface brightness")
#
cum<- cumsum(sb$y)
cum<- cum/max(cum)
plot(ras,cum,type="l",xlab="radius arc secs",ylab="flux fraction")
#go_on<- locator(1)
dev.off()
