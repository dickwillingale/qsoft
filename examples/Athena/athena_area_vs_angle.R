#!/usr/bin/env Rscript
# Athena vignetting function
args <- commandArgs(trailingOnly = TRUE)
if(length(args)!=1) {
	cat("Innvocation: athena_area_vs_angle.R kev\n")
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
# Set position of aperture centre
a2j<-max(lm)
pcen<-c(flen+a2j,0,0)
# set source to cover aperture
ssize<-rmax*1.05
smin<-rmin*0.95
slim<-c(smin,ssize)
spos<-c(flen+300,0,0)
nray=500000
# Detector same size as source aperture
drad=ssize
dpos<-c(0,0,0)
dlim<-c(0,drad)
radet<-0.0
# set deformation index
id<-1
# set image parameters - pixel 0.25 arc secs
dpix<- (0.25/3600)*pi/180*flen
irad<-3
nx<- trunc(irad/dpix*2)
ny<- nx
irad<- nx*dpix/2
print(nx)
# Loop for offaxis angles
offarcmin<-seq(from=0,to=30,by=2)
na<-length(offarcmin)
area<- 1:na
hew<- 1:na
fwhm<- 1:na
dshft<- 1:na
rmsa<- 1:na
rmsb<- 1:na
tha<- 1:na
denx<- 1:na
deny<- 1:na
#
cat("offarcmin area hew fwhm dshft rmsa rmsb tha denx deny\n")
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
 dshft[i]=results$dshft
 detpos<-read.table("detected.dat",header=TRUE)
 #rad<- flen*max(gr)*offrad/8
 #xcen<- offrad*flen+rad
 xcen<- offrad*flen
 ycen<- 0
 xmin<- xcen-irad
 xmax<- xcen+irad
 ymin<- ycen-irad
 ymax<- ycen+irad
 a<-qri_binxy(detpos$YDET,detpos$ZDET,0,detpos$AREA,xmin,xmax,nx,ymin,ymax,ny)
 bb<-qri_beam(a$data_array,nx,0,0)
 area[i]<-bb$flux/100.
 hew[i]<-(bb$hewc*a$xsam/flen)*(180/pi)*3600
 fwhm[i]<-(bb$fwhmc*a$xsam/flen)*(180/pi)*3600
 rmsa[i]<-(bb$rmsa*a$xsam/flen)*(180/pi)*3600
 rmsb[i]<-(bb$rmsb*a$xsam/flen)*(180/pi)*3600
 tha[i]<- bb$tha
 denx[i]<- (bb$cen[1]-nx/2)*a$xsam
 deny[i]<- (bb$cen[2]-ny/2)*a$xsam
 cat(offarcmin[i],area[i],hew[i],fwhm[i],dshft[i],rmsa[i],rmsb[i],tha[i],
 denx[i],deny[i],"\n")
}
# Plot results
pdfile<- paste("ladapt_area_vs_angle_",kev,"_kev.pdf",sep="")
pdf(file=pdfile)
header<- paste("Energy",kev,"keV",sep=" ")
par(mfrow=c(2,2))
plot(offarcmin,area,type='l',ylim=c(0,max(area)),xlab="arc mins",
ylab=expression(area~~cm^2),main=header)
plot(offarcmin,hew,type='l',ylim=c(0,max(hew)),xlab="arc mins",
ylab="arc sec",main=header)
#plot(offarcmin,dshft,type='l',ylim=c(min(dshft),max(dshft)),xlab="arc mins",
#ylab="detector shift mm")
dev.off()
#go_on<- locator(1)
# tabulate as ascii file
dfile<- paste("ladapt_area_vs_angle_",kev,"_kev.dat",sep="")
write.table(
list(arcmin=offarcmin,area=area,
hew=hew,rmsa=rmsa,rmsb=rmsb,
tha=tha,denx=denx,deny=deny),
file=dfile,row.names=F,quote=F)
