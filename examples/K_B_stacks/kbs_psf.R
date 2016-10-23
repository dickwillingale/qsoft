#!/usr/bin/env Rscript
# Point spread function of Kirpatrick-Baez stacks
# useful vectors
sn<-c(1,0,0)
nn<-c(0,0,0)
rx<-c(0,1,0)
# reflecting surface definition
rough<-0
rind<-1.4
brk<-1.0
ekev<-1.0
is<-1
Ir<-qrt_xopt("Ir",22.65,ekev,1)
alpha<-Ir$alpha
gamma<-Ir$gamma
# Deformations
id<-0
# K-B stack parameters
fl=5000
#ipack<-1
#rmin<-244
ipack<-3
rmin<-0
rmax<-1900
pnor<-c(1,0,0)
raxi<-c(0,1,0)
pd<-0.610
wall<-0.15
pitch<-pd+wall
csize<-(rmax-rmin)/10
grdeg<-1.0
gr<-grdeg*pi/180
pl<-pd/gr
cat("grazing angle",grdeg,"degrees, axial length",pl,"mm\n")
pcen1<-c(fl+pl,0,0)
pcen2<-c(fl,0,0)
# set off-axis direction
offarcmin<-0.0
offrad<-(offarcmin/60)*pi/180.
ccy<-sin(offrad)
di<-c(-sqrt(1.0-ccy**2),-ccy,0.0)
rlim<-c(rmin-10,rmax+10)
spos<-c(fl+pl+3,0,0)
nray=1000000
# Detector - spherical with radius of curvature equal focal length
dnorm<--di
dpos<-c(-fl,0,0)+dnorm*fl
dlim<-c(0,rmax)
radet<-fl
halfd<-4
#
qrt_reset()
s1<-qrt_kbs(pcen1,pnor,raxi,ipack,rmin,rmax,fl,csize,pitch,wall,pl,pl,id,is)
s2<-qrt_kbs(pcen2,pnor,raxi,ipack,rmin,rmax,-fl,csize,pitch,wall,pl,pl,id,is)
qrt_source(1,di,nn,spos,sn,rx,rlim,0,nray,0)
qrt_detector(3,dpos,dnorm,rx,dlim,radet)
qrt_surface(is,1,ekev,rough,brk,rind,alpha,gamma,0,0,0,0,0)
results<-qrt_trace(0,100,2)
cat("detector axial shift",results$dshft,"\n")
# Get detected postions
detpos<-read.table("detected.dat",header=TRUE)
#plot(detpos$YDET,detpos$ZDET,type="p",pch=20)
# create image
xmin<- -halfd
xmax<- halfd
ymin<- -halfd
ymax<- halfd
nx<-100
ny<-100
a<-qri_binxy(detpos$YDET,detpos$ZDET,0,detpos$AREA,xmin,xmax,nx,ymin,ymax,ny)
#
X11()
qri_displayimage(a)
#,0,0,c(0,a$zlim/10))
#persp(x,y,a$data_array,asp=1)
qri_setfield(nx,xmin,xmax,ny,ymin,ymax)
rbeam<-pd*5
ibeam<-rbeam/a$xsam
b<-qri_beam(a$data,ibeam,0,0)
hew<-(b$hew*a$xsam/fl)*(180/pi)*3600
cat("HEW arc sec",hew,"\n")
fwhm<-(b$fwhm*a$xsam/fl)*(180/pi)*3600
cat("FWHM arc sec",fwhm,"\n")
area<-b$flux/100
cat("area",area,"cm2\n")
psfarea<-(b$hew*a$xsam/10)^2
fgain<-area/psfarea
cat("focusing gain",fgain,"\n")
go_on<- locator(1)
