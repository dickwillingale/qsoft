#!/usr/bin/env Rscript
# Athena area vs. energy
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
ist<-2
# Set optical constants for pure Ir (used if ist=1)
Ir<-qrt_xopt(mspec,rho,ekev,1)
alpha<-Ir$alpha
gamma<-Ir$gamma
#ist<- 1
# Set position of aperture centre
a2j<-max(lm)
pcen<-c(flen+a2j,0,0)
# set source
ssize<-rmax*1.05
smin<-rmin*0.95
slim<-c(smin,ssize)
spos<-c(flen+300,0,0)
nray=50000
# Detector
drad<- ssize
dpos<-c(0,0,0)
dlim<-c(0,drad)
radet<-0.0
# set deformation index
id<-1
# set image parameters
irad<-5
nx<-120
ny<-120
# Loop for energies
ne<-length(ekev)
area<- 1:ne
hew<- 1:ne
#
layout(1)
for(i in 1:ne) {
 offrad<-0
 ccy<-sin(offrad)
 di<-c(-sqrt(1.0-ccy^2),ccy,0.0)
 qrt_reset()
 qrt_surface(is,ist,ekev[i],nrough,brk,rind,alpha[i],gamma[i],angs,
 refs[1:nag,i:i],0,0,0)
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
 qr_rseed(88)
 results<-qrt_trace(0,drad,-2)
 detpos<-read.table("detected.dat",header=TRUE)
 xmin<- offrad*flen-irad
 xmax<- offrad*flen+irad
 ymin<- -irad
 ymax<- irad
 a<-qri_binxy(detpos$YDET,detpos$ZDET,0,detpos$AREA,xmin,xmax,nx,ymin,ymax,ny)
 bb<-qri_beam(a$data_array,nx,0,0)
 area[i]<-bb$flux/100.
 hew[i]<-(bb$hew*a$xsam/flen)*(180/pi)*3600
 cat(ekev[i],area[i],hew[i],"\n")
}
# Plot results
pdf(file="ladapt_area_vs_kev.pdf")
plot(ekev,area,type='l',log="y",ylim=c(50,max(area)),xlab="keV",
ylab=expression(area~~cm^2))
dev.off()
#go_on<- locator(1)
# tabulate as ascii file
write.table(list(ekev=ekev,area=area,hew=hew),
file="ladapt_area_vs_kev.dat",row.names=F,quote=F)
