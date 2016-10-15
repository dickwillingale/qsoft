#!/usr/bin/env Rscript
# Generate baseline sector modules for Athena
rmin<-258.0
rmax<- 1450.0
# rib spacing, radial and azimuthal pore pitch
wall<-0.17
pr<-0.605
rpitch<-pr+wall
apitch<- 2.3
# Axial curvature and principle plane
# Conical approximation
#iax<-0
# Wolter I
#iax<- 1
# Wolter I with spherical principal surface
#iax<- 4
# Conical approximation with spherical principal surface
#iax<- 5
# Wolter I with spherical principal surface and grid
#iax<- 6
# Spherical principal plane with curved-plane curvature
#iax<- 7
# Spherical principal plane with plane-curved curvature
iax<- 8
# Spherical principal plane with -1/2, 5/2 curvature
#iax<- 9
# Spherical principal plane with 1/2, 3/2 curvature
#iax<- 10
# Spherical principal plane with -1, 3 curvature
#iax<- 11
# Spherical principal plane with -1.5, 3.5 curvature
#iax<- 12
# Spherical principal plane with -2, 4 curvature
#iax<- 13
# reflecting surface roughness definition
rough<- 5.0
rind<-1.4
brk<- 1.
gam<- rind-1.0
nrough<- rough^2/((1+1/gam)*brk^(-gam))
mspec<- "Ir"
rho<- 22.65
# radial gap and module aperture height
gapr<- 8.0
mr<-54.1
# step size in radius
dring<-mr+gapr
# azimuthal gap
gapw<- 15.0
# width of spider arms
gaps<- 20.0
# module aperture width outer
dang<- 100.0
# module aperture width inner
dangi<- 60.0
# module frame width (around module aperture)
wfr<- 1.0
# focal length
flen<-12000.0
# ratio of grazing angles theta_h=theta_p*gra
gra<- 1.0
# Set up modules in rings
nring<-trunc((rmax-rmin)/dring)
rring<- double(length=nring)
hring<- double(length=nring)
graz<- double(length=nring)
graz1<- double(length=nring)
graz2<- double(length=nring)
imod<-0
nmod<-1
im<-1:nmod
ir<-1:nmod
is<-1:nmod
xm<-1:nmod
ym<-1:nmod
wm<-1:nmod
hm<-1:nmod
lm<-1:nmod
tm<-1:nmod
am<-1:nmod
cm<-1:nmod
gm<-1:nmod
rsplit<- 1:nring
for(i in 1:nring) {
	ring<-rmin+(i-0.5)*dring
	rring[i]<- ring
# Work in 60 degree sectors
	if(ring<500) {
		na<-trunc((ring*pi/3-gaps)/dangi)+1
	} else {
		na<-trunc((ring*pi/3-gaps)/dang)+1
	}
	dth<-(pi/3-gaps/ring)/na
	mw<- dth*ring-gapw
	graz[i]=atan(ring/flen)/2.0/(1.0+gra)
	graz1[i]=atan((ring-mr/2)/flen)/2.0/(1.0+gra)
	graz2[i]=atan((ring+mr/2)/flen)/2.0/(1.0+gra)
	ml<-pr/graz[i]
	hring[i]<- ml
	for(k in 1:na) {
		for(l in 0:5) {
			theta<-(k-1)*dth+l*pi/3+gaps/ring/2+mw/ring/2
			imod<-imod+1
			im[imod]<-imod
			ir[imod]<-i
			is[imod]<-k
			xm[imod]<-ring*cos(theta)
			ym[imod]<-ring*sin(theta)
			wm[imod]<-mr
			hm[imod]<-mw
			lm[imod]<-ml
			tm[imod]<-theta
			cm[imod]<- iax
			gm[imod]<- gra
		}
	}
	rsplit[i]<- rring[i]+dring/2
	cat(i,rsplit[i],graz[i]*180/pi*60,graz1[i]*180/pi*60,
	graz2[i]*180/pi*60,"\n")
}
# Radius of module aperture centres
rm<-sqrt(xm^2+ym^2)
# azimuthal position of module centres
pm<-atan2(ym,xm)
# Add frame width
wm<-wm+2*wfr
hm<-hm+2*wfr
# grazing angle at module centres
gr<-atan(rm/flen)/4
#
nm<-imod
cat(nm,"modules\n")
cat("maximum radius",max(rm)+mr/2,"\n")
cat("minimum radius",min(rm)-mr/2,"\n")
# set up deformation vectors that give module alignment and pore figure errors
xd<-1:nm
yd<-1
# module rotation about optical axis radians (3 arc secs)
rrms<-(3/3600)*pi/180
# shift of module in aperture plane rms mm (applied to x and y in aperture)
xyrms<-0.03
# axial module shift rms mm (includes error in kink angle)
zrms<-1.0
# in-plane figure errors rms radians (1.2 arc secs)
brms<-(1.2/3600)*pi/180
# out-of-plane figure errors rms radians (1.5 arc secs)
crms<-(1.5/3600)*pi/180
# set up deformation vectors for modules
te<-rnorm(nm,mean=0,sd=rrms)
dx<-rnorm(nm,mean=0,sd=xyrms)
dy<-rnorm(nm,mean=0,sd=xyrms)
dz<-rnorm(nm,mean=0,sd=zrms)
db<-rep(brms,nm)
dc<-rep(crms,nm)
# save parameters
save(rmin,rmax,flen,pr,iax,nrough,rind,brk,nm,rpitch,apitch,wall,wfr,
im,ir,is,xm,ym,wm,hm,lm,tm,am,cm,gm,rm,pm,gr,xd,yd,te,dx,dy,dz,db,dc,mspec,rho,
rring,dring,hring,graz,graz1,graz2,
file="athena_ladapt_modules.RData")
#
write.table(list(im=im,ir=ir,is=is,xm=xm*0.001,ym=ym*0.001,wm=wm,hm=hm,lm=lm,
rm=rm),file="athena_ladapt_modules.dat",row.names=F,quote=F)
# estimate volume
vol<-sum(wm*hm*lm)/1000
cat("volume cm3",vol,"\n")
cat("min max width",min(wm),max(wm),"\n")
cat("min max height",min(hm),max(hm),"\n")
cat("min max length",min(lm),max(lm),"\n")
#
pdf("athena_ladapt_modules.pdf")
plot.new()
hwid<-rmax*1.05
plot.window(c(-hwid,hwid),c(-hwid,hwid),asp=1)
qr_rectangles(xm[1:imod],ym[1:imod],wm[1:imod]/2,hm[1:imod]/2,tm[1:imod])
qr_circles(c(0,0),c(0,0),c(rmin,rmax))
par(col="red")
cen<- rep(0,nring)
qr_circles(cen,cen,rsplit)
axis(1)
axis(2)
dev.off()
# 
pdf("axial_length_ladapt.pdf")
par(mfrow=c(4,1))
hwid<-rmax*1.05
lwid<- max(hring)*1.05
plot(rmax*2,0,type="p",xlim=c(-hwid,hwid),ylim=c(-lwid,lwid),asp=1,
xlab="mm",ylab="mm")
qr_rectangles(rring,0,mr/2,hring,0)
qr_rectangles(-rring,0,mr/2,hring,0)
segments(-hwid,0,hwid,0)
dev.off()
# 
delta<- flen-sqrt(flen^2-rring^2)
pdf("ws_principal_plane_ladapt.pdf")
par(mfrow=c(4,1))
hwid<-rmax*1.05
lwid<- max(hring)*1.05
plot(rmax*2,0,type="p",xlim=c(-hwid,hwid),ylim=c(-lwid,lwid),asp=1,
xlab="mm",ylab="mm")
qr_rectangles(rring,-delta,mr/2,hring,0)
qr_rectangles(-rring,-delta,mr/2,hring,0)
rr<- c(-hwid,-rring[nring:1]-rmin/2,0,rmin/2,rring,hwid)
dd<- flen-sqrt(flen^2-rr^2)
lines(rr,-dd,type="l")
dev.off()
