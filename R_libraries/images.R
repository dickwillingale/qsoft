# qR routines for images
if(!exists("qri_init")) {
qri_init<-function() {
	.Fortran("qri_init")
	invisible()
}
qri_init()
# utilities
qr_arrows<-function(x,y,da,t) {
        x1<-x+cos(t)*da
        y1<-y+sin(t)*da
        for(i in 1:length(x))
                lines(c(x[i],x1[i]),c(y[i],y1[i]),type='l')
        invisible()
}
qr_rectangles<-function(x,y,w,h,t,lty=par("lty"),fg=par("fg")) {
        cth<-cos(t)
        sth<-sin(t)
        xa<-w*cth
        xb<-w*sth
        ya<-h*cth
        yb<-h*sth
        x1<-x+xa+yb
        x2<-x+xa-yb
        x3<-x-xa-yb
        x4<-x-xa+yb
        y1<-y+xb-ya
        y2<-y+xb+ya
        y3<-y-xb+ya
        y4<-y-xb-ya
        for(i in 1:length(x))
                polygon(c(x1[i],x2[i],x3[i],x4[i]),c(y1[i],y2[i],y3[i],y4[i]),
                lty=lty,fg=fg)
        invisible()
}
qr_circles<-function(x,y,r,lty=par("lty"),lwd=par("lwd"),col=par("col")) {
        usr<-par()$usr
        del<-min(usr[2]-usr[1],usr[4]-usr[3])/100
        for(i in 1:length(x)) {
                nth<-max(8,trunc((2*pi*r[i]/del)))
                th<-seq(from=-pi,to=pi,length=nth)
                polygon(x[i]+r[i]*cos(th),y[i]+r[i]*sin(th),lty=lty,
                lwd=lwd,border=col)
        }
        invisible()
}
qr_collut<-function(fname) {
        t<-read.table(fname,header=T)
        rgb(t$r,t$g,t$b,1)
}
qr_xyztrans<-function(x,y,z,cen,aer) {
        a<- .Fortran("qr_xyztrans",
        as.integer(length(x)),
        as.double(x),
        as.double(y),
        as.double(z),
        as.double(cen),
        as.double(aer),
        xp=double(length(x)),
        yp=double(length(y)),
        zp=double(length(z)))
        return(list(xp=a$xp,yp=a$yp,zp=a$zp))
}
qr_rawread<-function(filename,scale,linearize,linfile) {
# Read .raw MCP event files
# Adrian Martindale's version 21st July 2015
	bfile<- file(filename,"rb")
	header<- readBin(bfile,raw(),n=512,size=1)
	header<- paste(rawToChar(header,multiple=T),sep="")
	nbuf<- 4000000
	idat<- integer()
	ndat<- 0
	end<- "little"
	a<- readBin(bfile,integer(), n=nbuf,size=2,endian=end)
	ni<- length(a)
	while(ni>0) {
		idat<- c(idat,a)
		ndat<- ndat+ni
		a<- readBin(bfile,integer(), n=nbuf,size=2,endian=end)
		ni<- length(a)
	}
	close(bfile)
	ne<- ndat/4
	dim(idat)<- c(4,ne)
#Dick's algorithm
#	x<- idat[3,1:ne]/(idat[1,1:ne]+idat[3,1:ne])-0.5
#	y<- idat[2,1:ne]/(idat[2,1:ne]+idat[4,1:ne])-0.5
#Image_Display algorithm
        x <- idat[2,1:ne]/(idat[2,1:ne]+idat[4,1:ne])-0.5
        y <- idat[3,1:ne]/(idat[1,1:ne]+idat[3,1:ne])-0.5  
# linearise the data?
        cat("linearise = ",linearize,"\n")
	if(linearize == 1){
	  print('Linearising data')
	  #format x and y as done by IDL (no shift and times number of pixels)
          x <- (x+0.5)*512.0
          y <- (y+0.5)*512.0
#read linearisation data file
	  lindata <- read.table(linfile,header=T)
	 
	  #format xfit and yfit as 4x4 matrices. NB R is reading as if it is row major 
	  #cf IDL which is column major no issue as the indices work the same way.
	  xfit <- lindata$xfit
          yfit <- lindata$yfit
          dim(xfit) <- c(4,4)
	  dim(yfit) <- c(4,4)
    
#reset scale from unlinearised estimate to real number
          mmperpix <- lindata$pitch[1]/lindata$scale[1]
          scal <- mmperpix*512.0/2.0
#NB 512/2 is half the width of the plot in image_display
          cat("scale set to = ",scal, " = ", mmperpix,
          "mm per pixel for an image_display pixel","\n")
	  nx = rep(0.0,ne)
	  ny = rep(0.0,ne)
	  for(k in 0:3){
	    for(j in 0:3){
	      tmp <- x^k * y^j
	      nx <- nx + xfit[j+1,k+1]*tmp
	      ny <- ny + yfit[j+1,k+1]*tmp
	    }
	  }
	  x <- (nx/512.0)-0.5
	  y <- (ny/512.0)-0.5
	} else {
          scal=scale
        }
	cs<- sqrt(2)
	ss<- cs
# need x and y to be the same for R and IDL. NB becaus of row vs column major the x-axis needs to be flipped 
# to maintain consistency with IDL definitions which are implicit in the linearisation file.
	xp<- (y*ss-x*cs)*scal
	yp<- (x*ss+y*cs)*scal
	ph<- (idat[1,1:ne]+idat[2,1:ne]+idat[3,1:ne]+idat[4,1:ne])
	events<- data.frame(xp=xp,yp=yp,ph=ph)
	return(list(header=header,events=events))
}
# Image processing routines - prefix qri_
qri_setfield<-function(nx,xleft,xright,ny,ybot,ytop) {
	.Fortran("qri_setfield",
	as.integer(nx),
	as.double(xleft),
	as.double(xright),
	as.integer(ny),
	as.double(ybot),
	as.double(ytop))
	xsam<-(xright-xleft)/nx
	ysam<-(ytop-ybot)/ny
	invisible(list(nx=nx,xleft=xleft,xright=xright,xsam=xsam,
	ny=ny,ybot=ybot,ytop=ytop,ysam=ysam,
	xp=seq(from=xleft+xsam/2,by=xsam,length.out=nx),
	yp=seq(from=ybot+ysam/2,by=ysam,length.out=ny)))
}
qri_setsky<-function(xtodeg,ytodeg,ipr,mjd,ra,dec,roll)
	.Fortran("qri_setsky",
	as.double(xtodeg),
	as.double(ytodeg),
	as.integer(ipr),
	as.double(mjd),
	as.double(ra),
	as.double(dec),
	as.double(roll))
qri_setpos<-function(ipos,p) {
	.Fortran("qri_setpos",
	as.integer(ipos),
	as.double(p))
	invisible()
}
qri_getpos<-function()
	a<-.Fortran("qri_getpos",
	pix=double(length=2),
	xyl=double(length=2),
	aes=double(length=2),
	equ=double(length=2),
	ecl=double(length=2),
	gal=double(length=2))
qri_printpos<-function() {
	a<-qri_getpos()
	print(
	c(sprintf("pixel       %10.5f %10.5f",a$pix[1],a$pix[2]),
	sprintf("xyloc       %10.5f %10.5f",a$xyl[1],a$xyl[2]),
	sprintf("sploc       %10.5f %10.5f",a$aes[1],a$aes[2]),
	sprintf("Equatorial  %10.5f %10.5f",a$equ[1],a$equ[2]),
	sprintf("Ecliptic    %10.5f %10.5f",a$ecl[1],a$ecl[2]),
	sprintf("Galactic    %10.5f %10.5f",a$gal[1],a$gal[2])))
}
qri_makeimage<-function(a,xleft,xright,ybot,ytop) {
	nx<-dim(a)[1]; ny<-dim(a)[2]
	b<-qri_setfield(nx,xleft,xright,ny,ybot,ytop)
	invisible(list(data_array=a,
	nx=nx,xleft=xleft,xright=xright,xsam=b$xsam,
	ny=ny,ybot=ybot,ytop=ytop,ysam=b$ysam,
	xp=b$xp,yp=b$yp,xlim=c(xleft,xright),ylim=c(ybot,ytop),
	zlim=range(a)))
}
qri_binxy<-function(x,y,iq,w,xleft,xright,nx,ybot,ytop,ny) {
	a<-.Fortran("qri_binxy",
	as.integer(length(x)),
	as.double(x),
	as.double(y),
	as.integer(length(iq)),
	as.integer(iq),
	as.integer(length(w)),
	as.double(w),
	as.double(xleft),
	as.double(xright),
	as.double(ybot),
	as.double(ytop),
	as.integer(nx),
	as.integer(ny),
	data_array=double(length=nx*ny))
	dim(a$data_array)<-c(nx,ny)
	invisible(qri_makeimage(a$data_array,xleft,xright,ybot,ytop))
}
qri_blur<- function(a,m) {
	nx<-dim(a)[1]; ny<-dim(a)[2]
	mx<-dim(m)[1]; my<-dim(m)[2]
	b<-.Fortran("qri_blur",
	as.integer(nx),
	as.integer(ny),
	as.double(a),
	as.integer(mx),
	as.integer(my),
	as.double(m),
	arr=double(length=nx*ny))
	dim(b$arr)<-c(nx,ny)
	invisible(b$arr)
}
qri_displayimage<- function(ima) {
	qri_zoomimage(ima,ima$xlim,ima$ylim,ima$zlim)
	invisible()
}
qri_zoomimage<-function(ima,xl,yl,zl) {
	if(is.null(ima$xlab)) xla<-"" else xla<-ima$xlab
	if(is.null(ima$ylab)) yla<-"" else yla<-ima$ylab
	plot.new()
	plot.window(xl,yl,asp=1,xaxs="i",yaxs="i")
#,xaxs="i",yaxs="i")
	#plot.window(xl,yl)
	clut<-qr_collut(paste(Sys.getenv("QSOFT"),"/data/lut6.dat",sep=""))
	image(ima$xp,ima$yp,ima$data_array,zl,add=TRUE,col=clut,
	useRaster=TRUE,asp=1)
	title(xlab=xla,ylab=yla)
	box(which="outer")
	axis(1,pos=yl[1]);axis(2,pos=xl[1])
	invisible()
}
qri_displayrays<-function(a,nmax,az,el,rl,lims=NA) {
	aer<-c(az,el,rl)*pi/180.
	b<-qr_xyztrans(a$RXP,a$RYP,a$RZP,c(0,0,0),aer)
	if(is.na(lims[1])) {
		xmin<-min(b$xp)
		xmax<-max(b$xp)
		ymin<-min(b$yp)
		ymax<-max(b$yp)
	} else {
		xmin<- lims[1]
		xmax<- lims[2]
		ymin<- lims[3]
		ymax<- lims[4]
	}
	plot.new()
        plot.window(c(xmin,xmax),c(ymin,ymax),asp=1)
	np<-length(a$IQU)
	nplot<-0
	is<-1
	i<-1
	while((nplot<nmax)&&(i<=np)) {
		ii<-i+1	
		if(a$IQU[ii]==-2 || i==np) {
			lines(b$xp[is:i],b$yp[is:i],type="l")
			nplot<-nplot+1
			is<-ii
		}
		i<-ii
	}
	title(xlab="x mm",ylab="y mm",line=2)
	axis(1);axis(2)
	invisible()
}
qri_beam<-function(arr,rbeam,blev,bvar) {
	a<-.Fortran("qri_beam",
	as.integer(dim(arr)[1]),
	as.integer(dim(arr)[2]),
	as.double(arr),
	as.double(rbeam),
	as.double(blev),
	as.double(bvar),
	nsam=integer(length=1),
	bflux=double(length=1),
	bsigma=double(length=1),
	flux=double(length=1),
	fsigma=double(length=1),
	peak=double(length=2),
	cen=double(length=2),
	tha=double(length=1),
	rmsa=double(length=1),
	rmsb=double(length=1),
	fwhm=double(length=1),
	hew=double(length=1),
	w90=double(length=1),
	fwhmp=double(length=1),
	hewp=double(length=1),
	w90p=double(length=1),
	fwhmc=double(length=1),
	hewc=double(length=1),
	w90c=double(length=1))
# construct list of results
	return(list(nsam=a$nsam,bflux=a$bflux,bsigma=a$bsigma,flux=a$flux,
	fsigma=a$fsigma,peak=a$peak,cen=a$cen,tha=a$tha,
	rmsa=a$rmsa,rmsb=a$rmsb,
	fwhm=a$fwhm,hew=a$hew,w90=a$w90,
	fwhmp=a$fwhmp,hewp=a$hewp,w90p=a$w90p,
	fwhmc=a$fwhmc,hewc=a$hewc,w90c=a$w90c))
}
qri_annulus<-function(arr,rmin,rmax) {
	a<-.Fortran("qri_annulus",
	as.integer(dim(arr)[1]),
	as.integer(dim(arr)[2]),
	as.double(arr),
	as.double(rmin),
	as.double(rmax),
	nsam=integer(length=1),
	blev=double(length=1),
	bvar=double(length=1))
	return(list(nsam=a$nsam,blev=a$blev,bvar=a$bvar))
}
qri_abeam<-function(arr,rbeam,counting) {
	rmin=rbeam
	rmax=rbeam*2
	a<-.Fortran("qri_annulus",
	as.integer(dim(arr)[1]),
	as.integer(dim(arr)[2]),
	as.double(arr),
	as.double(rmin),
	as.double(rmax),
	nsam=integer(length=1),
	blev=double(length=1),
	bvar=double(length=1))
#
	if(counting) var<--1 else var<-a$bvar
#
	b<-.Fortran("qri_beam",
	as.integer(dim(arr)[1]),
	as.integer(dim(arr)[2]),
	as.double(arr),
	as.double(rbeam),
	as.double(a$blev),
	as.double(var),
	nsam=integer(length=1),
	bflux=double(length=1),
	bsigma=double(length=1),
	flux=double(length=1),
	fsigma=double(length=1),
	peak=double(length=2),
	cen=double(length=2),
	tha=double(length=1),
	rmsa=double(length=1),
	rmsb=double(length=1),
	fwhm=double(length=1),
	hew=double(length=1),
	w90=double(length=1))
	return(b)
}
}
