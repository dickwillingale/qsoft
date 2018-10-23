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
qr_arcs<-function(x,y,r,tl,th,lty=par("lty"),lwd=par("lwd"),col=par("col")) {
        usr<-par()$usr
        del<-min(usr[2]-usr[1],usr[4]-usr[3])/100
        for(i in 1:length(x)) {
		delt<- (th[i]-tl[i])*del/r[i]
                th<-seq(from=tl,to=th,by=delt)
                lines(x[i]+r[i]*cos(th),y[i]+r[i]*sin(th),lty=lty,
                lwd=lwd,col=col)
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
qr_quaderr<-function(x0,y0,x1,y1,y) {
# Quadratic estimator for confidence limit
# x0,y0		parameter and statistic at minimum
# x1,y1		parameter and statistic near minimum
# y		required statistic value
# returns	estimate of parameter corresponding to y
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
qri_peakchisq<-function(fpars) {
	npf<-length(fpars)
        .Fortran("qri_peakchisq",
        as.integer(npf),
        as.double(fpars),
        chisq=double(length=1))$chisq
}
qri_lorentzian<-function(x,p) {
	return(p[1]/(1+((x-p[2])*2.0/p[3])^2))
}
qri_sqbeam<-function(arr,hbeam,blev,bvar) {
	np<-1000
	a<-.Fortran("qri_sqbeam",
	as.integer(dim(arr)[1]),
	as.integer(dim(arr)[2]),
	as.double(arr),
	as.double(hbeam),
	as.double(blev),
	as.double(bvar),
	as.integer(np),
	xpi=double(length=np),
	xpr=double(length=np),
	ypi=double(length=np),
	ypr=double(length=np),
	bug=double(length=np),
	nx=integer(length=1),
	ny=integer(length=1),
	bflux=double(length=1),
	bsigma=double(length=1),
	flux=double(length=1),
	fsigma=double(length=1),
	peak=double(length=2),
	cen=double(length=2),
	rmsx=double(length=1),
	rmsy=double(length=1),
	pi5=double(length=2),
	pi25=double(length=2),
	med=double(length=2),
	pi75=double(length=2),
	pi95=double(length=2),
	hewx=double(length=1),
	hewy=double(length=1),
	w90x=double(length=1),
	w90y=double(length=1))
# Do peak fit
	if(bvar!=0) {
		qri_xchisq<-function(fpars) {
			xm<- qri_lorentzian(a$xpi[1:a$nx],fpars)
			return(sum((a$xpr[1:a$nx]-xm)^2/bvar/a$nx))
		}
		qri_ychisq<-function(fpars) {
			ym<- qri_lorentzian(a$ypi[1:a$ny],fpars)
			return(sum((a$ypr[1:a$ny]-ym)^2/bvar/a$ny))
		}
		delstat<- qchisq(0.9,3)
		derr<- c(F,F,T)
		pval<- max(a$xpr[1:a$nx])
		spars<- c(pval,a$med[1],a$hewx)
		lpars<- c(pval/2,a$med[1]-a$hewx/2,a$hewx/2)
		upars<- c(pval*2,a$med[1]+a$hewx/2,a$hewx*2)
		fx<- qr_srchmin(spars,lpars,upars,qri_xchisq,delstat,derr)
		pval<- max(a$ypr[1:a$ny])
		spars<- c(pval,a$med[2],a$hewy)
		lpars<- c(pval/2,a$med[2]-a$hewy/2,a$hewy/2)
		upars<- c(pval*2,a$med[2]+a$hewy/2,a$hewy*2)
		fy<- qr_srchmin(spars,lpars,upars,qri_ychisq,delstat,derr)
	} else {
		fx<-F
		fy<-F
	}
# Return list of results
	return(list(nx=a$nx,ny=a$ny,
	xpi=a$xpi[1:a$nx],xpr=a$xpr[1:a$nx],ypi=a$ypi[1:a$ny],ypr=a$ypr[1:a$ny],
	bflux=a$bflux,bsigma=a$bsigma,flux=a$flux,
	fsigma=a$fsigma,peak=a$peak,cen=a$cen,
	rmsx=a$rmsx,rmsy=a$rmsy,pi5=a$pi5,pi25=a$pi25,med=a$med,pi75=a$pi75,
	pi95=a$pi95,hewx=a$hewx,hewy=a$hewy,w90x=a$w90x,w90y=a$w90y,
	fitx=fx,fity=fy))
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
# Do peak fit
	if(bvar!=0) {
	  delstat<- qchisq(0.9,4)
	  pval<- arr[a$cen[1],a$cen[2]]
	  spars<- c(pval,a$cen[1],a$cen[2],a$fwhmc/2.)
	  lpars<- c(pval/2,a$cen[1]-a$fwhmc/2,a$cen[2]-a$fwhmc/2,a$fwhmc/2./2)
	  upars<- c(pval*2,a$cen[1]+a$fwhmc/2,a$cen[2]+a$fwhmc/2,a$fwhmc/2.*2)
	  derr<- c(F,T,T,T)
	  f<- qr_srchmin(spars,lpars,upars,qri_peakchisq,delstat,derr)
	} else {
	  f<- F
	}
# The fit parameters are saved in the list fit returned
#	1	peak value (no error range calculated)
#	2	peak X pixel position including 90% upper and lower bounds
#	3	peak Y pixel position including 90% upper and lower bounds
#	4	Gaussian sigma including 90% upper and lower bounds
# construct list of results
	return(list(nsam=a$nsam,bflux=a$bflux,bsigma=a$bsigma,flux=a$flux,
	fsigma=a$fsigma,peak=a$peak,cen=a$cen,tha=a$tha,
	rmsa=a$rmsa,rmsb=a$rmsb,
	fwhm=a$fwhm,hew=a$hew,w90=a$w90,
	fwhmp=a$fwhmp,hewp=a$hewp,w90p=a$w90p,
	fwhmc=a$fwhmc,hewc=a$hewc,w90c=a$w90c,
	fit=f))
}
qri_lecbeam<-function(arr,s,h,blev,bvar,nt) {
	a<-.Fortran("qri_lecbeam",
	as.integer(dim(arr)[1]),
	as.integer(dim(arr)[2]),
	as.double(arr),
	as.double(s),
	as.double(h),
	as.double(blev),
	as.double(bvar),
	as.integer(nt),
	qua=double(length=nt*nt),
	nqua=double(length=nt*nt),
	nsam=integer(length=1),
	bflux=double(length=1),
	bsigma=double(length=1),
	flux=double(length=1),
	fsigma=double(length=1),
	peak=double(length=2),
	cen=double(length=2),
	hew=double(length=1),
	w90=double(length=1),
	ahew=double(length=1),
	aw90=double(length=1),
	fpeak=double(length=1))
# functions for fitting the lobster eye quadrant
        quafun<- function(xp,yp,par) {
	 #par[1] peak value
         #par[2] G width of Lorentzian central spot
         #par[3] eta strength of cross-arm wrt central spot
         eg<- par[3]*par[2]/hqua
	 x1<-1/(1+(xp*2.0/par[2])^2)
	 x2<-eg*(1-(xp/hqua)^2)
         y1<-1/(1+(yp*2.0/par[2])^2)
	 y2<-eg*(1-(yp/hqua)^2)
         return((x1*y1+x1*y2+x2*y1+x2*y2)*par[1]/(1+eg)^2)
        }
        quastat<- function(p) {
	   mm<-outer(xqua,yqua,quafun,p)
           return(sum((zqua-mm)^2/mm*(nqua>0)))
        }
# Do fit to quadrant profile
        delstat<- qchisq(0.9,2)
        zqua<- a$qua
	nqua<- a$nqua
        xqua<- (1:nt)-0.5
        yqua<- xqua
        hqua<- h
        spars<- c(a$qua[1],a$hew,1.0)
        lpars<- c(a$qua[1]/10,a$hew/10,0.1)
        upars<- c(a$qua[1]*10,a$hew*10,10.0)
        derr<- c(F,F,F)
        f<- qr_srchmin(spars,lpars,upars,quastat,delstat,derr)
        mod<- outer(xqua,yqua,quafun,f$par)
	dim(mod)<- c(nt,nt)
# construct list of results
	dim(a$qua)<-c(nt,nt)
	return(list(qua=a$qua,
	nsam=a$nsam,bflux=a$bflux,bsigma=a$bsigma,flux=a$flux,
	fsigma=a$fsigma,peak=a$peak,cen=a$cen,
	hew=a$hew,w90=a$w90,ahew=a$ahew,aw90=a$aw90,fpeak=a$fpeak,
	norm=f$par[1],G=f$par[2],eta=f$par[3],xqua=xqua,yqua=yqua,mod=mod))
}
qri_lecimage<-function(s,h,b,xcen,ycen,nx,ny) {
	b<-.Fortran("qri_lecimage",
	as.double(s),
	as.double(h),
	as.double(b),
	as.double(xcen),
	as.double(ycen),
	as.integer(nx),
	as.integer(ny),
	arr=double(length=nx*ny))
	dim(b$arr)<-c(nx,ny)
	invisible(b$arr)
}
qri_lepsf<-function(s,h,g,eta,xcen,ycen,nx,ny) {
	b<-.Fortran("qri_lepsf",
	as.double(s),
	as.double(h),
	as.double(g),
	as.double(eta),
	as.double(xcen),
	as.double(ycen),
	as.integer(nx),
	as.integer(ny),
	arr=double(length=nx*ny))
	dim(b$arr)<-c(nx,ny)
	invisible(b$arr)
}
qri_lebin<-function(xe,ye,s,h,g,eta,nx,ny) {
	ne<-length(xe)
	b<-.Fortran("qri_lebin",
	as.integer(ne),
	as.double(xe),
	as.double(ye),
	as.double(s),
	as.double(h),
	as.double(g),
	as.double(eta),
	as.integer(nx),
	as.integer(ny),
	arr=double(length=nx*ny))
	dim(b$arr)<-c(nx,ny)
	invisible(b$arr)
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
