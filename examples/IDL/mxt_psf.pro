PRO mxt_PSF

  model= 'stm'
  ekev= 1.49
  psca= 3
  ; useful vectors
  sn=[1,0,0]
  nn=[0,0,0]
  rx=[0,1,0]
  ; get plate data
  RESTORE,'/Users/Adrian/IDLWorkspace/am136/qIDL/examples/'+model+"_mpos.sav"
  ;PRINT,plt
  idef= 1
  ;
  flen= rocf/2.
  ; set detector twice size of MXT detector
  dpix= 0.075
  np= 512*3
  hdet= dpix*np/2.
  dlim=[-hdet,-hdet,hdet,hdet]
  hmxt= hdet/2
  ; surface optical properties
  xop=QRT_XOPT(mspec,rho,ekev,1)
  is=1
  ; set up source
  ccy=SIN(0.0*!dpi/180)
  ccz=SIN(0*!dpi/180)
  di=[-SQRT(1.0-ccy^2-ccz^2),ccy,ccz]
  slim=[-hwid,-hwid,hwid,hwid]
  nray=1000000
  psrc=[flen*1.05,0,0]
  ; set position of centre of square pore MCP array
  pmcp=[flen,0,0]
  ; rotate detector axes wrt the plate reference axes
  thetadet=-0.*!dpi/180.
  rxdet=[0,cos(thetadet),sin(thetadet)]
  ;
  qrt_q=QRT_RESET()
  qrt_q=QRT_SURFACE(is,1,ekev,nrough,brk,rind,xop.alpha,xop.gamma,0,0,0,0,0)
  qrt_q=QRT_SOURCE(2,di,nn,psrc,sn,rx,slim,0,nray,0)
  qrt_q=QRT_DEFORM(idef,1,5,nd,1)
  qrt_q=QRT_DEFMAT(idef,1,xd,yd,zx)
  qrt_q=QRT_DEFMAT(idef,2,xd,yd,zy)
  qrt_q=QRT_DEFMAT(idef,3,xd,yd,za)
  qrt_q=QRT_DEFMAT(idef,4,xd,yd,zd)
  qrt_q=QRT_DEFMAT(idef,5,xd,yd,fd)
  qrt_q=QRT_SQMPOARR(pmcp,sn,rx,rocf,hwid,idef,plt)
  qrt_q=QRT_DETECTOR(2,nn,sn,rxdet,dlim,0)
  ;qrt_list()
  results=QRT_TRACE(0,2,-1)
  ; get detected positions
  detpos=READ_QRT('detected.dat',/det)
  ; create image
  a=QRI_BINXY(detpos.ydet,detpos.zdet,0,detpos.area,-hdet,hdet,np,-hdet,hdet,np)
  
  ;pdf(paste(model,"_psf_",ekev,".pdf",sep=""))
  a.xlab="mm"
  a.ylab="mm"
  QRI_ZOOMIMAGE,a,[-hdet,hdet],[-hdet,hdet],[0,a.zlim[1]/psca]
  
  ;lines(c(-hmxt,hmxt,hmxt,-hmxt,-hmxt),c(-hmxt,-hmxt,hmxt,hmxt,-hmxt),col="red")
  ;title(main=paste(model,"Ray Tracing Simulation",ekev,"keV"))
  ;dev.off()
  QRI_s=QRI_SETFIELD(np,-hdet,hdet,np,-hdet,hdet)
  QRI_sp=QRI_SETPOS(2,[0,0])
  ; beam size of plate
  rpix=hdet/dpix
  b=QRI_BEAM(a.data_array,rpix,0,0)
  ;fwhm=b.fit.par[4]*2.0*dpix
  ;amfwhm=(fwhm/flen)*180*60/pi
  ;cat("Lorentzian fit FWHM mm arcmins",fwhm,amfwhm,"\n")
  ; PSF Lobster eye cross-beam analysis - assume L/d=60
  ;sb= hdet/dpix*2
  ;hb= 2.0/60*flen/dpix
  ;nqua= trunc(hb)
  ;blecp=qri_lecbeam(a$data_array,sb,hb,0,0,nqua)
  ;alGp= blecp$G*dpix
  ;lhewp= blecp$hew*dpix
  ;amalGp=(alGp/flen)*180*60/pi
  ;amlhewp=(lhewp/flen)*180*60/pi
  ;albeamp=blecp$flux/blecp$fpeak*dpix*dpix
  ;amlbeamp= sqrt(albeamp)/flen*180*60/pi
  ;cat("PSF Lobster eye beam pixels sb hb hew",sb,hb,blecp$hew,"\n")
  ;cat("PSF area cm2 in lobster eye patch",blecp$flux*0.01,"\n")
  ;cat("PSF fitted norm G eta",blecp$norm,blecp$G,blecp$eta,"\n")
  ;cat("PSF fitted G and HEW mm",alGp,lhewp,"\n")
  ;cat("PSF fitted G and HEW arcmins",amalGp,amlhewp,"\n")
  ;cat("PSF beam area mm2",albeamp,"\n")
  ;cat("PSF sqrt(beam area) arcmins",amlbeamp,"\n")
  ;; Check total area in detector array
  ;totarea= sum(a$data_array)*0.01
  ;cat("total detected area cm2",totarea,"\n")
  ;
  ;dfile= paste(model,"_psf_",ekev,".dat",sep="")
  ;write(a$data_array,file=dfile,ncolumns=1)
  ; Set up analytical model
  ;mod= qri_lepsf(sb,hb,blecp$G,blecp$eta,np/2,np/2,np,np)
  ;md= qri_makeimage(mod,-hdet,hdet,-hdet,hdet)
  ;qri_zoomimage(md,c(-hdet,hdet),c(-hdet,hdet),c(0,md$zlim[2]/psca))
  ;lines(c(-hmxt,hmxt,hmxt,-hmxt,-hmxt),c(-hmxt,-hmxt,hmxt,hmxt,-hmxt),col="red")
  ;title(main=paste(model,"Analytical Model",ekev,"keV"))
  ;
  ;bb=qri_beam(md$data_array,rpix,0,-1)
  ;mfwhm=bb$fit$par[4]*2.0*dpix
  ;mamfwhm=(mfwhm/flen)*180*60/pi
  ;cat("Mod Lorentzian fit FWHM mm arcmins",mfwhm,mamfwhm,"\n")
  ; Model Lobster eye cross-beam analysis - assume L/d=60
  ;mlecp=qri_lecbeam(md$data_array,sb,hb,0,0,nqua)
  ;malGp= mlecp$G*dpix
  ;mlhewp= mlecp$hew*dpix
  ;mamalGp=(malGp/flen)*180*60/pi
  ;mamlhewp=(mlhewp/flen)*180*60/pi
  ;malbeamp=mlecp$flux/mlecp$fpeak*dpix*dpix
  ;mamlbeamp= sqrt(malbeamp)/flen*180*60/pi
  ;cat("Mod Lobster eye beam pixels sb hb hew",sb,hb,mlecp$hew,"\n")
  ;cat("Mod area cm2 in lobster eye patch",mlecp$flux*0.01,"\n")
  ;cat("Mod fitted norm G eta",mlecp$norm,mlecp$G,mlecp$eta,"\n")
  ;cat("Mod fitted G and HEW mm",malGp,mlhewp,"\n")
  ;cat("Mod fitted G and HEW arcmins",mamalGp,mamlhewp,"\n")
  ;cat("Mod beam area mm2",malbeamp,"\n")
  ;cat("Mod sqrt(beam area) arcmins",mamlbeamp,"\n")


RETURN
END