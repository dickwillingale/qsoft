FUNCTION QRT_KBS,pcen,pnor,raxi,ipack,rmin,rmax,flen,csize,pitch,wall,plmin,plmax,idf,iq

;          *PCEN    input      centre of telescope aperture
;          *PNOR    input      normal to aperture
;          *RAXI    input      reference axis at centre of aperture
;          *IPACK   input      packing 0 1 module, 1 sunflower, 2 cartesian, 3 wide
;          *                    field cartesian
;          *RAPMIN  input      minimum radius for aperture of constellation
;          *RAPMAX  input      maximum radius for aperture of constellation
;          *FLEN    input      focal length (-ve for 2nd stack)
;          *CSIZE   input      size of each module in constellation
;          *PITCH   input      pitch of K-B slots
;          *WALL    input      wall thickness of K-B slots
;          *PLMIN   input      minimum axial length of K-B slots
;          *PLMAX   input      maximum axial length of K-B slots
;          *IDF     input      deformation index
;          *IQ      input      surface quality index
;          *NMAX    input      maximum number of module coordinate returned
;          *RC      output     radius of each module
;          *PC      output     azimuth of each module
;          *TC      output     rotation of each module
;          *AC      output     axial length of each module
;          *NSET    output     number of module coordinates returned
;          * IF IPACK=0 then single module with centre at x=y=(RAPMIN+RAPMAX)/2/SQRT(2)
;   
; note call_external writes to stdout not the IDL command line unit. 
; the output may be in th eterminal from which you started idlde if you are working in the IDE

rc=DBLARR(2000)
pc=DBLARR(2000)
tc=DBLARR(2000)
ac=DBLARR(2000)
nset=0l

qsoft = GETENV('QSOFT')
  A=CALL_EXTERNAL(qsoft + '/R_libraries/xsrtR.so', "qrt_kbs_", $
      DOUBLE(pcen),$
      DOUBLE(pnor),$
      DOUBLE(raxi),$
      LONG(ipack),$
      DOUBLE(rmin),$
      DOUBLE(rmax),$
      DOUBLE(flen),$
      DOUBLE(csize),$
      DOUBLE(pitch),$
      DOUBLE(wall),$
      DOUBLE(plmin),$
      DOUBLE(plmax),$
      LONG(idf),$
      LONG(iq),$
      LONG(2000),$
      DOUBLE(rc),$
      DOUBLE(pc),$
      DOUBLE(tc),$
      DOUBLE(ac),$
      LONG(nset),$
    /AUTO_GLUE)
    
    kbs = {rc:rc,pc:pc,tc:tc,ac:ac,nset:nset}

RETURN, kbs
END