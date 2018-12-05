  ;+
  ;
  ; NAME:
  ;   QR_SRCHMIN
  ;
  ; PURPOSE:
  ;   This Library implememnts th QSOFT analysis package in IDL using calls to sharable objects
  ;   external to IDL. Original routines were in general written in Fortan 77 and are compiled
  ;   into sharable objects.
  ;
  ;   The purpose of this routine is to search for minimum statistic and return best fit parameters and
  ;   confidence limits of the parameters
  ;
  ;
  ; :Categories:
  ;    Image analysis
  ;
  ; :Keywords:
  ;
  ; :Examples:
  ;   A=QR_SRCHMIN(pars,pl,ph,stat,delstat,derr)
  ;     INPUTS:
  ;       # pars    initial parameter values
  ;       # pl    hard lower limit of parameter values
  ;       # ph    hard upper limit of parameter values
  ;       # stat    the statistic function to be minimised
  ;       # delstat the change in statistic for confidence limits
  ;       # derr    logical array to specify parameters for confidence estimates
  ;       # returns list from optim() plus confidence limits parlo and parhi
  ;
  ; :Author:
  ;       Adrian Martindale based on R. Willingale's code in R
  ;
  ; :History:
  ;     Change History::
  ;       Written by: Adrian Martindale, 15/10/2018 to implement QSOFT in IDL.
  ;
  ;-

FUNCTION QR_SRCHMIN, pars,pl,ph,stat,delstat,derr

  np = N_ELEMENTS(pars)
; make local copies of hard limits
  prl=pl
  prh=ph
  
; find statistic minimum and estimate hessian matrix
  ;ft=optim(pars,stat,gr=NULL,method="L-BFGS-B",lower=prl,upper=prh,hessian=T)
  
  P = REPLICATE({fixed:0, limited:[1,1], limits:[0.D,0.D]},4)
  P.limits(0,*) = pl
  P.limits(1,*) = ph
  
  ft = TNMIN(stat,PARINFO=P,/autoderivative)

;NEED TO CONVERT ft to a structure
;how to get hessian???


; estimate errors on parameters
  ft.delstat = delstat
  dig = ABS(DIAG_MATRIX(ft.hessian))
  delpar = DBLARR(N_ELEMENTS(np))
; fix hard limits of parameters for which we don't want error estimate
  FOR k = 0, np-1 DO BEGIN
  ; Trap zero on hessian diagonal
    IF dig(k) EQ 0 THEN BEGIN
      delpar(k) = 0
      IF ~derr(k) THEN BEGIN
        prl(k) = ft.par(k)-ft.par(k)/200
        prh(k) = ft.par(k)+ft.par(k)/200 
      ENDIF
    ENDIF ELSE BEGIN  
      delpar(k) = sqrt(2*delstat/dig(k))
      IF ~!derr(k) THEN BEGIN
        prl(k) = ft.par(k)-delpar(k)/200
        prh(k) = ft.par(k)+delpar(k)/200
      ENDIF
    ENDELSE
  ENDFOR 

PRINT, "Min statistic",ft.value,"\n"
PRINT, "best fit parameters",ft.par,"\n"
PRINT, "diag",dig,"\n"
PRINT, "delpar",delpar,"\n"
nfr = TOTAL(derr)
PRINT, "Find limits for ",nfr," parameters\n"
parhi = DBLARR(N_ELEMENTS(np))
parlo = DBLARR(N_ELEMENTS(np))

FOR k = 0,np-1 DO BEGIN
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
ENDFOR

ft$parlo<-parlo
ft$parhi<-parhi
return(ft)
}
RETURN
END