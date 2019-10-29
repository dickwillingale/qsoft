*+QRA_GRBSETFIT set parameters for start of fit
        SUBROUTINE QRA_GRBSETFIT(NPF,ITOFIT,NTBAT,FBAT,NTXRT,FXRT,SPARS,
     +        LPARS,UPARS,BGOOD,XGOOD,NDOF)
        IMPLICIT NONE
        INTEGER NPF,ITOFIT(NPF),NTBAT,NTXRT
        INTEGER FBAT(NTBAT),FXRT(NTXRT),BGOOD,XGOOD,NDOF
        DOUBLE PRECISION SPARS(NPF),LPARS(NPF),UPARS(NPF)
*NPF        input        number of parameters to be fitted
*ITOFIT        input        indices of parameters to be fitted
*NTBAT        input        number of BAT time samples
*FBAT        input        array BAT flags 1 include, 0 exclude
*NTXRT        input        number of XRT time samples
*FXRT        input        array XRT flags 1 include, 0 exclude
*SPARS        output        starting value of parameters to be fitted
*LPARS        output        lower limit of parameters to be fitted
*UPARS        output        upper limit of parameters to be fitted
*BGOOD        output        number of good BAT time samples
*XGOOD        output        number of good XRT time samples
*NDOF        output        number of degrees of freedom for the fit
*-Author Dick Willingale 2013-Jan-24
        include 'QR_COM'
        include 'GRBS_COM'
        integer j,jj,npul,jpul,jpar
C
        IF(ISTAT.NE.0) RETURN
C
        IF(NTBAT.ne.m1) THEN
                write(*,*) 'QRA_GRBSETFIT - wrong number of BAT flags',
     +                NTBAT,m1
                ISTAT=1
                RETURN
        ENDIF
        IF(NTXRT.ne.m2) THEN
                write(*,*) 'QRA_GRBSETFIT - wrong number of XRT flags',
     +                NTXRT,m2
                ISTAT=1
                RETURN
        ENDIF
C Get fit parameters values from common
        npul=pars(1)
        nfit=NPF
        DO J=1,NPF
                jj=ITOFIT(J)
                if(jj.gt.npars.or.jj.lt.2) then
                     write(*,*)
     +               'QRA_GRBSETFIT - parameter index out of range',jj
                            return
                endif
                ifit(J)=jj
                SPARS(J)=pars(jj)
C Set allowed range for this parameter
                jpul=int((jj-2)/9)+1
                if(jpul.le.npul) then
C prompt pulse
                        jpar=mod(jj-2,9)+1
                        if(jpar.eq.1) then
C b1
                                LPARS(J)=-2.5
                                UPARS(J)=4.0
                        elseif(jpar.eq.2) then
C amb
                                LPARS(J)=0.0
                                UPARS(J)=10.0
                        elseif(jpar.eq.3) then
C d
                                LPARS(J)=-10.0
                                UPARS(J)=0.0
                        elseif(jpar.eq.4) then
C a
                                LPARS(J)=0.0
                                UPARS(J)=10.0
                        elseif(jpar.eq.5) then
C fr 
                                LPARS(J)=0.01
                                UPARS(J)=0.99
                        elseif(jpar.eq.6) then
C Tf 
                                LPARS(J)=1.e-2
                                UPARS(J)=1.e5
                        elseif(jpar.eq.7) then
C Epk
                                LPARS(J)=0.3
                                UPARS(J)=1.e5
                        elseif(jpar.eq.8) then
C Tpk
                                LPARS(J)=-x1(1)
                                UPARS(J)=1.e7
                        elseif(jpar.eq.9) then
C Fpk
                                LPARS(J)=0.0
                                UPARS(J)=1.e7
                        endif
                else
C Get afterglow type
                        if(pars(npul*9+2).eq.100) then
                                jpar=jj-npul*9-3
                                if(jpar.eq.0) then
C vv
                                        LPARS(J)=-2.0
                                        UPARS(J)=3.0
                                endif
                        elseif(pars(npul*9+2).eq.200) then
                                jpar=jj-npul*9-10
                                if(jpar.eq.-7) then
C od - optical depth of dust at 1 keV a=0.1 microns
                                        LPARS(J)=0.0
                                        UPARS(J)=100.0
                                elseif(jpar.eq.-6) then
C am - minimum grain size microns
                                        LPARS(J)=0.0001
                                        UPARS(J)=0.1
                                elseif(jpar.eq.-5) then
C ap - maximum grain size microns
                                        LPARS(J)=0.1
                                        UPARS(J)=3.0
                                elseif(jpar.eq.-4) then
C dd - grain size distribution index
                                        LPARS(J)=1.0
                                        UPARS(J)=6.0
                                elseif(jpar.eq.-3) then
C zd - redshift for dust
                                        LPARS(J)=0
                                        UPARS(J)=0
                                elseif(jpar.eq.-2) then
C rl - minimum dust distance pc
                                        LPARS(J)=1
                                        UPARS(J)=10000
                                elseif(jpar.eq.-1) then
C rd - distance range of dust layers pc
                                        LPARS(J)=1
                                        UPARS(J)=10000
                                elseif(jpar.eq.0) then
C rn  - number of dust layers starting at rl
                                        LPARS(J)=1
                                        UPARS(J)=20
                                endif
                        else
                                jpar=jj-npul*9-1
                        endif
                        if(jpar.eq.1) then
C ba
                                LPARS(J)=-3.5
                                UPARS(J)=1.0
                        elseif(jpar.eq.2) then
C Ea
                                LPARS(J)=0.3
                                UPARS(J)=1000.0
                        elseif(jpar.eq.3) then
C Tr
                                LPARS(J)=1.0
                                UPARS(J)=1.e4
                        elseif(jpar.eq.4) then
C Ta
                                LPARS(J)=10.0
                                UPARS(J)=1.e6
                        elseif(jpar.eq.5) then
C Fm
                                LPARS(J)=0.0
                                UPARS(J)=1.e6
                        elseif(jpar.eq.6) then
C aa
                                LPARS(J)=0.0
                                UPARS(J)=10.0
                        elseif(jpar.eq.7) then
C Tb
                                LPARS(J)=10.0
                                UPARS(J)=1.e8
                        elseif(jpar.eq.8) then
C ab
                                LPARS(J)=-3.0
                                UPARS(J)=3.0
                        elseif(jpar.eq.9) then
C Tb1
                                LPARS(J)=10.0
                                UPARS(J)=1.e8
                        elseif(jpar.eq.10) then
C ab1
                                LPARS(J)=-3.0
                                UPARS(J)=3.0
                        elseif(jpar.eq.11) then
C Tb2
                                LPARS(J)=10.0
                                UPARS(J)=1.e8
                        elseif(jpar.eq.12) then
C ab2
                                LPARS(J)=-3.0
                                UPARS(J)=3.0
                        endif
                endif
        ENDDO
C
        ncall=0
C Set BAT flags
        BGOOD=0
        DO J=1,m1
                ix1(j)=FBAT(J)
                if(FBAT(J).ne.0) then
                        BGOOD=BGOOD+1
                endif
        ENDDO
C Set XRT flags
        XGOOD=0
        DO J=1,m2
                ix2(j)=FXRT(J)
                if(FXRT(J).ne.0) then
                        XGOOD=XGOOD+1
                endif
        ENDDO
C Set NDOF
        NDOF=BGOOD*mc1+XGOOD*mc2-NPF
        END
