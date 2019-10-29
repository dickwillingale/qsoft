*+QRA_GRBCHISQ calculate Chi-squared for given parameters
        SUBROUTINE QRA_GRBCHISQ(NPF,FPARS,CHISQ)
        IMPLICIT NONE
        INTEGER NPF
        DOUBLE PRECISION FPARS(NPF),CHISQ
*NPF        input        number of parameters to be fitted
*FPARS        input        fitting parameter values
*CHISQ        output        Chi-squared value
*-Author Dick Willingale 2013-Jan-24
        include 'QR_COM'
        include 'GRBS_COM'
C
        integer nr
        PARAMETER (nr=8)
        double precision rates(nr),delc,ylog,rlog,elog
        double precision flux,elo,ehi,tim
        integer j,jj,i,ii,k,npul
C
        IF(ISTAT.NE.0) RETURN
C
        if(NPF.ne.nfit) then
                write(*,*)
     +    'QRA_GRBCHISQ - incorrect number of parameters',NPF,nfit
                ISTAT=1
                return
        endif
C 
        ncall=ncall+1
C Poke input parameter values into common
        DO J=1,NPF
                JJ=ifit(J)
                pars(JJ)=FPARS(J)
C        write(*,*) 'fn',jj,pars(jj)
        ENDDO

        chisq=0.0
C loop for BAT data points
        elo=15.0
        ehi=350.0
        do j=1,m1
            if(ix1(j).ne.0) then
                call qra_grbrates(x1(j),npars,pars,-1,
     +          elo,ehi,nr,rates,flux)
                ncall=ncall+1
                do i=1,mc1
                        if(il.eq.0.or.il.eq.1) then
                                if(ye1(i,j).ne.0.0) then
                           delc=(y1(i,j)-rates(i+nc2))**2/ye1(i,j)**2
                                else
                                        delc=0.0
                                endif
                        else
                              if(y1(i,j).le.0.0.or.ye1(i,j).le.0.0.or.
     +                                rates(i+nc2).le.1.e-10) then
                                      delc=0.0
                                else
                                     ylog=log10(y1(i,j))
                                     rlog=log10(rates(i+nc2))
                                     elog=log10(y1(i,j)+ye1(i,j))-ylog
                                     delc=(ylog-rlog)**2/elog**2
                                endif
                        endif
                          chisq=chisq+delc
                enddo
            endif
        enddo
C loop for XRT data points
        elo=0.3
        ehi=10.0
        do j=1,m2
            if(ix2(j).ne.0) then
                call qra_grbrates(x2(j),npars,pars,-1,
     +          elo,ehi,nr,rates,flux)
                ncall=ncall+1
                do i=1,mc2
                        ii=mx2(j)+(i-1)*2
                        if(il.eq.0) then
                                delc=(y2(i,j)-rates(ii))**2/ye2(i,j)**2
                        else
                                if(y2(i,j).le.0.0.or.
     +                          rates(ii).le.1.e-10) then
                                        delc=0.0
                                else
                                      ylog=log10(y2(i,j))
                                      rlog=log10(rates(ii))
                                      elog=log10(y2(i,j)+ye2(i,j))-ylog
                                      delc=(ylog-rlog)**2/elog**2
                                endif
                        endif
                          chisq=chisq+delc
                enddo
            endif
        enddo
C        write(*,*) 'chi',chisq
        END
