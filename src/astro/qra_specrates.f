*+QRA_SPECRATES        use lookup table to predict count rates for sprectral function
        SUBROUTINE QRA_SPECRATES(p1,p2,p3,nr,rates,istat)
        implicit none
        integer nr,istat
        double precision p1,p2,p3,rates(nr)
*p1                parameter 1
*p2                parameter 2
*p3                parameter 3
*nr                number of energy channels (8, 1-4 XRT, 5-8 BAT)
*rates                predicted count rates
*istat                return status
*-Author Dick Willingale 2008-Oct-11
        include 'GRBS_COM'
        character*(80) tablefile
        integer i,ierr,j1,j2,nsam,j,iu
        integer len_trim
        external len_trim
        if(na.eq.0) then
C Get lookup table
                CALL SYS_GETLUN(IU,ISTAT)
                tablefile=grbname(1:len_trim(grbname))//'_table.dat'
                OPEN(UNIT=IU,FILE=tablefile(1:len_trim(tablefile)),
     +                STATUS='OLD',IOSTAT=IERR)
                if(ierr.ne.0) then
                        write(*,*) 'IOERR - failed to open '//tablefile
                        istat=1
                        return
                endif
                READ(IU,*) is,na,nb,nc
                READ(IU,*) amin,amax,bmin,bmax,cmin,cmax
                nsam=na*nb*nc
                if(nsam*6.gt.maxr) then
                        CLOSE(UNIT=IU)
                        write(*,*) 'DIMERR - GRBRATES too many rates'
                        istat=1
                        return
                endif
                do i=1,nsam
                        j1=(i-1)*nr+1
                        j2=j1+nr-1
                        READ(IU,*,IOSTAT=IERR) (table(j),j=j1,j2)
                        if(ierr.ne.0) then
                               write(*,*) 'IOERR - reading '//tablefile
                                istat=1
                                return
                        endif
                enddo
                CLOSE(UNIT=IU)
        endif
        if(is.eq.1) then
C is=1 Band function 
C p1        alpha                low photon index  (with sign)
C p2        amb                change in index alpha-beta (with sign)
C p3        ecut                exponential cut-off energy
               call qra_specinterp(p1,p2,log10(p3),amin,amax,bmin,bmax,
     +                cmin,cmax,na,nb,nc,table,nr,rates)
             elseif(is.eq.2) then
C is=2 broken power law
C p1        gamma1                low photon index (no sign)
C p2        gamma2-gamma1        change in photon index (no sign)
C p3        ebr                break energy
          call qra_specinterp(-p1,-p1+p2,log10(p3),amin,amax,bmin,bmax,
     +                cmin,cmax,na,nb,nc,table,nr,rates)
             elseif(is.eq.3) then
C is=3 truncated  power law
C p1        gamma                photon index (with sign)
C p2        et                low trunction energy keV
C p3        ec                high cut-off energy keV
            call qra_specinterp(p1,log10(p2),log10(p3),amin,amax,bmin,
     +                 bmax,cmin,cmax,na,nb,nc,table,nr,rates)
             elseif(is.eq.4) then
C is=4 smoothed broken power law
C p1        gamma1                low photon index (with sign)
C p2        si                smoothing index
C p3        eb                break energy
C Note gamma2=1.0 fixed
             call qra_specinterp(p1,p2,log10(p3),amin,amax,bmin,bmax,
     +                cmin,cmax,na,nb,nc,table,nr,rates)
        endif
            end
        SUBROUTINE QRA_SPECINTERP(a,b,c,amin,amax,bmin,bmax,cmin,cmax,
     +        na,nb,nc,table,nr,rates)
        integer na,nb,nc,nr
        double precision a,b,c,amin,amax,bmin,bmax,cmin,cmax
        double precision table(nr,na,nb,nc),rates(nr)
        double precision asam,bsam,csam,afr,bfr,cfr,r1,r2,r3,r4,r5,r6
        integer i,ia,ib,ic,ia1,ib1,ic1,ilow
C Calculate cell size
        asam=(amax-amin)/(na-1)
        bsam=(bmax-bmin)/(nb-1)
        csam=(cmax-cmin)/(nc-1)
C Calculate position of cell in lookup table
        ia=max(min(int((a-amin)/asam),na-2),0)+1
        ib=max(min(int((b-bmin)/bsam),nb-2),0)+1
        ic=max(min(int((c-cmin)/csam),nc-2),0)+1
C Calculate fractional distances across cell
        afr=(a-(amin+asam*(ia-1)))/asam
        bfr=(b-(bmin+bsam*(ib-1)))/bsam
        cfr=(c-(cmin+csam*(ic-1)))/csam
        if(afr.le.-0.001.or.afr.gt.1.001) then
C                write(*,*) 'WARNING - a interpolation outside table',a
C                afr=max(min(afr,1.0),0.0)
        endif
        if(bfr.le.-0.001.or.bfr.gt.1.001) then
C                write(*,*) 'WARNING - b interpolation outside table',b
C                bfr=max(min(bfr,1.0),0.0)
        endif
        ilow=0
        if(cfr.le.-0.001.or.cfr.gt.1.001) then
C                write(*,*) 'WARNING - c interpolation outside table',c,cfr
                ilow=1
C                cfr=max(min(cfr,1.0),0.0)
        endif
C Linear interpolate using corners of cell
        ia1=ia+1
        ib1=ib+1
        ic1=ic+1
        do i=1,nr
                if(afr.gt.bfr) then
                        r1=table(i,ia,ib,ic)
                        r2=table(i,ia,ib,ic1)
                        r3=table(i,ia1,ib,ic)
                        r4=table(i,ia1,ib1,ic)
                        r5=table(i,ia1,ib,ic1)
                        r6=table(i,ia1,ib1,ic1)
                        if(r1.gt.0.0.and.r2.gt.0.0.and.r3.gt.0.0.
     +                  and.r4.gt.0.0.and.r5.gt.0.0.and.r6.gt.0.0) then
                                     r1=log10(r1)
                                r2=log10(r2)
                                r3=log10(r3)
                                r4=log10(r4)
                                r5=log10(r5)
                                r6=log10(r6)
                                r1=r1+(r3-r1)*afr+(r4-r3)*bfr
                                r2=r2+(r5-r2)*afr+(r6-r5)*bfr
                                rates(i)=10.0**(r1+(r2-r1)*cfr)
                        else
                                r1=r1+(r3-r1)*afr+(r4-r3)*bfr
                                r2=r2+(r5-r2)*afr+(r6-r5)*bfr
                                rates(i)=max(r1+(r2-r1)*cfr,0.0)
                        endif
                else
                        r1=table(i,ia,ib,ic)
                        r2=table(i,ia,ib,ic1)
                        r3=table(i,ia,ib1,ic)
                        r4=table(i,ia1,ib1,ic)
                        r5=table(i,ia,ib1,ic1)
                        r6=table(i,ia1,ib1,ic1)
                        if(r1.gt.0.0.and.r2.gt.0.0.and.r3.gt.0.0.
     +                and.r4.gt.0.0.and.r5.gt.0.0.and.r6.gt.0.0) then
                                     r1=log10(r1)
                                r2=log10(r2)
                                r3=log10(r3)
                                r4=log10(r4)
                                r5=log10(r5)
                                r6=log10(r6)
                                r1=r1+(r3-r1)*bfr+(r4-r3)*afr
                                r2=r2+(r5-r2)*bfr+(r6-r5)*afr
                                rates(i)=10.0**(r1+(r2-r1)*cfr)
                        else
                                r1=r1+(r3-r1)*bfr+(r4-r3)*afr
                                r2=r2+(r5-r2)*bfr+(r6-r5)*afr
                                rates(i)=max(r1+(r2-r1)*cfr,0.0)
                        endif
                endif
C                if(ilow.eq.1) then
C                        rates(i)=0.0
C                endif
        enddo
        END
