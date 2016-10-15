*+QRA_NHTOT Estimate total Galactic NH for given line of sight
        subroutine qra_nhtot(equinox,np,radeg,decdeg,ssize,disio,
     +        nha,nhw,ebva,ebvw,nhma,nhmw,nhtota,nhtotw)
        implicit none
        integer np
        double precision equinox,radeg(np),decdeg(np),ssize,disio
        double precision nha(np),nhw(np),ebva(np),ebvw(np)
        double precision nhma(np),nhmw(np),nhtota(np),nhtotw(np)
*equinox        input        Equinox - default 2000 - same for all positions
*np             input        number of positions
*radeg          input        R.A. degrees array
*decdeg         input        Declination degrees array
*ssize          input        size of submap degrees - default 3
*disio          input        search distance degrees - default 1
*nha            output       NH1 LAB average
*nhw            output       NH1 LAB weighted average
*ebva           output       E(B-V) Schlegel et al. average
*ebvw           output       E(B-V) Schlegel et al. weighted average
*nhma           output       NH2 average
*nhmw           output       NH2 weighted average
*nhtota         output       NH total average
*nhtotw         output       NH total weighted average
*column densities in units of H atoms cm-2
*-Dick Willingale 2013-Jan-18
        include 'QR_COM'
        character*(255) nhmap,dumap,error
        integer lend,inhmap,idumap,ierr,LEN_TRIM
        integer block,i
        EXTERNAL LEN_TRIM
C
        if(ISTAT.NE.0) RETURN
C
        CALL SYS_GETENV('QSOFT',nhmap,lend,ISTAT)
        dumap=nhmap(1:lend)//'/data/dust_IREbv.fits'
        nhmap=nhmap(1:lend)//'/data/h1_nh_LAB.fits'
c
        CALL SYS_GETLUN(inhmap,ISTAT)
        lend=LEN_TRIM(nhmap)
        CALL ftopen(inhmap, nhmap(1:lend), 0, block, ISTAT)
        IF (ISTAT.NE.0) THEN
                write(*,*) 'Error opening input file ' //nhmap(1:lend),
     +                istat
                call FTCLOS(nhmap,ISTAT)
        ENDIF
        idumap=inhmap-1
c        CALL SYS_GETLUN(idumap,ISTAT)
        lend=LEN_TRIM(dumap)
        CALL ftopen (idumap, dumap(1:lend), 0, block, ISTAT)
        IF (ISTAT.NE.0) THEN
                write(*,*) 'Error opening input file ' //dumap(1:lend),
     +                istat
                call FTCLOS(idumap,ISTAT)
        ENDIF
C
        DO i=1,np
                CALL nhtot(equinox, radeg(i), decdeg(i), ssize,
     +          disio, inhmap,
     +          idumap, nha(i), nhw(i), ebva(i), ebvw(i), nhma(i),
     +                nhmw(i), nhtota(i), nhtotw(i), error, ISTAT)
        ENDDO
C
        CALL FTCLOS(inhmap,ISTAT)
        CALL FTCLOS(idumap,ISTAT)
        if(ISTAT.NE.0) THEN
                write(*,*) error(1:LEN_TRIM(error))
        endif
        end
