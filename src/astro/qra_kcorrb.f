*+QRA_KCORRB        K-correction using numerical integration of Band function
        SUBROUTINE QRA_KCORRB(e1src,e2src,n,e1obs,e2obs,z,gamma1,
     +  gamma2,ecobs,c,oi,si)
             implicit none
        integer n
        double precision z(n),e1obs(n),e2obs(n),e1src,e2src
        double precision gamma1(n),gamma2(n),ecobs(n)
        double precision c(n),oi(n),si(n)
*e1src        input                lower source frame energy
*e2src        input                upper source frame energy
*n        input                number of objects
*e1obs        input                lower observed energies
*e2obs        input                upper observed energies
*z        input                redshift of objects
*gamma1        input                photon indices below Ec
*gamma2        input                photon indices above Ec
*ecobs        input                observed Ec energies for each object
*c        output                correction factor for each object
*oi        output                integral over observed band
*si        output                integral over Eiso band in observed frame
*-Author Dick Willingale 2013-Feb-24
              include 'QR_COM'
        integer i
        double precision e1,e2,es1,es2,eb,g1,g2,f1,f2
C check QR status
        if(istat.ne.0) then
                write(*,*) 'qra_kcorr status',istat
                return
        endif
c loop for all GRBs
        do i=1,n
c observed energy range
                e1=e1obs(i)
                e2=e2obs(i)
c Eiso energy range in observed frame
                es1=e1src/(z(i)+1.0)
                es2=e2src/(z(i)+1.0)
c set photon spectrum parameters
                g1=-gamma1(i)
                g2=-gamma2(i)
                eb=ecobs(i)
                call qra_bandint(es1,es2,g1,g2,eb,f2)
                call qra_bandint(e1,e2,g1,g2,eb,f1)
c set k-correction factor
                oi(i)=f1
                si(i)=f2
                c(i)=f2/f1
        enddo
        end
