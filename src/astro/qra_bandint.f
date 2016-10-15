*+QRA_BANDINT        integrate the Band function
        SUBROUTINE QRA_BANDINT(e1,e2,g1,g2,eb,obsi)
             implicit none
        double precision e1,e2,g1,g2,eb,obsi
*e1        input                lower observed energy
*e2        input                upper observed energy
*g1        input                photon index below Ec
*g2        input                photon index above Ec
*eb        input                observed Ec 
*obsi        output                integral over observed band
*-Author Dick Willingale 2013-March-22
        integer j,ne
        parameter (ne=100)
        double precision x(ne),y(ne),xmin,xmax,xsam,ff
        external qra_band
              double precision qra_band
c set log10 energy range for integration
        xmin=log10(e1)
        xmax=log10(e2)
c set array of energies using logarithmic increments
        xsam=(xmax-xmin)/(ne-1)
        do j=1,ne
                x(j)=10.0**(xmin+xsam*(j-1))
C convert photons to flux
                y(j)=qra_band(x(j),-g1,-g2,eb)*x(j)
        enddo
c integrate flux
        obsi=0.0
        do j=2,ne
C add flux in interval
                ff=(y(j)+y(j-1))*(x(j)-x(j-1))*0.5
                obsi=obsi+ff
        enddo
        end
