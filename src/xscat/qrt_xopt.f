*+QRT_XOPT Calculate X-ray properties of a material
        SUBROUTINE QRT_XOPT(NC,MSPEC,RHO,NE,EKEV,ITYPE,
     +  ALPHA,GAMM,ABSL,F1,F2)
*NC        input                number of characters in specification
*MSPEC     input                specification of composition
*RHO       input                density gm/cm**3
*NE        input                number of photon energies
*EKEV      input                array of photon energies in keV
*ITYPE     input                data source (0=Cromer, 1=Henke)
*ALPHA     output               array of real parts dielectric constant
*GAMM      output               array of imaginary parts dielectric constant
*ABSL      output               array of absortion lengths cm-1
*F1        output               array of real parts scattering factor
*F2        output               array of imaginary parts scattering factor
        IMPLICIT NONE
        INTEGER NC,NE,ITYPE
        CHARACTER MSPEC*(*)
        DOUBLE PRECISION RHO,EKEV(NE)
        DOUBLE PRECISION ALPHA(NE),GAMM(NE),ABSL(NE),F1(NE),F2(NE)
Cf2py  intent(in) NC,MSPEC,RHO,NE,EKEV,ITYPE
Cf2py  intent(out) ALPHA,GAMM,ABSL,F1,F2
*-Author Dick Willingale 2012-May-7
        INTEGER MAT,LND
        PARAMETER (MAT=20)
        CHARACTER ATOM(MAT)*2
        CHARACTER DFILE*100,QSOFT*100
        INTEGER ICOMP(MAT),NAT,ID,IO,J
        DOUBLE PRECISION AKEV,WL
        PARAMETER (AKEV=12.397639)
        INTEGER ISTAT
C Parse material specification
        ISTAT=0
        CALL XX_DECODE(MSPEC(1:NC),MAT,ATOM,ICOMP,NAT,ISTAT)
        IF(ISTAT.NE.0) THEN
                WRITE(*,*) 'QRT_XOPT - material specification error'
                ISTAT=1
                RETURN
        ELSE
C Get atomic data
                CALL SYS_GETLUN(ID,ISTAT)
                IF(ITYPE.EQ.0) THEN
                        CALL SYS_GETENV('QSOFT',QSOFT,LND,ISTAT)
                        DFILE=QSOFT(1:LND)//'/data/cromer_liberman.dat'
                        OPEN(UNIT=ID,FILE=DFILE,STATUS='OLD')
                        CALL SYS_GETLUN(IO,ISTAT)
                        OPEN(UNIT=IO,STATUS='SCRATCH')
                        CALL XX_GET_CROMER(ID,IO,ATOM,ICOMP,NAT)
                        CLOSE(UNIT=ID)
                ELSEIF(ITYPE.EQ.1) THEN
                        CALL SYS_GETENV('QSOFT',QSOFT,LND,ISTAT)
                        DFILE=QSOFT(1:LND)//'/data/henke.dat'
                        OPEN(UNIT=ID,FILE=DFILE,STATUS='OLD')
                        CALL SYS_GETLUN(IO,ISTAT)
                        OPEN(UNIT=IO,STATUS='SCRATCH')
                        CALL XX_GET_HENKE(ID,IO,ATOM,ICOMP,NAT)
                        CLOSE(UNIT=ID)
                ELSE
                        CLOSE(UNIT=ID)
                        CLOSE(UNIT=IO)
                        WRITE(*,*) 'QRT_XOPT - unrecognised data srce'
                        ISTAT=1
                        RETURN
                ENDIF
                IF(ISTAT.EQ.0) THEN
                        DO J=1,NE
                                WL=AKEV/EKEV(J)
                                REWIND IO
                                CALL XX_OPTL(ITYPE,IO,0,ICOMP,RHO,
     +                          NAT,WL,
     +                          ALPHA(J),GAMM(J),ABSL(J),F1(J),F2(J))
                             ENDDO
                        CLOSE(UNIT=IO)
                ENDIF
        ENDIF
        END
