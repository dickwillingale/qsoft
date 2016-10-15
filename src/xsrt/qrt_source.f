*+QRT_SOURCE Set source parameters
        SUBROUTINE QRT_SOURCE(IT,SD,SP,AP,AN,AR,AL,APRY,NR,IDEF)
*IT        input        source type 
*                1 point source at infinity               radial limits
*                2 point source at infinity               cartesian limits
*                3 point source at finite distance        radial limits
*                4 point source at finite distance        cartesian limits
*SD        input        source direction cosines
*SP        input        source position
*AP        input        aperture position
*AN        input        aperture normal
*AR        input        aperture reference axis
*AL        input        aperture limits
*APRY      input        area per ray - if 0 then use NRAY
*NR        input        number of rays
*IDEF      input        deformation index
Cf2py  intent(in) IT,SD,SP,AP,AN,AR,AL,APRY,NR,IDEF
        IMPLICIT NONE
        INTEGER IT,NR,IDEF
        DOUBLE PRECISION SD(3),SP(3),AP(3),AN(3),AR(3),AL(6),APRY
*-Author Dick Willingale 2012-Apr-30
        INCLUDE 'SRT_COM'
        DOUBLE PRECISION AA,SA
        DOUBLE PRECISION ALIM(6)
        INCLUDE 'QR_COM'
C Check status
        IF(ISTAT.NE.0) return
C
        IF(MOD(IT,2).EQ.0) THEN
C Find aperture area cartesian limits
                AA=(AL(3)-AL(1))*(AL(4)-AL(2))
                IF(AA.LE.0.0) THEN
                   WRITE(*,*) 'QRT_SOURCE error - bad cartesiam limits'
                   ISTAT=1
                   RETURN
                ENDIF
                ALIM(1)=AL(1)
                ALIM(2)=AL(2)
                ALIM(3)=AL(3)
                ALIM(4)=AL(4)
                ALIM(5)=0.0
                ALIM(6)=0.0
        ELSE
C Find aperture area radial limits
                AA=PI*(AL(2)**2-AL(1)**2)
                IF(AA.LE.0.0) THEN
                   WRITE(*,*) 'QRT_SOURCE error - bad radial limits'
                   ISTAT=1
                   RETURN
                ENDIF
                ALIM(1)=AL(1)
                ALIM(2)=AL(2)
                ALIM(3)=0.0
                ALIM(4)=0.0
                ALIM(5)=0.0
                ALIM(6)=0.0
        ENDIF
        SA=APRY
        IF(SA.EQ.0.0) THEN
C Specify number of rays and calculate area per ray
                IF(NR.GT.0) THEN
                        SA=AA/NR
                        NRAYS=NR
                ELSE
                    WRITE(*,*) 'QRT_SOURCE error - number of rays <= 0'
                    ISTAT=1
                    RETURN
                ENDIF
        ELSEIF(SA.GT.0.0) THEN
                NRAYS=NINT(AA/APRY)
                SA=AA/NRAYS
        ELSE
                WRITE(*,*) 'QRT_SOURCE error - negative area per ray'
                ISTAT=1
                RETURN
        ENDIF
C Finally actually set the parameters in common
        CALL SRT_SETS(IT,IDEF,SD,SP,SA,AP,AN,AR,ALIM,ISTAT)
        END
