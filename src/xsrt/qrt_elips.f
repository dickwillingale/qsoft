*+QRT_ELIPS        Set elliptical grazing incidence mirror
        SUBROUTINE TRS_ELIPS(ORG,AXS,CEN,XMIN,XMAX,AMIN,AMAX,SMB,RAB,
     +  IDE,IQ)
        IMPLICIT NONE
        DOUBLE PRECISION CEN(3),ORG(3),AXS(3)
        DOUBLE PRECISION XMIN,XMAX,AMIN,AMAX,SMB,RAB
        INTEGER IDE,IQ
Cf2py  intent(in) ORG,AXS,CEN,XMIN,XMAX,AMIN,AMAX,SMB,RAB,IDE,IQ
* Converted to qsoft 2017-May-3 RW
*-Author Charly Feldman 2006-Nov-14 based on trs_sphg created by
* Dick Willingale 1997-Nov-14
        INTEGER IDEF(2),J
        DOUBLE PRECISION RNM(3),TTT(3),RAD,PP(16)
        INCLUDE 'QR_COM'
C
        IF(ISTAT.NE.0) RETURN
C Calculate parameters
C Normal to surface at reference position and radius
        CALL SRT_DIDI(ORG,CEN,RNM,RAD)
C Reference axis in direction of optic axis
        CALL SRT_VCRS(AXS,RNM,TTT)
C check unit vectors
        CALL SRT_VNRM(RNM,ISTAT)
        CALL SRT_VNRM(TTT,ISTAT)
        CALL SRT_VNRM(AXS,ISTAT)
C Pack parameters into single array
        DO J=1,3
                PP(J)=AXS(J)
                PP(J+3)=RNM(J)
                PP(J+6)=CEN(J)
        ENDDO
        PP(10)=-(RAB**2)
        PP(11)=0
        PP(12)=SMB**2
        PP(13)=XMIN
        PP(14)=AMIN
        PP(15)=XMAX
        PP(16)=AMAX
C Get surface deformation and quality indices
        IDEF(1)=IDE
        IDEF(2)=1
C Set parameters in common for ellipse
        CALL SRT_SETF(0,19,16,PP,IDEF,IQ,-1,-1,ISTAT)
        END        
