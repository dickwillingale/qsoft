*+QRT_PRTATHENA Trace proton through Athena telescope with magnetic diverter
        SUBROUTINE QRT_PRTATHENA(NPOLE,DM,PDX,PDY,PDZ,DDX,DDY,DDZ,
     +  EKV,TSIG,XAPER,XDIV,NRINGS,RRINGS,TRINGS,DRINGS,RDET,
     +  MAXST,XP,YP,ZP,NPATH,IQUAL)
        IMPLICIT NONE
        INTEGER NPOLE,NRINGS,MAXST,NPATH,IQUAL
        DOUBLE PRECISION DM(NPOLE),EKV,TSIG,XAPER,XDIV
        DOUBLE PRECISION RRINGS(NRINGS),TRINGS,DRINGS,RDET
        DOUBLE PRECISION PDX(NPOLE),PDY(NPOLE),PDZ(NPOLE)
        DOUBLE PRECISION DDX(NPOLE),DDY(NPOLE),DDZ(NPOLE)
        DOUBLE PRECISION XP(MAXST),YP(MAXST),ZP(MAXST)
Cf2py  intent(in) NPOLE,DM,PDX,PDY,PDZ,DDX,DDY,DDZ,EKV,TSIG
Cf2py  intent(in) XAPER,XDIV,NRINGS,RRINGS,TRINGS,DRINGS,RDET,MAXST
Cf2py  intent(out) XP,YP,ZP,NPATH,IQUAL
*NPOLE  input      number of dipoles
*DM     input      dipole moments Gauss cm3 
*PDX    input      x positions of dipoles cm
*PDY    input      y positions of dipoles cm
*PDZ    input      z positions of dipoles cm
*DDX    input      x directions of dipoles cm
*DDY    input      y directions of dipoles cm
*DDZ    input      z directions of dipoles cm
*EKV    input      proton energy keV
*TSIG   input      rms width of beam degrees
*XAPER  input      x position of mirror aperture cm (focal length)
*XDIV   input      x position of diverter input aperture cm
*NRINGS input      number of rings in diverter
*RRINGS input      radius of rings in diverter cm
*TRINGS input      radial thickness of rings in diverter cm
*DRINGS input      axial depth of rings in diverter cm
*RDET   input      radius of detector cm (axial position XDET=0.0)
*MAXST  input      maximum number of steps along path
*XP     output     x positions in path cm
*YP     output     y positions in path cm
*ZP     output     z positions in path cm
*NPATH  output     number of points along path
*IQUAL  output     path quality returned
C 0 hits active detector, 1 too close to dipole, 2 hits telescope tube
C 3 hits mirror aperture, 4 hits diverter, 5 hits focal plane beyond detector
C 6 maximum number of steps
*-Author Dick Willingale 2015-Jun-15
        LOGICAL STEPPING
        DOUBLE PRECISION DTOR,RWALL,RAMIN,RAMAX,BF
        DOUBLE PRECISION EV,XX,YY,ZZ,RAD,FX,FY,FZ,SM,RR,RRP,RD,RMM
        DOUBLE PRECISION DFX,DFY,DFZ,BX,BY,BZ,VDT,ETA,PSI,VETA
        EXTERNAL SYS_DRAND
        DOUBLE PRECISION RAN(2),SYS_DRAND,STEP_MAX,XDET,SS
        INTEGER I,K,NARC
C
        INCLUDE 'QR_COM'
        IF(ISTAT.NE.0) RETURN
C Coodinate system :-
C x-axis is optical axis, +ve from detector to mirrors
C centre of detector at XDET=0.0
C Input/output dimensions in cm
        XDET=0.0
C degrees to radians
        DTOR = 4.D0*ATAN(1.D0)/180.D0
C Set radii of mirror aperture and radius of aft body wall
C Note allow for converging beam using ratio XAPER/XDIV
        RAMIN=RRINGS(1)
        RAMAX=RRINGS(1)
        DO I=2,NRINGS
                RAMIN=MIN(RRINGS(I),RAMIN)
                RAMAX=MAX(RRINGS(I),RAMAX)
        ENDDO
        RWALL=(RAMAX+TRINGS*0.5)*XAPER/XDIV
        RAMAX=(RAMAX-TRINGS*0.5)*XAPER/XDIV
        RAMIN=RAMIN+TRINGS*0.5
C Set energy in eV
        EV=EKV*1000.0
C Maximum step 1cm
        STEP_MAX=10.0
C Number of points around arc for each step
        NARC=10
C Radius of curvature of electron path is r=sqrt(2 E m/e)/B
C where B is in Tesla, r is in metres, E in eV
C So for electrons of energy Eev require constant sqrt(2 m/e)=3.37E-6
C Proton mass 1837 times electron mass
C 1 tesla=1E4 Gauss so if B in Gauss we require factor 3.37*SQRT(1837)
C to get radius in cm
        ETA =3.37*SQRT(EV)*SQRT(1837.)
        PSI = 0.01*ETA*2.0*4.0*ATAN(1.D0)
C Choose random position within rear mirror aperture
        RAD=0.D0
        YY=-RAMAX+2.D0*RAMAX*SYS_DRAND()
        ZZ=-RAMAX+2.D0*RAMAX*SYS_DRAND()
        DO WHILE(RAD.LT.RAMIN.OR.RAD.GT.RAMAX)
              YY=-RAMAX+2.D0*RAMAX*SYS_DRAND()
              ZZ=-RAMAX+2.D0*RAMAX*SYS_DRAND()
              RAD=SQRT(YY**2+ZZ**2)
        ENDDO
        XX=XAPER
C direction towards the centre of the detector
        FX=XDET-XAPER
        FY=-YY
        FZ=-ZZ
        SM=SQRT(FX**2+FY**2+FZ**2)
        FX=FX/SM
        FY=FY/SM
        FZ=FZ/SM        
C Perturb using random numbers about central direction
        CALL SYS_GAUSS(2,RAN,0.0D0,TSIG*DTOR,ISTAT)
        FY = FY+RAN(1)
        FZ = FZ+RAN(2)
        SM=SQRT(FX**2+FY**2+FZ**2)
        FX=FX/SM
        FY=FY/SM
        FZ=FZ/SM        
C Loop to max number of steps
        STEPPING=.TRUE.
        NPATH=1
        XP(NPATH)=XX
        YP(NPATH)=YY
        ZP(NPATH)=ZZ
        DO WHILE(STEPPING)
C Field calculation
                CALL QRT_BFIELD(NPOLE,DM,PDX,PDY,PDZ,DDX,DDY,DDZ,1,
     +          XX,YY,ZZ,BX,BY,BZ,RMM)
                BF=SQRT(BX*BX+BY*BY+BZ*BZ)
                VDT =PSI/BF/NARC
                IF(VDT.GT.STEP_MAX/NARC) VDT=STEP_MAX/NARC
                RRP=SQRT(YY*YY+ZZ*ZZ)
                VETA = -VDT/ETA
C Make step by moving around arc at gyro radius
                DO K=1,NARC
C Direction of step given by cross-product F * B
                        DFX = VETA*(FY*BZ-FZ*BY)
                        DFY = VETA*(FZ*BX-FX*BZ) 
                        DFZ = VETA*(FX*BY-FY*BX)
                        XX = XX + FX*VDT
                        YY = YY + FY*VDT
                        ZZ = ZZ + FZ*VDT
                        FX=FX+DFX
                        FY=FY+DFY
                        FZ=FZ+DFZ
                ENDDO
C Save latest position
                NPATH=NPATH+1
                XP(NPATH)=XX
                YP(NPATH)=YY
                ZP(NPATH)=ZZ
                RR=SQRT(YY*YY+ZZ*ZZ)
C Test to see if hit stops or walls etc.
                IF(RMM.LT.0.1) THEN
C Too close to a dipole
                        STEPPING=.FALSE.
                        IQUAL=1
                ELSEIF(RR.GT.RWALL) THEN
C Outside telescope tube
                        STEPPING=.FALSE.
                        IQUAL=2
                ELSEIF(XP(NPATH).GT.XAPER) THEN
C Back to mirror aperture
                        STEPPING=.FALSE.
                        IQUAL=3
                ELSEIF(XP(NPATH-1).GE.XDIV.AND.
     +          XP(NPATH).LT.XDIV) THEN
C Crossed diverter input plane
                        RD=RRP+(RR-RRP)*
     +                  (XDIV-XP(NPATH-1))/(XP(NPATH)-XP(NPATH-1))
                        I=0
                        DO WHILE(STEPPING.AND.I.LT.NRINGS)
                              I=I+1
                              IF(RD.GT.(RRINGS(I)-TRINGS*0.5).AND.
     +                        RD.LT.(RRINGS(I)+TRINGS*0.5)) THEN
                                      STEPPING=.FALSE.
                                      IQUAL=4
                              ENDIF
                         ENDDO
                ELSEIF(XP(NPATH-1).GE.(XDIV-DRINGS).AND.
     +          XP(NPATH).LT.(XDIV-DRINGS)) THEN
C Crossed diverter output plane
                        RD=RRP+(RR-RRP)*
     +                  (XDIV-DRINGS-XP(NPATH-1))/
     +                  (XP(NPATH)-XP(NPATH-1))
                        I=0
                        DO WHILE(STEPPING.AND.I.LT.NRINGS)
                              I=I+1
                              IF(RD.GT.(RRINGS(I)-TRINGS*0.5).AND.
     +                        RD.LT.(RRINGS(I)+TRINGS*0.5)) THEN
                                      STEPPING=.FALSE.
                                      IQUAL=4
                              ENDIF
                         ENDDO
                ELSEIF(XP(NPATH).LE.XDET.AND.XP(NPATH-1).GT.XDET) THEN
C Crossed focal plane
                        STEPPING=.FALSE.
                        SS=(XDET-XP(NPATH-1))/(XP(NPATH)-XP(NPATH-1))
                        RD=RRP+(RR-RRP)*SS
                        XP(NPATH)=XDET
                        YP(NPATH)=YP(NPATH-1)+(YP(NPATH)-YP(NPATH-1))*SS
                        ZP(NPATH)=ZP(NPATH-1)+(ZP(NPATH)-ZP(NPATH-1))*SS
                        IF(RD.LE.RDET) THEN
C Hits active detector area
                                IQUAL=0
                         ELSE
                                IQUAL=5
                         ENDIF
                ELSEIF(NPATH.GE.MAXST) THEN
C Too many steps
                        STEPPING=.FALSE.
                        IQUAL=6
                ENDIF 
        ENDDO
        END
