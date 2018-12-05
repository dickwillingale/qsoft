*+SRT_PLNS         Intersection of a ray with a plane
*Specific deformations for Schmidt Lobster (SLE)
* Derived form SRT_PLNE by Vladimir Tichy
        SUBROUTINE SRT_PLNS(DIR,XA,PPT,PLN,XRF,P,IDEF,IKON,IHIT,ANS,RN,
     +        ISTAT)
        IMPLICIT NONE
        DOUBLE PRECISION DIR(3),XA(3),PPT(3),PLN(3),XRF(3),P(4)
        DOUBLE PRECISION ANS(3),RN(3)
        INTEGER IDEF(2),IKON,ISTAT,IHIT
        INTEGER NDI
*DIR    input   direction of ray
*XA     input   origin of ray
*PPT    input   point on plane
*PLN    input   normal to plane
*XRF        input        reference axis in plane
*P        input        limits
*IDEF        input        deformation indices
*IKON        input         configuration
*IHIT   output  annulus number hit (0 if not hit)
*ANS    output  intersection point on plane
*RN        output        normal to plane at intersection point (including deformation)
* The configuration is given by:
* IKON                deformation        limits
*  1                normal                cartesian
*  2                normal                radial
*  3                radial                radial
*  -ve                radial                nested radial
* If limits cartesian then P(1)=xmin, P(2)=ymin, P(3)=xmax, P(4)=ymax
* If limits radial then P(1)=rmin, P(2)=rmax, P(3)=rmin, P(4)=rmax, etc.
* If deformation is normal then this is used during search
* If deformation is radial then we perturb the radial limits
* If IKON -ve then ABS(IKON) is the number of nested annuli. In this
* case search for annulus which contains the ray. The radii must
* start with the largest (outer) first. 2 calls are required to SRT_DFRM,
*        Use IDEF(ISHELL) for outer radius
*        Use IDEF(ISHELL)+1 for inner radius
*ISTAT        in/out        returned status
        DOUBLE PRECISION DR(3),S,A,G,VP(3),R,X,Y,VY(3),DN,DM
        DOUBLE PRECISION DHDX,DHDY,DH,DL
        INTEGER J,JJ,IDL(3),NPASS
        LOGICAL SEARCH
        DOUBLE PRECISION SMALL
        PARAMETER (SMALL=1.E-4)
        INTEGER MAXPASS
C
        MAXPASS=10
        NDI=2
C
        IF(ISTAT.NE.0.0) RETURN
c        WRITE(*,*) 'DEBUG SRT_PLNS IDEF:',IDEF(1),IDEF(2)
c
c
C Calculate cos of angle between ray and normal
c
        CALL SRT_VDOT(PLN,DIR,A)
        IF(A.EQ.0.0) THEN
                IHIT=0
        ELSE
C Calculate direction and distance to point on plane
                CALL SRT_DIDI(XA,PPT,DR,S)
C Calculate cos of angle between line to point on plane and normal
                CALL SRT_VDOT(PLN,DR,G)
C Find y reference axis in plane
                CALL SRT_VCRS(PLN,XRF,VY)
C Start with no deformation in the normal direction
                DHDX=0.0
                DHDY=0.0
                DH=0.0
C Search for intersection
                SEARCH=.TRUE.
                NPASS=0
                DM=0.0
                DO WHILE(SEARCH)
C Calculate intersection position
                        DN=S*G+DH
C                        WRITE(*,*)  'SRT_PLNS DEBUG: DN=', DN
                        DO J=1,3
                                ANS(J)=DIR(J)*DN/A+XA(J)
                        ENDDO
C Find vector from origin on plane and intersection
                        CALL SRT_DIDI(PPT,ANS,VP,R)
C Calculate local x and y on plane
                        CALL SRT_VDOT(VP,XRF,X)
                        X=X*R
                        CALL SRT_VDOT(VP,VY,Y)
                        Y=Y*R
                        IF((IKON.EQ.1.OR.IKON.EQ.2).AND.IDEF(1).GT.0.AND.
     +                        ABS(DN-DM).GT.SMALL) THEN
                                DM=DN
C Get deformation
                                CALL SRT_MULTIDFRM(IDEF,X,Y,DH,DHDX,DHDY,ISTAT)
C                                CALL SRT_DFRM(IDEF,X,Y,DH,DHDX,DHDY,ISTAT)
                        ELSE
C Converged
                                SEARCH=.FALSE.
                        ENDIF
                        NPASS=NPASS+1
C                        IF(NPASS.EQ.10) THEN
C                                WRITE(*,*)  'SRT_PLNS warning: NPASS=10'
C     +                                 ,IDEF(1),IDEF(2)
C                        ENDIF
                        IF(NPASS.GT.MAXPASS) THEN
C                        IF(NPASS.LE.-1) THEN
                                WRITE(*,*)  'SRT_PLNS error - failed to converge'
     +                                 ,NPASS,IDEF(1),IDEF(2)
                                WRITE(*,*) 'local position',X,Y
                                SEARCH=.FALSE.
                        ENDIF
                ENDDO
C                WRITE(*,*)  'SRT_PLNS DEBUG ',IDEF(1),IDEF(2),NPASS
C Check limits
                IF(IKON.EQ.1) THEN
C                        IF (P(2).NE.-22.0.AND.IDEF(2).NE.0) THEN
C                                WRITE(*,*) 'DEBUG SRT_PLNS P()=',P(1),P(2),P(3),P(4)
C                                WRITE(*,*) 'IDEF=',IDEF(1),IDEF(2)
C                        ENDIF
                        IHIT=1
                        IF(X.LT.P(1).OR.X.GT.P(3)
     +                        .OR.Y.LT.P(2).OR.Y.GT.P(4)) THEN
                                IHIT=0
                        ENDIF
                        IF(DN/A.LT.SMALL) IHIT=0
                ELSE
                        WRITE(*,*) 'SRT_PLNS ERROR: UNKNOWN CONFIGURATION'
                ENDIF
                IF(IHIT.GT.0) THEN
C Calculate normal including deformation
                        DO J=1,3
                                RN(J)=PLN(J)-XRF(J)*DHDX-VY(J)*DHDY
                        ENDDO
                        CALL SRT_VNRM(RN,ISTAT)
                ENDIF
        ENDIF
        END
