*+XX_GAUSS		Function used by XX_CROMER
*-Author Dick Willingale
      FUNCTION XX_GAUSS(Y,M,LTBL)
C
C     GAUSSIAN-LEGENDRE QUADRATURE FORMULA FOR M POINTS
C     LTBL PROVIDES THE WEIGHTS AND ABSCISSAE
C     Y IS THE FUNCTION TO BE INTEGRATED
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LL=1
      G=0.
17    CALL LTBL(M,LL,A,Z)
      G=G+A*Y(Z)
      LL=LL+1
      IF(LL.LE.M) GOTO 17
      XX_GAUSS=G
      RETURN
      END
