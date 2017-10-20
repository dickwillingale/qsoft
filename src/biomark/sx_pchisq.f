*+SX_PCHISQ  Returns upper-tail probability for Chi-squared function
      REAL FUNCTION SX_PCHISQ(CHISQ, NFREE)
      REAL CHISQ
      INTEGER NFREE
*CHISQ  input  Non-reduced chi-squared value
*NFREE  number of degrees of freedom
*-Author  Clive Page  1993-OCT-27  (after routine in Bevington)
      DOUBLE PRECISION Z, TERM, SUM
      REAL CHISQR, FREE, PWR
      INTEGER NEVEN, IMAX, I
      REAL SX_GAMMA
      EXTERNAL SX_GAMMA
*
      IF(NFREE .LE. 0) THEN
         SX_PCHISQ = 0.0
      ELSE
         CHISQR = CHISQ / NFREE
         FREE = NFREE
         Z = CHISQR * FREE / 2.0
         NEVEN = 2 * (NFREE/2)
         IF(NFREE .LE. NEVEN .OR. Z .GT. 25.0) THEN
            IF(Z .GT. 25.0) Z = CHISQR * NEVEN/2.0
            IMAX = NFREE / 2
            TERM = 1.0
            SUM  = 0.0
            DO I = 1,IMAX
               SUM = SUM + TERM
               TERM = TERM * Z / REAL(I)
            END DO
            SX_PCHISQ = REAL(SUM * EXP(-Z))
         ELSE
*odd number of degress of freedom
            PWR = FREE / 2.0
            TERM = 1.0
            SUM = TERM / PWR
            DO I = 1,1000
               TERM = -TERM * Z / REAL(I)
               SUM = SUM + TERM/(PWR + REAL(I))
               IF(ABS(TERM/SUM) .LE. 1.0E-5) GO TO 100
            END DO
100         CONTINUE
            SX_PCHISQ = REAL(1.0 - (Z**PWR) * SUM / SX_GAMMA(PWR))
         END IF
      END IF
      END

