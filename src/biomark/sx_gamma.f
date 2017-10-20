*+SX_GAMMA  Computes Gamma function for integers and half-integers
      REAL FUNCTION SX_GAMMA(X)
      REAL X
*X  input  argument will be an integer or half-integer
*-Author Clive Page 1993-OCT-27  (after routine in Bevington).
      INTEGER N, I
      REAL PROD, SUM, SX_FACTOR
*
      N = INT(X)
      IF(X .LE. N + 0.25) THEN
*argument is integer
         SX_GAMMA = SX_FACTOR(N)
      ELSE
*argument is half-integer
         PROD = 1.772455385
         IF(N .LE. 0) THEN
            SX_GAMMA = PROD
         ELSE IF(N .LE. 10) THEN
            DO I = 1,N
               PROD = PROD * (REAL(I) - 0.5)
            END DO
            SX_GAMMA = PROD
         ELSE
            SUM = 0.0
            DO I = 11,N
               SUM = SUM + LOG(REAL(I)-0.5)
            END DO
            SX_GAMMA = PROD * 639383.8623 * EXP(SUM)
         END IF
      END IF
      END

