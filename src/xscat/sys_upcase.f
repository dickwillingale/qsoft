*+SYS_UPCASE  Converts any lower-case letters in string to upper-case
      SUBROUTINE SYS_UPCASE(STRING)
      CHARACTER STRING*(*)
*STRING  in/out  String, returned in UPPER-case
*-Author  Clive Page  1991-Jan-1
*Restrictions: assumes ASCII code is in use, otherwise portable.
      INTEGER K
*
      DO 200,K = 1, LEN(STRING)
         IF(STRING(K:K) .GE. 'a' .AND. STRING(K:K) .LE. 'z')
     $       STRING(K:K) = CHAR(ICHAR(STRING(K:K))-32)
200   CONTINUE
      END
