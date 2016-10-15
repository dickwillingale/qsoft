*+SYS_RAN0 Numerical Recipes random number generator
      DOUBLE PRECISION FUNCTION SYS_DRAN0(idum)
* Minimal random number generator of Park and  Miller. Returns
* a uniform random deviate between 0.0 and 1.0. Set or reset
* idum to any integer value (except MASK) to initialise the
* sequence; idum must not be altered between calls for 
* successive deviates in a sequence.
*-from Mike Denby 1997-Nov-28
      INTEGER idum,IA,IM,IQ,IR,MASK
      DOUBLE PRECISION  AM
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     *MASK=123459876)
      INTEGER k
      idum=ieor(idum,MASK)
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      sys_dran0=AM*idum
      idum=ieor(idum,MASK)
      return
      END
