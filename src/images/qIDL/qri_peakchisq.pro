;+
;
; NAME:
;   QRI_PEAKCHISQ
;
; PURPOSE:
;   This Library implememnts th QSOFT analysis package in IDL using calls to sharable objects
;   external to IDL. Original routines were in general written in Fortan 77 and are compiled
;   into sharable objects.
;
;   The purpose of this routine is to fit a lorentzian to the peak of the PSF
;
;
; :Categories:
;    Image analysis
;
; :Keywords:
;
; :Examples:
;     A = QRI_PEAKCHISQ(fpars)
;       *FPARS        input        fitting parameter values
;       *                parameter 1 normalisation (value at peak)
;       *                parameter 2 x-centre (pixel position)
;       *                parameter 3 y-centre (pixel position)
;       *                parameter 4 Lorentian width (pixels) 24th Sept. 2017 RW
;       * Deprecated     parameter 4 Gaussian sigma (pixels)
;       *CHISQ        output        Chi-squared value
;
; :Author:
;       Adrian Martindale
;
; :History:
;     Change History::
;       Written by: Adrian Martindale, 15/10/2018 to implement QSOFT in IDL.
;
;-
;

FUNCTION QRI_PEAKCHISQ, fpars

chisq = 0.d

  qsoft = GETENV('QSOFT')
  A=CALL_EXTERNAL(qsoft + '/R_libraries/imagesR.so', "qri_peakchisq_", $
      LONG(N_ELEMENTS(fpars)),$
      DOUBLE(FPARS),$
      DOUBLE(chisq),$
    /AUTO_GLUE)


RETURN,chisq
END
