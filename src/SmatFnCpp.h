#include <R.h>
#include <Rinternals.h>

extern SEXP SmatFnCpp(SEXP nXbasisw, SEXP nWbasisw, SEXP nUbasisj, 
                      SEXP nAbasisj, SEXP BVECW, SEXP AVECJ, SEXP UCOEFJ, 
                      SEXP BATENS);
