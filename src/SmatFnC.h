#include <R.h>
#include <Rinternals.h>

extern SEXP SmatFnC(SEXP nXbasisw, SEXP nWbasisw, SEXP nUbasisj, 
                    SEXP nAbasisj, SEXP nrep,     
                    SEXP BVECW, SEXP AVECJ, SEXP UCOEFJ, SEXP BATENS);
