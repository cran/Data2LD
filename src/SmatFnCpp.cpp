#include <R.h>
#include <Rinternals.h>

extern "C" {
  SEXP SmatFnCpp(SEXP nXbasisw, SEXP nWbasisw, SEXP nUbasisj, 
                 SEXP nAbasisj, SEXP BVECW, SEXP AVECJ, SEXP UCOEFJ, 
                 SEXP BATENS)
  {
    int nSmat;
    SEXP Smat;
    int nXw = *INTEGER(nXbasisw);
    int nWw = *INTEGER(nWbasisw);
    int nUj = *INTEGER(nUbasisj);
    int nAj = *INTEGER(nAbasisj);
    double* Bvecw  = REAL(BVECW);
    double* Avecj  = REAL(AVECJ);
    double* Ucoefj = REAL(UCOEFJ);
    double* BAtens = REAL(BATENS);
    double* rSmat;
    nSmat = nXw;
    Smat = PROTECT(allocVector(REALSXP, nSmat));
    rSmat = REAL(Smat);
    /* assign 0 to all positions */
    for (int p = 0; p < nSmat; p++) rSmat[p]=0;
    /* loop through all columns of all four input matrices */
    for (int i = 0; i < nXw; i++)
    {
      for (int j = 0; j < nWw; j++)
      {
        for (int k = 0; k < nUj; k++)
        {
          for (int l = 0; l < nAj; l++)
          {
            int ijkl = i*nAj*nUj*nWw + j*nAj*nUj + k*nAj + l;
            /*  increment value at index by the product */
            rSmat[i] += Bvecw[j] * Avecj[l] * Ucoefj[k] * BAtens[ijkl];
          }
        }
      }
    }
    UNPROTECT(1);
    return Smat;
  }
}

