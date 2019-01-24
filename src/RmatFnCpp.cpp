#include <R.h>
#include <Rinternals.h>

extern "C" {
  SEXP RmatFnCpp(SEXP nXbasisw, SEXP nWbasisw, SEXP nXbasisx, SEXP nWbasisx,
                 SEXP BVECW,    SEXP BVECX,    SEXP BTENS) 
  {
    int nRmat;
    SEXP Rmat;
    int     nXw   = *INTEGER(nXbasisw);
    int     nWw   = *INTEGER(nWbasisw);
    int     nXx   = *INTEGER(nXbasisx);
    int     nWx   = *INTEGER(nWbasisx);
    double* Bvecw = REAL(BVECW);
    double* Bvecx = REAL(BVECX);
    double* Btens = REAL(BTENS);
    double* rRmat;
    nRmat = nXw*nXx;
    Rmat  = PROTECT(allocVector(REALSXP, nRmat));
    rRmat = REAL(Rmat);
    /* assign 0 to all positions */
    for (int p = 0; p < nRmat; p++) rRmat[p] = 0;
     /* loop through all columns of all four input matrices */
    for (int i = 0; i < nXw; i++)
    {
      for (int k = 0; k < nXx; k++)
      {
        for (int j = 0; j < nWw; j++)
        {
          for (int l = 0; l < nWx; l++)
          {
            int ijkl = i*nWx*nXx*nWw + j*nWx*nXx + k*nWx + l;
            /*  increment value at index by the product */
            rRmat[k*nXw+i] += Bvecw[j] * Bvecx[l] * Btens[ijkl];
          }
        }
      }
    }
    UNPROTECT(1);
    return Rmat;
  }
}
