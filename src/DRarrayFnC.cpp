#include <R.h>
#include <Rinternals.h>

extern "C" {
  SEXP DRarrayFnC(SEXP nXbasisw, SEXP nWbasisw, SEXP nXbasisx, SEXP nWbasisx, 
                  SEXP BVECX,    SEXP BTENS)
  {
    int nDRarray;
    SEXP DRarray;
    int nXw = *INTEGER(nXbasisw);
    int nWw = *INTEGER(nWbasisw);
    int nXx = *INTEGER(nXbasisx);
    int nWx = *INTEGER(nWbasisx);
    double* Bvecx = REAL(BVECX);
    double* Btens = REAL(BTENS);
    double* rDRarray;
    nDRarray = nXw*nXx*nWw;
    DRarray = PROTECT(allocVector(REALSXP, nDRarray));
    rDRarray = REAL(DRarray);
     /* assign 0 to all positions */
    for (int p = 0; p < nDRarray; p++)
    {
      rDRarray[p]=0;
    }
    /* loop through all columns of all four input matrices */
    for (int i = 0; i < nXw; i++)
    {
      for (int j = 0; j < nWw; j++)
      {
        for (int k = 0; k < nXx; k++)
        {
          for (int l = 0; l < nWx; l++)
          {
            int ijkl = i*nWx*nXx*nWw + j*nWx*nXx + k*nWx + l;
            rDRarray[i+k*nXw+j*nXw*nXx] += Bvecx[l] * Btens[ijkl];
          }
        }
      }
    }
    UNPROTECT(1);
    return DRarray;
  }
}
