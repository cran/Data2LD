#include <R.h>
#include <Rinternals.h>

extern "C" {
  SEXP DASarrayFnCpp(SEXP nXbasisw, SEXP nWbasisw, SEXP nUbasisj, SEXP nAbasisj, 
                     SEXP nrep,     SEXP BVECW,    SEXP UCOEFJ,   SEXP BATENS)
  {
    int nDASarray;
    SEXP DASarray;
    int nXw = *INTEGER(nXbasisw);
    int nWw = *INTEGER(nWbasisw);
    int nUj = *INTEGER(nUbasisj);
    int nAj = *INTEGER(nAbasisj);
    int nRp = *INTEGER(nrep);
    double* Bvecw  = REAL(BVECW);
    double* Ucoefj = REAL(UCOEFJ);
    double* BAtens = REAL(BATENS);
    double* rDASarray;
    nDASarray = nXw*nRp*nAj;
    DASarray  = PROTECT(allocVector(REALSXP, nDASarray));
    rDASarray = REAL(DASarray);
    /* assign 0 to all positions */
    for (int p = 0; p < nDASarray; p++)
    {
      rDASarray[p]=0;
    }
    /* loop through all columns of all four input matrices */
    for (int r = 0; r < nRp; r++)
    {
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
              rDASarray[i + r*nXw + l*nXw*nRp] +=
              Bvecw[j] * Ucoefj[k+r*nUj] * BAtens[ijkl];
            }
          }
        }
      }
    }
    UNPROTECT(1);
    return DASarray;
  }
}
