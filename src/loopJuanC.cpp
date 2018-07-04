#include <R.h>
#include <Rinternals.h>

extern "C" { 
  SEXP loopJuanC(SEXP BVEC1, SEXP BVEC2, SEXP BVEC3, SEXP BVEC4, SEXP WVEC)
  {
    int nBtensor;
    SEXP Btensor;
    
    int nrow = length(WVEC);
    int nbasis1 = length(BVEC1)/nrow;
    int nbasis2 = length(BVEC2)/nrow;
    int nbasis3 = length(BVEC3)/nrow;
    int nbasis4 = length(BVEC4)/nrow;
    double* Bvec1 = REAL(BVEC1); 
    double* Bvec2 = REAL(BVEC2); 
    double* Bvec3 = REAL(BVEC3); 
    double* Bvec4 = REAL(BVEC4);
    double* Wvec  = REAL(WVEC);
    double* rBtensor;
    
    // allocate memory for the four-way tensor //
    
    nBtensor = nbasis1*nbasis2*nbasis3*nbasis4;
    
    Btensor = PROTECT(allocVector(REALSXP, nBtensor));
    rBtensor = REAL(Btensor);
    
    
    /* assign 0 to all positions */
    
    for (int p = 0; p < (nBtensor); p++)
    {
      rBtensor[p]=0;
    }
    
    /* loop through all columns of all four input matrices */
    
    for (int i = 0; i < nbasis1; i++)
    {
      int mn1 = i*nrow;
      
      for (int j = 0; j < nbasis2; j++)
      {
        int mn2 = j*nrow;
        
        for (int k = 0; k < nbasis3; k++)
        {
          int mn3 = k*nrow;
          
          for (int l = 0; l < nbasis4; l++)
          {
            int mn4 = l*nrow;
            
            /* loop through rows of matrices */
            
            for (int m = 0; m < nrow; m++)
            {
              
              /* compute index in output vector to increment  */
              
              int index = i*nbasis4*nbasis3*nbasis2
              + j*nbasis4*nbasis3
              + k*nbasis4
              + l;
              
              /*  increment value at index by the product */
              
              rBtensor[index] += Bvec1[m + mn1]*Bvec2[m + mn2]*
              Bvec3[m + mn3]*Bvec4[m + mn4]*Wvec[m];
              
            }
          }
          
        }
      }
    }
    UNPROTECT(1);
    return Btensor;
  }
}
