//
//  four-way tensor of inner products of four basis matrices
//  called by inprod.TPbasis3 within each iteration
//
//  Last updated by Jim Ramsay 6 March 2018
//

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

/* The gateway function */
SEXP loopJuan(SEXP BMAT1, SEXP BMAT2, SEXP BMAT3, SEXP BMAT4, SEXP WVEC)
{
    /* get numbers of columns */
    
    int ncols1 = ncols(BMAT1);
    int ncols2 = ncols(BMAT2);
    int ncols3 = ncols(BMAT3);
    int ncols4 = ncols(BMAT4);
    
    /*  get common number of rows from first matrix */
    
    int nrows1 = nrows(BMAT1);
    
    /* set up copies of input matrices  */
    /* do we need to do this? see p. 130 and Michelle's code in loop.c  */
    
    /* define points for copies of input matrices */
    
    SEXP BMAT1NEW;
    SEXP BMAT2NEW;
    SEXP BMAT3NEW;
    SEXP BMAT4NEW;
    SEXP WVECNEW;
    
    /* protect copies of input matrices */
    
    PROTECT(BMAT1NEW = coerceVector(BMAT1,REALSXP));
    PROTECT(BMAT2NEW = coerceVector(BMAT2,REALSXP));
    PROTECT(BMAT3NEW = coerceVector(BMAT3,REALSXP));
    PROTECT(BMAT4NEW = coerceVector(BMAT4,REALSXP));
    PROTECT( WVECNEW = coerceVector( WVEC,REALSXP));
    
    /* allocate memory for the four-way tensor */
    
    SEXP PTR;
    PROTECT(PTR = allocVector(REALSXP,ncols1*ncols2*ncols3*ncols4));
    
    /* assign 0 to all positions */
    
    for (int p = 0; p < (ncols1*ncols2*ncols3*ncols4); p++)
    {
        REAL(PTR)[p]=0;
    }
    
    /* loop through all columns of all four input matrices */
    
    for (int i = 0; i < ncols1; i++)
    {
        int mn1 = i*nrows1;
        
        for (int j = 0; j < ncols2; j++)
        {
            int mn2 = j*nrows1;
            
            for (int k = 0; k < ncols3; k++)
            {
                int mn3 = k*nrows1;
                
                for (int l = 0; l < ncols4; l++)
                {
                    int mn4 = l*nrows1;
                    
                    /* loop through rows of matrices */
                    
                    for (int m = 0; m < nrows1; m++)
                    {
                        
                        /* compute index in output vector to increment  */
                        
                        int index = i*ncols4*ncols3*ncols2
                                  + j*ncols4*ncols3
                                  + k*ncols4
                                  + l;
                        
                        /*  get values to be multiplied  */
                        
                        double wv   = REAL( WVECNEW)[m];
                        double arri = REAL(BMAT1NEW)[m + mn1];
                        double arrj = REAL(BMAT2NEW)[m + mn2];
                        double arrk = REAL(BMAT3NEW)[m + mn3];
                        double arrl = REAL(BMAT4NEW)[m + mn4];
                        
                        /*  increment value at index by the product */
                        
                        REAL(PTR)[index] += arri*arrj*arrk*arrl*wv;
                        
                    }
                }
                
            }
        }
    }
    
    UNPROTECT(6);
    return PTR;
    
}
