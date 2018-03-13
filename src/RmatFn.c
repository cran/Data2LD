//
//  Last updated by Jim Ramsay 6 March 2018
//
//

#include <stdio.h>

void RmatFn(int    *nXbasisw, int    *nWbasisw, int    *nXbasisx, int    *nWbasisx,
            double *Bvecw,    double *Bvecx,    double *Btenswx,  double *Rmat)
{
    int nXw = *nXbasisw;
    int nWw = *nWbasisw;
    int nXx = *nXbasisx;
    int nWx = *nWbasisx;
    int nRmat = nXw*nXx;
    
    /* assign 0 to all positions */
    
    for (int p = 0; p < nRmat; p++) Rmat[p]=0;
    
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
                    
                    Rmat[k*nXw+i] += Bvecw[j] * Bvecx[l] * Btenswx[ijkl];
                    
                }                
            }
        }
    }
}
