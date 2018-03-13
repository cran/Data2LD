//
//  DRarrayFn.c
//  
//
//  Created by James Ramsay on 2018-3-6.
//
//

#include <stdio.h>

void DRarrayFn(int   *nXbasisw, int    *nWbasisw, int    *nXbasisx, int *nWbasisx,
              double *Bvecx,    double *Btenswx,  double *DRarray)
{
    
    int nXw = *nXbasisw;
    int nWw = *nWbasisw;
    int nXx = *nXbasisx;
    int nWx = *nWbasisx;
    
    int nDRarray = nXw*nXx*nWw;
    
    /* assign 0 to all positions */
    
    for (int p = 0; p < nDRarray; p++)
    {
        DRarray[p]=0;
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
                    
                    /*  increment value at index by the product */
                    
                    DRarray[i+k*nXw+j*nXw*nXx] += Bvecx[l] * Btenswx[ijkl];
                    
                }
                
            }
        }
    }
    
}
