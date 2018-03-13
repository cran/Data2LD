//  
//
//  Created by James Ramsay on 2018-3-6.
//
//

#include <stdio.h>

void DBSarrayFn(int    *nXbasisw, int    *nWbasisw, int    *nUbasisj, int *nAbasisj,
                int    *nrep,     double *Avecj,    double *Ucoefj,
                double *BAtenswj, double *DBSarray)
{
    
    int nXw = *nXbasisw;
    int nWw = *nWbasisw;
    int nUj = *nUbasisj;
    int nAj = *nAbasisj;
    int nRp = *nrep;
    
    int nDBSarray = nXw*nRp*nWw;
    
    /* assign 0 to all positions */
    
    for (int p = 0; p < nDBSarray; p++)
    {
        DBSarray[p]=0;
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
                        
                        DBSarray[i+r*nXw+j*nRp*nXw] +=
                               Avecj[l] * Ucoefj[k+r*nUj] * BAtenswj[ijkl];
                        
                    }
                }
            }
        }
    }
}
