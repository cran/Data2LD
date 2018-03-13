//
//  Last updated by Jim Ramsay 6 March 2018
//

#include <stdio.h>

void SmatFn(int    *nXbasisw, int    *nWbasisw, int    *nUbasisj, int *nAbasisj,
            int    *nrep,     double *Bvecw,    double *Avecj,    double *Ucoefj,
            double *BAtenswj, double *Smat)
{
    int nXw = *nXbasisw;
    int nWw = *nWbasisw;
    int nUj = *nUbasisj;
    int nAj = *nAbasisj;
    int nRp = *nrep;
    int nSmat = nXw*nRp;
    /* assign 0 to all positions */
    for (int p = 0; p < nSmat; p++) Smat[p]=0;
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
                        Smat[r*nXw+i] +=
                        Bvecw[j] * Avecj[l] * Ucoefj[k+r*nUj] * BAtenswj[ijkl];
                    }
                }
            }
        }
    }
}
