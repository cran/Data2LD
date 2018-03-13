//
//
//  Created by James Ramsay on 2018-3-6.
//
//

#include <stdlib.h>
#include <R_ext/Rdynload.h>

#define CALLDEF(name, n) {#name, (DL_FUNC) &name, n}

extern void DASarrayFn(int    *nXbasisw, int    *nWbasisw, int    *nUbasisj, int *nAbasisj,
                       int    *nrep,     double *Bvecw,    double *Ucoefj,
                       double *BAtens,   double *DASarray);
extern void DBSarrayFn(int    *nXbasisw, int    *nWbasisw, int    *nUbasisj, int *nAbasisj,
                       int    *nrep,     double *Avecj,    double *Ucoefj,
                       double *BAtenswj, double *DBSarray);
extern void DRarrayFn(int   *nXbasisw, int    *nWbasisw, int    *nXbasisx, int *nWbasisx,
                      double *Bvecx,    double *Btenswx,  double *DRarray);
extern void RmatFn(int    *nXbasisw, int    *nWbasisw, int    *nXbasisx, int    *nWbasisx,
                   double *Bvecw,    double *Bvecx,    double *Btenswx,  double *Rmat);
extern void SmatFn(int    *nXbasisw, int    *nWbasisw, int    *nUbasisj, int *nAbasisj,
                   int    *nrep,     double *Bvecw,    double *Avecj,    double *Ucoefj,
                   double *BAtenswj, double *Smat);

static const R_CMethodDef R_CDef[] = {
  CALLDEF(DASarrayFn, 9),
  CALLDEF(DBSarrayFn, 9),
  CALLDEF(DRarrayFn, 7),
  CALLDEF(RmatFn, 8),
  CALLDEF(SmatFn, 10),
  {NULL, NULL, 0, NULL}
};

void
R_init_Data2LD(DllInfo *dll)
{
  R_registerRoutines(dll, R_CDef, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}

