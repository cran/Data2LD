#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include "DASarrayFnCpp.h"
#include "DBSarrayFnCpp.h"
#include "DRarrayFnCpp.h"
#include "loopJuanCpp.h"
#include "RmatFnCpp.h"
#include "SmatFnCpp.h"

static const R_CallMethodDef CallEntries[] = {
     {"DASarrayFnCpp", (DL_FUNC) &DASarrayFnCpp, 8},
     {"DBSarrayFnCpp", (DL_FUNC) &DBSarrayFnCpp, 8},
     {"DRarrayFnCpp",  (DL_FUNC) &DRarrayFnCpp,  6},
     {"loopJuanCpp",   (DL_FUNC) &loopJuanCpp,   5},
     {"RmatFnCpp",     (DL_FUNC) &RmatFnCpp,     7},
     {"SmatFnCpp",     (DL_FUNC) &SmatFnCpp,     9},
     {NULL, NULL, 0}
};

void R_init_Data2LD(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
