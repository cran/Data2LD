#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include "DASarrayFnC.h"
#include "DBSarrayFnC.h"
#include "DRarrayFnC.h"
#include "loopJuanC.h"
#include "RmatFnC.h"
#include "SmatFnC.h"

static const R_CallMethodDef CallEntries[] = {
     {"DASarrayFnC", (DL_FUNC) &DASarrayFnC, 8},
     {"DBSarrayFnC", (DL_FUNC) &DBSarrayFnC, 8},
     {"DRarrayFnC",  (DL_FUNC) &DRarrayFnC,  6},
     {"loopJuanC",   (DL_FUNC) &loopJuanC,   5},
     {"RmatFnC",     (DL_FUNC) &RmatFnC,     7},
     {"SmatFnC",     (DL_FUNC) &SmatFnC,     9},
     {NULL, NULL, 0}
};

void R_init_Data2LD(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
