 
/*-------------------------------------------------------------------------*/

#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

#include "Tinflex.h"

/*-------------------------------------------------------------------------*/

static const R_CallMethodDef CallEntries[] = {
    {"Tinflex_sample",   (DL_FUNC) &Tinflex_sample,   2},
    {"make_guide_table", (DL_FUNC) &make_guide_table, 3},
    {NULL, NULL, 0}
};

/*-------------------------------------------------------------------------*/

void R_init_Tinflex(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}

/*-------------------------------------------------------------------------*/
