 
/*-------------------------------------------------------------------------*/

#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

#include "Tinflex_RC.h"
#include "Tinflex_C.h"
#include "Tinflex_lib.h"

/*-------------------------------------------------------------------------*/

static const R_CallMethodDef CallEntries[] = {
    {"Tinflex_RC_make_guide_table", (DL_FUNC) &Tinflex_RC_make_guide_table, 3},
    {"Tinflex_RC_sample",   (DL_FUNC) &Tinflex_RC_sample,           2},
    {"Tinflex_C_setup",     (DL_FUNC) &Tinflex_C_setup,             9},
    {"Tinflex_C_sample",    (DL_FUNC) &Tinflex_C_sample,            2},
    {"Tinflex_C_2_R",       (DL_FUNC) &Tinflex_C_2_R,               1},
    {NULL, NULL, 0}
};

/*-------------------------------------------------------------------------*/

void R_init_Tinflex(DllInfo *dll)
{

  /* Declare some C routines to be callable from other packages */
  R_RegisterCCallable("Tinflex", "Tinflex_lib_setup",  (DL_FUNC) Tinflex_lib_setup);
  R_RegisterCCallable("Tinflex", "Tinflex_lib_sample", (DL_FUNC) Tinflex_lib_sample);
  R_RegisterCCallable("Tinflex", "Tinflex_lib_free",   (DL_FUNC) Tinflex_lib_free);

  /* Register native routines */
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}

/*-------------------------------------------------------------------------*/
