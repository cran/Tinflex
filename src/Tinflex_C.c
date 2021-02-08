/*****************************************************************************/
/*                                                                           */
/*  Tinflex_C                                                                */
/*  Setup and sampling routine for generator                                 */
/*  (C only version)                                                         */
/*                                                                           */
/*****************************************************************************/

/* #define DEBUG 1 */

/*---------------------------------------------------------------------------*/
/* header files */

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>

#include "Tinflex_lib.h"
#include "Tinflex_RC_arrays.h"
#include "Tinflex_C.h"

/*---------------------------------------------------------------------------*/
/* Prototypes of private functions                                           */

/* Structure for storing R functions for calling in Tinflex_lib functions    */
typedef struct {
  SEXP lpdf;
  SEXP dlpdf;
  SEXP d2lpdf;
  SEXP env;
} TinflexC_lpdf ;

static void Tinflex_C_free (SEXP sexp_gen);
/*---------------------------------------------------------------------------*/
/* Free R object and all allocated memory blocks from C structure            */
/*---------------------------------------------------------------------------*/

static double Tinflex_C_eval_funct( int what, double x, const void *params);
static double eval_lpdf( double x, const void *params);
static double eval_dlpdf( double x, const void *params);
static double eval_d2lpdf( double x, const void *params);
/*---------------------------------------------------------------------------*/
/* Evaluate log-density and its derivatives                                  */
/*---------------------------------------------------------------------------*/

static SEXP Tinflex_C_tag(void); 
/*---------------------------------------------------------------------------*/
/* Create tag for R object                                                   */
/*---------------------------------------------------------------------------*/

/*****************************************************************************/

SEXP
Tinflex_C_setup (SEXP sexp_obj, SEXP sexp_env,
		 SEXP sexp_lpdf, SEXP sexp_dlpdf, SEXP sexp_d2lpdf,
		 SEXP sexp_ib, SEXP sexp_c,
		 SEXP sexp_rho, SEXP sexp_max_intervals)
/*---------------------------------------------------------------------------*/
/* Create Tinflex_C object                                                   */
/*                                                                           */
/* Parameters:                                                               */
/*   obj    ... R list                                                       */ 
/*   env    ... R environment                                                */
/*   lpdf   ... log-density                                                  */
/*   dlpdf  ... derivative of log-density                                    */
/*   d2lpdf ... 2nd derivative of log-density                                */
/*   ib     ... interval boundaries of decomposition                         */
/*   c      ... parameter for transformation (global or for each interval)   */
/*   rho    ... performance parameter: requested upper bound for ratio       */
/*              between area below hat and area below squeeze                */
/*   max_intervals ... maximal numbers of intervals                          */
/*---------------------------------------------------------------------------*/
/* Return:                                                                   */
/*   R object that contains pointer to Tinflex_lib C structure               */
/*---------------------------------------------------------------------------*/
{
  const double *ib;
  int n_ib;
  const double *c;
  int n_c;
  double rho;
  int max_intervals;

  TinflexC_lpdf *params;
  TINFLEX_GEN *gen;
  SEXP sexp_gen;

  /* extract arguments */
  ib = REAL(sexp_ib);
  n_ib = length(sexp_ib);
  c = REAL(sexp_c);
  n_c = length(sexp_c);
  rho = NUMERIC_VALUE(sexp_rho);
  max_intervals = INTEGER_VALUE(sexp_max_intervals);

  params = R_Calloc(1, TinflexC_lpdf);
  params->lpdf = sexp_lpdf;
  params->dlpdf = sexp_dlpdf;
  params->d2lpdf = sexp_d2lpdf;
  params->env = sexp_env;

  gen = Tinflex_lib_setup (eval_lpdf, eval_dlpdf, eval_d2lpdf, params,
  			   n_ib, ib, n_c, c, rho, max_intervals);

  /* make R external pointer and store pointer to structure */
  PROTECT(sexp_gen = R_MakeExternalPtr(gen, Tinflex_C_tag(), sexp_obj));
  
  /* register destructor as C finalizer */
  R_RegisterCFinalizer(sexp_gen, Tinflex_C_free);

  /* return pointer to R */
  UNPROTECT(1);

  return (sexp_gen);
  
} /* end of Tinflex_C_setup() */

/*---------------------------------------------------------------------------*/
  
SEXP
Tinflex_C_sample (SEXP sexp_gen, SEXP sexp_n)
/*---------------------------------------------------------------------------*/
/* Draw sample from Tinflex_C generator object.                               */
/*                                                                           */
/* Parameters:                                                               */
/*   gen ... generator object (list)                                         */
/*   n   ... sample size (positive integer)                                  */
/*                                                                           */
/* Return:                                                                   */
/*   random sample of size 'n'                                               */
/*---------------------------------------------------------------------------*/
{
  int n;                      /* sample size                                 */
  TINFLEX_GEN *gen;

  /* Extract and check sample size. */
  n = *(INTEGER (AS_INTEGER (sexp_n)));
  if (n<=0) {
    error("sample size 'n' must be positive integer");
  }

  gen = R_ExternalPtrAddr(sexp_gen);
  return Tinflex_lib_sample (gen, n);
  
} /* end of Tinflex_C_sample() */

/*---------------------------------------------------------------------------*/

void
Tinflex_C_free (SEXP sexp_gen)
/*----------------------------------------------------------------------*/
/* Free Tinflex_C object.                                                */
/*----------------------------------------------------------------------*/
{
  TINFLEX_GEN *gen;

#ifdef DEBUG
  Rprintf("Tinflex_C_free() called!\n");
#endif

  gen = R_ExternalPtrAddr(sexp_gen);
  R_Free(gen->params);
  Tinflex_lib_free (gen);
  R_ClearExternalPtr(sexp_gen);

} /* end of Tinflex_C_free() */

/*---------------------------------------------------------------------------*/

double
Tinflex_C_eval_funct( int what, double x, const void *params)
/*---------------------------------------------------------------------------*/
/* Evaluate log-density and its derivatives                                  */
/*---------------------------------------------------------------------------*/
{
  const TinflexC_lpdf *f = params;
  SEXP R_fcall, funct, arg;
  double y = 0.;

  switch(what) {
  case 0L:  /* log-density */
    funct = f->lpdf; break;
  case 1L:  /* 1st derivative */
    funct = f->dlpdf; break;
  case 2L:  /* 2nd derivative */
    funct = f->d2lpdf; break;
  default:
    error("internal error");
  }
  
  PROTECT(arg = NEW_NUMERIC(1));
  NUMERIC_POINTER(arg)[0] = x;
  PROTECT(R_fcall = lang2(funct, arg));
  y = NUMERIC_POINTER(eval(R_fcall, f->env))[0];
  UNPROTECT(2);
  
  return y;
} /* end of Tinflex_C_eval_funct() */

double eval_lpdf  ( double x, const void *params) { return Tinflex_C_eval_funct(0L, x, params); }
double eval_dlpdf ( double x, const void *params) { return Tinflex_C_eval_funct(1L, x, params); }
double eval_d2lpdf( double x, const void *params) { return Tinflex_C_eval_funct(2L, x, params); }

/*---------------------------------------------------------------------------*/

SEXP Tinflex_C_tag(void) 
/*---------------------------------------------------------------------------*/
/* Make tag for R Tinflex_C object                                           */
/*                                                                           */
/* Parameters: none                                                          */
/*                                                                           */
/* Return:                                                                   */
/*   tag (R object)                                                          */ 
/*---------------------------------------------------------------------------*/
{
  static SEXP tag = NULL;

  /* make tag for R object */
  if (!tag) tag = install("R_TINFLEX_C_TAG");

  return tag;
} /* end Tinflex_C_tag() */

/*---------------------------------------------------------------------------*/

#define MAX_LIST  (11)       /* maximum number of list entries */

struct Rlist {
  int len;                   /* length of list (depends on method) */
  char *names[MAX_LIST];     /* names of list elements */
  SEXP values;               /* pointer to R list */
};

void add_numeric(struct Rlist *list, char *key, double num)
{
  list->names[list->len] = key;
  SET_VECTOR_ELT(list->values, list->len, ScalarReal(num));
  ++list->len;
} /* end of add_numeric() */

void add_numeric_vec(struct Rlist *list, char *key, double *num, int n_num)
{
  int i;
  SEXP val;

  list->names[list->len] = key;

  val = NEW_NUMERIC(n_num);
  for (i=0; i<n_num; i++)
    REAL(val)[i] = num[i];
  SET_VECTOR_ELT(list->values, list->len, val);

  ++list->len;
} /* end of add_numeric_list() */

void add_integer_vec(struct Rlist *list, char *key, int *inum, int n_num)
{
  int i;
  SEXP val;

  list->names[list->len] = key;

  val = NEW_INTEGER(n_num);
  for (i=0; i<n_num; i++)
    INTEGER(val)[i] = inum[i];
  SET_VECTOR_ELT(list->values, list->len, val);

  ++list->len;
} /* end of add_integer_list() */

void add_ivs_data(struct Rlist *list, TINFLEX_GEN *gen)
{
  int i;
  TINFLEX_IV *iv;              /* pointer to current and next interval */
  SEXP val;

    /* data for intervals */
  list->names[list->len] = "ivs";
  val = allocMatrix(REALSXP, L_IVS, gen->n_ivs+1);
  for (i=0; i<gen->n_ivs+1; i++) {
    iv = gen->ivs + i;
    REAL(val)[i*L_IVS+IV_X] = iv->x;
    REAL(val)[i*L_IVS+IV_C] = iv->c;
    REAL(val)[i*L_IVS+IV_HT_A] = iv->ht_a;
    REAL(val)[i*L_IVS+IV_HT_B] = iv->ht_b;
    REAL(val)[i*L_IVS+IV_HT_Y] = iv->ht_y;
    REAL(val)[i*L_IVS+IV_SQ_A] = iv->sq_a;
    REAL(val)[i*L_IVS+IV_SQ_B] = iv->sq_b;
    REAL(val)[i*L_IVS+IV_SQ_Y] = iv->sq_y;
    REAL(val)[i*L_IVS+IV_A_HT] = iv->A_ht;
    REAL(val)[i*L_IVS+IV_A_SQ] = iv->A_sq;
    REAL(val)[i*L_IVS+IV_TYPE] = iv->type;
    REAL(val)[i*L_IVS+IV_TFX] = iv->Tfx;
    REAL(val)[i*L_IVS+IV_DTFX] = iv->dTfx;
    REAL(val)[i*L_IVS+IV_D2TFX] = iv->d2Tfx;
    REAL(val)[i*L_IVS+IV_NEXT] = iv->next + 1;
  }
  SET_VECTOR_ELT(list->values, list->len, val);
  ++list->len;
  
} /* end of add_ivs_data() */

/*...........................................................................*/

SEXP Tinflex_C_2_R (SEXP sexp_gen)
/*---------------------------------------------------------------------------*/
/* Copy entries from C structure into R list                                 */
/*---------------------------------------------------------------------------*/
{
  SEXP sexp_list;              /* pointer to R list */
  SEXP sexp_names;             /* array of keywords */
  struct Rlist list;           /* array of list elements */
  TINFLEX_GEN *gen;
  int i;
  
  /* extract Tinflex_C structure */
  gen = R_ExternalPtrAddr(sexp_gen);
  
  /* create temporary R list */
  PROTECT(list.values = allocVector(VECSXP, MAX_LIST));
  list.len = 0;

  /* copy data */
  add_numeric(&list,"A.ht.tot",gen->A_ht_tot);
  add_numeric(&list,"A.sq.tot",gen->A_sq_tot);
  add_numeric_vec(&list, "Acum", gen->Acum, gen->n_ivs);
  add_integer_vec(&list, "gt", gen->gt, gen->n_ivs);
  add_ivs_data(&list, gen);

  /* create final list */
  PROTECT(sexp_list = allocVector(VECSXP, list.len));
  for(i = 0; i < list.len; i++)
    SET_VECTOR_ELT(sexp_list, i, VECTOR_ELT(list.values, i));
    
  /* an array of the "names" attribute of the objects in our list */
  PROTECT(sexp_names = allocVector(STRSXP, list.len));
  for(i = 0; i < list.len; i++)
    SET_STRING_ELT(sexp_names, i,  mkChar(list.names[i]));
 
  /* attach attribute names */
  setAttrib(sexp_list, R_NamesSymbol, sexp_names);

  /* return list */
  UNPROTECT(3);
  return sexp_list;

} /* end of Tinflex_C_2_R() */

/*---------------------------------------------------------------------------*/
