/*****************************************************************************/
/*                                                                           */
/*  Tinflex_RC                                                               */
/*  Setup and sampling routine for generator                                 */
/*  (C version using R structure)                                            */
/*                                                                           */
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
/* header files */

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>

#include "Tinflex_RC.h"
#include "Tinflex_RC_arrays.h"
#include "Tinflex_lib_shared_private.h"

/*---------------------------------------------------------------------------*/
/* Prototypes of private functions                                           */

static double logpdf( SEXP lpdf, double x, SEXP env );
/*---------------------------------------------------------------------------*/
/* Evaluate log-density and its derivatives                                  */
/*---------------------------------------------------------------------------*/

/*****************************************************************************/

SEXP
Tinflex_RC_make_guide_table (SEXP sexp_ivs, SEXP sexp_Acum, SEXP sexp_gt)
/*---------------------------------------------------------------------------*/
/* Create guide table for draw interval at random.                           */
/* The result is stored in 'Acum' and 'gt'.                                  */
/*                                                                           */
/* Parameters:                                                               */
/*   ivs  ... data for hat and squeeze                                       */
/*   Acum ... cumulated areas (double vector)                                */
/*   gt   ... guide table (integer vector)                                   */
/*                                                                           */
/* Return:                                                                   */
/*   total area below hat.                                                   */
/*---------------------------------------------------------------------------*/
{
  double *ivs;                /* data for hat and squeeze for each interval  */
  int n_ivs;                  /* number of intervals                         */
  double *Acum;               /* array for storing cumulative areas          */
  int *gt;                    /* array for storing guide table               */ 
  SEXP sexp_A_tot = R_NilValue; /* R object for total area                   */
  double sum;                 /* variables for computing guide table         */
  int i,j;
  double Astep;

  /* Number of intervals. */
  n_ivs = length(sexp_ivs) / L_IVS - 1;

  /* Check arguments. */
  if (n_ivs<=0 || 
      length(sexp_Acum) != n_ivs || length(sexp_gt) != n_ivs ||
      ! IS_NUMERIC(sexp_ivs) || ! IS_NUMERIC(sexp_Acum) || ! IS_INTEGER(sexp_gt)) 
    error("Interval error. Please report.");

  /* Get array with data. */
  ivs = NUMERIC_POINTER(sexp_ivs);

  /* Get arrays for storing cumulative areas and guide table, resp. */
  Acum = NUMERIC_POINTER(sexp_Acum);
  gt = INTEGER_POINTER(sexp_gt);

  /* Compute cumulated areas (probabilities). */
  for ( i=0, sum = 0.; i <n_ivs; i++ )
    Acum[i] = (sum += *(ivs + i*L_IVS + IV_A_HT));

  /* Compute guide table. */
  Astep = Acum[n_ivs-1] / n_ivs;
  sum = 0.;
  for( j=0, i=0; j<n_ivs; j++ ) {
    while (Acum[i] < sum) 
      i++;
    if (i >= n_ivs) break;
    gt[j] = i;
    sum += Astep;
  }
  /* If there has been an round off error, */
  /* we have to complete the guide table.  */
  for( ; j<n_ivs; j++ )
    gt[j] = n_ivs - 1;

  /* Return total area to R. */
  PROTECT(sexp_A_tot = NEW_NUMERIC(1));
  NUMERIC_POINTER(sexp_A_tot)[0] = Acum[n_ivs-1];
  UNPROTECT(1);
  return sexp_A_tot;

} /* end of make_guide_table() */

/*---------------------------------------------------------------------------*/

SEXP
Tinflex_RC_sample (SEXP sexp_gen, SEXP sexp_n)
/*---------------------------------------------------------------------------*/
/* Draw sample from Tinflex generator object.                                */
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
  double *ivs;                /* data for hat and squeeze for each interval  */
  int n_ivs;                  /* number of intervals                         */
  double A_ht_tot;            /* total area below hat                        */
  SEXP sexp_lpdf;             /* log-density                                 */
  SEXP sexp_env;              /* R environment where log-density is evaluated */
  SEXP sexp_res = R_NilValue; /* R object for storing random sample          */
  int i;                      /* auxiliary loop variable                     */

  double *Acum;               /* cumulated areas                             */
  int *gt;                    /* guide table                                 */

  /* Acceptance-rejection loop: */
  double U,V;                 /* uniform random numbers                      */
  int J;                      /* index of random interval                    */
  double *iv;                 /* data for particular interval                */
  double a,b,y,x0,cT;         /* parameters for hat (tangent)                */
  double z;                   /* auxiliary variable                          */
  double t;                   /* hat (tangent) at construction point         */
  double X;                   /* random point                                */
  double hx, sx;              /* hat and squeeze at X                        */

  /* Extract and check sample size. */
  n = *(INTEGER (AS_INTEGER (sexp_n)));
  if (n<0) {
    error("sample size 'n' must be non-negative integer");
  }

  /* Get parameters generator:        */
  /*   parameters for hat and squeeze */
  ivs = NUMERIC_POINTER(VECTOR_ELT(sexp_gen, GEN_IVS));
  /*   number of intervals */
  n_ivs = length(VECTOR_ELT(sexp_gen, GEN_IVS)) / L_IVS - 1;
  /*   total area below hat */
  A_ht_tot = NUMERIC_POINTER(VECTOR_ELT(sexp_gen, GEN_A_HT_TOT))[0];
  /*   log-density */
  sexp_lpdf = VECTOR_ELT(sexp_gen, GEN_LPDF);
  /*   environment for computing log-density */
  sexp_env = VECTOR_ELT(sexp_gen, GEN_ENV);
  /*   parameters for hat and squeeze */
  ivs = NUMERIC_POINTER(VECTOR_ELT(sexp_gen, GEN_IVS));
  /*   cumulated areas */
  Acum = NUMERIC_POINTER(VECTOR_ELT(sexp_gen, GEN_ACUM));
  /*   guide table */
  gt = INTEGER_POINTER(VECTOR_ELT(sexp_gen, GEN_GT));

  /* Allocate memory for storing random sample. */
  PROTECT(sexp_res = NEW_NUMERIC(n));

  /* Get state for the R built-in URNG. */
  GetRNGstate();

  /* Draw sample of size 'n'. */
  for (i=0; i<n; i++) {
    
    while (1) {
      /* Acceptance-rejection loop. */
      
      /* Draw interval. */
      
      /* Sequential search:
       *   double sum = 0.;
       *   U = A_ht_tot * unif_rand();
       *   sum = 0.;
       *   for (J=0; J<n_ivs; J++) {
       * 	sum += *(ivs + J*L_IVS + IV_A_HT);
       * 	if (sum >= U) break;
       * }
       */

      /* Use guide table: */
      U = unif_rand();
      J = gt[(int)(U*n_ivs)];
      U *= A_ht_tot;
      for (; J<n_ivs; J++) {
      	if (Acum[J] >= U) break;
      }

      /* Get parameters for hat in interval. */
      iv = ivs + J*L_IVS;
      
      a  = iv[IV_HT_A];
      b  = iv[IV_HT_B];
      y  = iv[IV_HT_Y];
      x0 = iv[IV_X];
      cT = iv[IV_C];
      
      /* Value of transformed hat (tangent) at left boundary of interval. */
      t = a+b*(x0-y);
      
      /* Need a uniform random number in interval (0, iv[IV_A_HT]). */
      /*    U = iv[IV_A_HT] * unif_rand();                          */
      /* However, we can recycle U for this task.                   */
      U = *(ivs + J*L_IVS + IV_A_HT) + U - Acum[J];
      
      /* Generate from "hat distribution":                  */  
      /*    X = y + ( FTinv(cT, FT(cT,t) + b*U) - a ) / b;  */
      
      /* For numerical reasons we have to distinguish */
      /* between different values of 'cT'.            */
      if (cT == 0.) {
	/* Case: T(x)=log(x) */
	
	double expt = exp(t);
	z = U * b / expt;
	if (fabs(z) > 1.e-6) {
	  X = y + ( log(expt + b*U) - a ) / b;
	  /* or for |z| < Inf:                                */
	  /*   X <- x0 + U / exp(a+b*(x0-y)) * 1/z * log(1+z) */
	} else {
	  /* We need approximation by Taylor polynomial to avoid */
	  /* severe round-off errors. */
	  X = x0 + U / expt * (1 - z/2 + z*z/3);
	}
      }
      
      else if (cT == -0.5) {
	/* Case: T(x) = -1/sqrt(x) */
	
	z = U * b * t;
	if (fabs(z) > 1.e-6)
	  X = y + ( 1./(1./t - b*U) - a ) / b;
	else
	  X = x0 + U * t*t * (1 + z + z*z);
      }
      
      else if (cT == 1.) {
	/* Case: T(x) = x */
	
	z = U * b / (t*t);
	if (fabs(z) > 1.e-6)
	  X = y + (FTinv(cT, FT(cT, t)+b*U)-a)/b;
	else
	  X = x0 + U / t * (1 - z/2 + z*z/2);
      }
      
      else {
	/* Case: T(x)=sign(c)*x^c */
	/* For all other cases we only use a rough approximation in */
	/* case of numerical errors. */
	
	if (fabs(b)>1e-10) {
	  X = y + ( FTinv(cT, FT(cT,t) + b*U) - a ) / b;
	}	else {
	  U = U / iv[IV_A_HT];
	  X = (1.-U) * iv[IV_X] + U * iv[L_IVS + IV_X];
	}
      }
      
      /* Compute hat and squeeze at X. */
      hx = Tinv(cT, a + b * (X-y));
      if (iv[IV_A_SQ] > 0) 
	sx = Tinv(cT, iv[IV_SQ_A]+iv[IV_SQ_B]*(X-iv[IV_SQ_Y]));
      else
	sx = 0;
      
      /* Accept or reject. */
      V = hx * unif_rand();
      
      if (V <= sx) break;
      if (V <= exp(logpdf(sexp_lpdf,X,sexp_env))) break;
    }

    /* Store random point. */
    NUMERIC_POINTER(sexp_res)[i] = X; 
  }

  /* Update state for the R built-in URNG. */
  PutRNGstate();

  /* Return random sample to R. */
  UNPROTECT(1);
  return sexp_res;

} /* end of Tinflex_RC_sample() */

/*---------------------------------------------------------------------------*/

double
logpdf( SEXP lpdf, double x, SEXP env )
/*---------------------------------------------------------------------------*/
/* Evaluate log-density.                                                     */
/*---------------------------------------------------------------------------*/
{
  /* const struct Runuran_distr_discr *Rdistr; */
  SEXP R_fcall, arg;
  double y = 0.;

  /* Rdistr = unur_distr_get_extobj(distr); */
  PROTECT(arg = NEW_NUMERIC(1));
  NUMERIC_POINTER(arg)[0] = x;
  PROTECT(R_fcall = lang2(lpdf, arg));
  y = NUMERIC_POINTER(eval(R_fcall, env))[0];
  UNPROTECT(2);
  return y;
} /* end of logpdf() */

/*---------------------------------------------------------------------------*/
