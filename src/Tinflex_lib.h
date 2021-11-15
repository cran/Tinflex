/*****************************************************************************/
/*                                                                           */
/*   Tinflex -- C library                                                    */
/*                                                                           */
/*****************************************************************************/

/*---------------------------------------------------------------------------*/
/* Structures for storing Tinflex objects (hat and squeezes)                 */

/* prototype for log-density function                                        */
typedef double TINFLEX_FUNCT (double x, const void *params);

/* data for one interval                                                     */
struct Tinflex_iv {
  /* data for hat and squeeze: */
  double x;        /* left boundary of interval */
  double c;        /* parameter for transformation */
  double ht_a;     /* intercept for hat in transformed scale */
  double ht_b;     /* slope for hat in transformed scale */
  double ht_y;     /* anchor point of linear function */
  double sq_a;     /* intercept for squeeze in transformed scale */
  double sq_b;     /* slope for squeeze in transformed scale */
  double sq_y;     /* anchor point of linear function */
  double A_ht;     /* area below hat */
  double A_sq;     /* area below squeeze */
  /* data required for setup only: */
  int type;        /* type of interval */
  double Tfx;      /* log-density at left boundary  */
  double dTfx;     /* derivative of log-density at left boundary  */
  double d2Tfx;    /* 2nd derivative of log-density at left boundary */
  int next;        /* poor man's linked list: */
                   /* index of next interval to the right (-1 for last interval) */
};

typedef struct Tinflex_iv TINFLEX_IV;

/* Tinflex object                                                            */
struct Tinflex_gen {
  double *c;              /* parameter for transformation */
  TINFLEX_FUNCT *lpdf;    /* log-density */
  TINFLEX_FUNCT *dlpdf;   /* derivative of log-density */
  TINFLEX_FUNCT *d2lpdf;  /* second derivative of log-density */
  const void *params;     /* parameters for log-density */
  TINFLEX_IV *ivs;        /* data for hat and squeeze */
  int n_ivs;              /* number of intervals */
  double A_ht_tot;        /* total area below hat */
  double A_sq_tot;        /* total area below squeeze */
  double *Acum;           /* cumulated areas */
  int *gt;                /* guide table */
};

typedef struct Tinflex_gen TINFLEX_GEN;

/*---------------------------------------------------------------------------*/
/* Tinflex API                                                               */

TINFLEX_GEN *Tinflex_lib_setup (TINFLEX_FUNCT *lpdf, TINFLEX_FUNCT *dlpdf, TINFLEX_FUNCT *d2lpdf,
				const void *params,
				int n_ib, const double *ib,
				int n_c, const double *c,
				double rho, int max_intervals);
/* Setup: compute hat and squeeze for density.                               */

void * Tinflex_lib_free (TINFLEX_GEN *gen);
/* Free allocated memory.                                                    */

SEXP Tinflex_lib_sample (TINFLEX_GEN *gen, int n);
/* Draw sample from Tinflex generator object.                                */

double Tinflex_lib_sample_double (TINFLEX_GEN *gen);
/* Draw one random number from Tinflex generator object.                     */

/*---------------------------------------------------------------------------*/

