/*****************************************************************************/
/*                                                                           */
/*  Tinflex C API                                                            */
/*  Setup, sampling and freeing                                              */
/*                                                                           */
/*****************************************************************************/

/* #define DEBUG 1 */

/*---------------------------------------------------------------------------*/
/* header files */

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>

#include "Tinflex_lib.h"
#include "Tinflex_lib_shared_private.h"

/*---------------------------------------------------------------------------*/
/* Prototypes of private functions                                           */
/*---------------------------------------------------------------------------*/

static int cmpiv (const void *p1, const void *p2);
/* auxiliary function for sorting intervals w.r.t. is 'x' value */

/* static double Tf (double x, TINFLEX_GEN *gen, TINFLEX_IV *iv); */
/* Compute transformed density by means of its log-density.                  */

static int Tfdd (double x, TINFLEX_GEN *gen, TINFLEX_IV *iv);
/* Compute transformed density and its derivatives                           */
/* by means of its log-density.                                              */

static double arcmean( double x0, double x1 );
/* compute "arctan mean" of two numbers.                                     */

static int type_iv (TINFLEX_IV *left, TINFLEX_IV *right);
/* Determine type of interval.                                               */

static double area (double c, double a, double b, double y, double from, double to);
/* Compute area underneath hat or squeeze in particular interval.            */

static double do_area (double c, double a, double b, double y, double from, double to);
/* Perform computation of area underneath hat or squeeze.                    */

static int hat_iv (TINFLEX_IV *left, TINFLEX_IV *right, int link);
/* Compute hat and squeeze for a paricular interval.                         */

static int make_guide_table (TINFLEX_GEN *gen);
/* Create guide table for drawing interval at random.                        */

#ifdef DEBUG
static void debug_iv (TINFLEX_IV *iv);
#endif

/*****************************************************************************/
/*                                                                           */
/* Transformed Density                                                       */
/*                                                                           */
/*****************************************************************************/

/* double */
/* Tf (double x, TINFLEX_GEN *gen, TINFLEX_IV *iv) */
/* /\*---------------------------------------------------------------------------*\/ */
/* /\* Compute transformed density by means of its log-density.                  *\/ */
/* /\*---------------------------------------------------------------------------*\/ */
/* /\*   x   ... argument                                                        *\/ */
/* /\*   gen ... Tinflex generator object                                        *\/ */
/* /\*---------------------------------------------------------------------------*\/ */
/* { */
/*   double lfx, Tfx; */
/*   double c = iv->c; */
/*   double s; */

/*   Tfx = lfx = gen->lpdf(x, gen->params); */
/*   if (c != 0.) { */
/*     s = (c > 0.) ? 1. : -1.; */
/*     Tfx = s * exp(c * lfx); */
/*   } */
    
/*   return (Tfx); */
/* } /\* end of Tf() *\/ */

/*---------------------------------------------------------------------------*/

int
Tfdd (double x, TINFLEX_GEN *gen, TINFLEX_IV *iv)
/*---------------------------------------------------------------------------*/
/* Compute transformed density and its derivatives                           */
/* by means of its log-density.                                              */
/*---------------------------------------------------------------------------*/
/*   x    ... argument                                                      */
/*   gen  ... Tinflex generator object                                      */
/*   iv   ... interval where the result is stored                            */
/*---------------------------------------------------------------------------*/
{
  double lfx, dlfx, d2lfx;
  double Tfx;
  double s;
  double c = iv->c;

  /* evaluate log density and its derivatives. */
  lfx = gen->lpdf(x, gen->params);
  dlfx = gen->dlpdf(x, gen->params);
  d2lfx = gen->d2lpdf(x, gen->params);

  if (c == 0.) {
    /* Case: T(x) = log(x) */
    iv->Tfx = lfx;
    iv->dTfx = dlfx;
    iv->d2Tfx = d2lfx;
  }

  else {
    /* Case: T(x) = sign(c) * x^c */
    s = (c > 0.) ? 1. : -1.;
    Tfx =  s * exp(c * lfx);
    iv->Tfx = Tfx;
    iv->dTfx = c * Tfx * dlfx;
    iv->d2Tfx = c * Tfx * (c * dlfx*dlfx + d2lfx);
  }

  /* success */
  return 0;
} /* end of Tfdd() */

/*---------------------------------------------------------------------------*/

double
Tinv (double c, double x)
/*---------------------------------------------------------------------------*/
/* Compute inverse transformation.                                           */
/*---------------------------------------------------------------------------*/
/*   c ... parameter for transformation                                      */
/*   x ... argument                                                          */
/*---------------------------------------------------------------------------*/
{
  if (c == 0.)
    /* Case: T(x) = log(x) */
    return (exp(x));

  if (c == -0.5)
    /* Case: T(x) = -1/sqrt(x) */
    return (1/(x*x));
  
  if (c == 1.)
    /* Case: T(x) = -1/x */
    return (x);

  else {
    /* Case: T(x) = sign(c) * x^c */
    double s = (c < 0.) ? -1. : 1; 
    return (R_pow(s*x, 1/c));
  }
} /* end of Tinv() */

/*---------------------------------------------------------------------------*/

double
FT (double c, double x)
/*---------------------------------------------------------------------------*/
/* Compute antiderivative of inverse transformation.                         */
/*---------------------------------------------------------------------------*/
/*   c ... parameter for transformation                                      */
/*   x ... argument                                                          */
/*---------------------------------------------------------------------------*/
{
  if (c == 0.)
    /* Case: T(x) = log(x) */
    return (exp(x));

  if (c == -0.5)
    /* Case: T(x) = -1/sqrt(x) */
    return (-1./x);
  
  if (c == -1.)
    /* Case: T(x) = -1/x */
    return (-log(-x));

  else {
    /* Case: T(x) = sign(c) * x^c */
    double s = (c < 0.) ? -1. : 1; 
    double xs = (s*x > 0.) ? s*x : 0.;
    return (s * c/(c+1) * R_pow(xs, (c+1.)/c));
    /* Remark: By construction s*x should always be non-negative. */
    /*         Thus if s*x < 0 we assume round-off errors.         */
  }
} /* end of FT() */

/*---------------------------------------------------------------------------*/

double
FTinv (double c, double x)
/*---------------------------------------------------------------------------*/
/* Compute inverse of antiderivative of inverse transformation.              */
/*---------------------------------------------------------------------------*/
/*   c ... parameter for transformation                                      */
/*   x ... argument                                                          */
/*---------------------------------------------------------------------------*/
{
  if (c == 0.)
    /* Case: T(x) = log(x) */
    return (log(x));

  if (c == -0.5)
    /* Case: T(x) = -1/sqrt(x) */
    return (-1./x);

  if (c == -1.)
    /* Case: T(x) = -1/x */
    return (-exp(-x));

  else {
    /* Case: T(x) = sign(c) * x^c */
    double s = (c < 0.) ? -1. : 1; 
    return (s * R_pow((c+1)/c * s * x, c/(c+1)));
  }
} /* end of FTinv() */

/*****************************************************************************/
/*                                                                           */
/* Auxiliary Functions.                                                      */
/*                                                                           */
/*****************************************************************************/

#define ARCMEAN_HARMONIC   (1.e3)   /* use harmonic mean when abs larger than this value */
#define ARCMEAN_ARITHMETIC (1.e-6)  /* use arithmetic mean when abs larger than this value */

double
arcmean( double x0, double x1 )
/*---------------------------------------------------------------------------*/
/* compute "arctan mean" of two numbers.                                     */
/*                                                                           */
/* parameters:                                                               */
/*   x0, x1 ... two numbers                                                  */
/*                                                                           */
/* return:                                                                   */
/*   mean                                                                    */
/*                                                                           */
/* comment:                                                                  */
/*   "arctan mean" = tan(0.5*(arctan(x0)+arctan(x1)))                        */
/*                                                                           */
/*   a combination of arithmetical mean (for x0 and x1 close to 0)           */
/*   and the harmonic mean (for |x0| and |x1| large).                        */
/*---------------------------------------------------------------------------*/
{
  double a0,a1;
  double r;

  /* we need x0 < x1 */
  if (x0>x1) {double tmp = x0; x0=x1; x1=tmp;}

  if (x1 < -ARCMEAN_HARMONIC || x0 > ARCMEAN_HARMONIC)
    /* use harmonic mean */
    return (2./(1./x0 + 1./x1));

  a0 = (x0<=R_NegInf) ? -M_PI/2. : atan(x0);
  a1 = (x1>=R_PosInf) ?  M_PI/2. : atan(x1);

  if (fabs(a0-a1) < ARCMEAN_ARITHMETIC)
    /* use arithmetic mean */
    r = 0.5*x0 + 0.5*x1;

  else
    /* use "arc mean" */
    r = tan((a0 + a1)/2.);

  return r;
} /* end of arcmean() */

/*---------------------------------------------------------------------------*/

/*****************************************************************************/
/*                                                                           */
/* Parameters for hat and squeeze in an interval.                            */
/*                                                                           */
/*****************************************************************************/

/* Remark:                                                                   */
/* Tangents and secants have the form                                        */
/*                                                                           */
/*    t(x) = a + b*(x-y)                                                     */
/*                                                                           */
/* where                                                                     */
/*    a ... intercept for hat in transformed scale                           */
/*    b ... slope for hat in transformed scale                               */
/*    y ... anchor point of linear function                                  */
/*                                                                           */
/*  If 'y' is either the left boundary or the right boundary.                */
/*  If possible then 'y' is chosen as the boundary point where the           */
/*  transformed density obtains higher valuation.                            */

/*---------------------------------------------------------------------------*/
/* Codes for types of interval (see paper):                                  */

enum {
      type_Ia   = -1,
      type_Ib   =  1,
      type_IIb  =  2,
      type_IIa  = -2,
      type_IIIa = -3,
      type_IIIb =  3,
      type_IVa  = -4,  /* concave */
      type_IVb  =  4   /* convex  */
};

/*---------------------------------------------------------------------------*/

int
type_iv (TINFLEX_IV *left, TINFLEX_IV *right)
/*---------------------------------------------------------------------------*/
/* Determine type of interval.                                               */
/*---------------------------------------------------------------------------*/
/*   left   ... data for boundary point to the left                          */
/*   right  ... data for boundary point to the right                         */
/*---------------------------------------------------------------------------*/
/* Return: type of interval;                                                 */
/*         0 if type cannot be determined.                                   */
/*---------------------------------------------------------------------------*/
{
  double c = left->c;   /* parameter 'c' for interval. */
  double R;             /* slope of secant */
  
  /* Case: Unbounded domain. */
  if ( ! R_FINITE(left->x)) {
    if (right->d2Tfx < 0. && right->dTfx > 0.)
      return (type_IVa);
    else
      return (0);
  }
  
  if ( ! R_FINITE(right->x)) {
    if (left->d2Tfx < 0. && left->dTfx < 0.)
      return (type_IVa);
    else
      return (0);
  }
  
  /* Case: Interval where the density vanishes at boundary */
  /*       (and thus the log-density is -Inf). */
  if ( (c > 0.  && left->Tfx == 0.)    ||
       (c <= 0. && left->Tfx <= R_NegInf) ) {
    if (right->d2Tfx < 0. && right->dTfx > 0.)
      return (type_IVa);
    if (right->d2Tfx > 0. && right->dTfx > 0.)
      return (type_IVb);
    else
      return (0);
  }
  
  if ((c > 0.  && right->Tfx == 0.)    ||
      (c <= 0. && right->Tfx <= R_NegInf) ) {
    if (left->d2Tfx < 0. && left->dTfx < 0.)
      return (type_IVa);
    if (left->d2Tfx > 0. && left->dTfx < 0.)
      return (type_IVb);
    else
      return (0);
  }
  
  /* Case: Domains where the density has a pole at boundary */
  /*       (and thus the transformed density equals 0 for c<0). */
  if (c < 0.) {
    if ((left->Tfx == 0.  && right->d2Tfx > 0.) ||
	(right->Tfx == 0. && left->d2Tfx > 0.)  )
      return (type_IVb);
  }
  
  /* Compute slope of secant. */
  R = (right->Tfx - left->Tfx) / (right->x - left->x);
  
  /* Check for all other possible cases. */
  if (left->dTfx >= R && right->dTfx >= R)
    return (type_Ia);
  
  if (left->dTfx <= R && right->dTfx <= R)
    return (type_Ib);
  
  if (left->d2Tfx < 0. && right->d2Tfx < 0.)
    return (type_IVa);
  
  if (left->d2Tfx > 0. && right->d2Tfx > 0.)
    return (type_IVb);
  
  if (left->d2Tfx < 0. || right->d2Tfx > 0.) {
    if (left->dTfx >= R && R >= right->dTfx)
      return (type_IIa);
  
    if (left->dTfx <= R && R <= right->dTfx)
      return (type_IIIa);
  }
  
  if (left->d2Tfx > 0. || right->d2Tfx < 0.) {
    if (left->dTfx >= R && R >= right->dTfx)
      return (type_IIb);
  
    if (left->dTfx <= R && R <= right->dTfx)
      return (type_IIIb);
  }
  
  /* Cannot estimate type of interval. */
  /* Do we need a warning? Probably not. */
  /* --> warning("cannot detect type of interval") */
  
  return (0);

} /* end of type_iv() */

/*---------------------------------------------------------------------------*/

double
area (double c, double a, double b, double y, double from, double to)
/*---------------------------------------------------------------------------*/
/* Compute area underneath hat or squeeze in particular interval.            */
/*---------------------------------------------------------------------------*/
/*   c        ... parameter for transformation                               */
/*   a, b, y  ... intercept, slope and anchor point of transformed hat or squeeze */
/*   from, to ... range of integration                                       */
/*---------------------------------------------------------------------------*/
/* Return: area;                                                             */
/*        0 in case of an error.                                             */
/*---------------------------------------------------------------------------*/
{
  double area;
  
  /* Compute area ... */
  area = do_area(c,a,b,y,from,to);
  
  /* ... and check result. */
  
  if (! R_FINITE(area)) {
    /* We return Inf in all cases where 'area' is not finite (e.g. NaN)      */
    area = R_PosInf;
  }
  
  if (area < 0.) {
    /* Area is strictly negative.                                            */
    /* Two possible reasons:                                                 */
    /* 1. The correct value of 'area' is extremely small or 0.               */
    /*    The negative number might result from cancelation.                 */
    /*    Then 'area' can be set to 0.                                       */
    /* 2. Due to sever round-off error, the result is numerical garbage.     */
    /*                                                                       */
    /* We observed (2) in our experiments. So we discard the result and      */
    /* return 'Inf' (which means that the interval must be split).           */

    area = R_PosInf;
  }

  return (area);
  
} /* end of area() */

/*...........................................................................*/

double
do_area (double c, double a, double b, double y, double from, double to)
/*---------------------------------------------------------------------------*/
/* Perform computation of area underneath hat or squeeze.                    */
/*---------------------------------------------------------------------------*/
/*   c        ... parameter for transformation                               */
/*   a, b, y  ... intercept, slope and anchor point of transformed hat or squeeze */
/*   from, to ... range of integration                                       */
/*---------------------------------------------------------------------------*/
/* Return: area.                                                             */
/*---------------------------------------------------------------------------*/
{
  double area;
  double s, sc, z;
  
  /* check for a new "empty" interval without data */
  if (ISNA(a)) {
    /* We have a new interval without data. */
    return (R_PosInf);
  }

  /* Test where the tangent is constructed: */
  /* s = +1 if tangent is constructed on lower boundary of interval; */
  /*     -1 if tangent is constructed on upper boundary of interval. */
  s = ((to-y) > (y-from)) ? 1. : -1.;

  /* Generally we have */
  /*    area <- (FT(c, a+b*(to-y)) - FT(c, a+b*(from-y))) / b */
  /* For numerical reasons we have to distinguish */
  /* between different values of 'c'. */

  if ( c == 0. ) {
    /* Case: T(x)=log(x) */

    z = s * b*(to-from);
    if ( fabs(z) > 1.e-6) {
      area = (exp(a+b*(to-y)) - exp(a+b*(from-y))) / b;
    } else {
      /* We need approximation by Taylor polynomial to avoid */
      /* severe round-off errors. */
      area = exp(a) * (to-from) * (1. + z/2. + z*z/6.);
    }
    return (area);
  }
  
  /* else: c!=0 */

  /* sign of c */
  sc = (c > 0.) ? 1. : -1.;
    
  /* The tangent to the transformed density must not vanish. */
  /* Otherwise, we cannot construct the hat function. */
  if (! (sc * (a+b*(from-y)) >= 0.) ||
      ! (sc * (a+b*(to-y))   >= 0.) ) {
    /* Case: numerical errors */
    if (!R_FINITE(a) && a < 0.) {
      /* underflow in computing Tf(x) (resulting in -Inf) */
      /* we assume that the area is 0 */
      return (0.);
    }
    /* else: */
    if (a < 1.e250 && ! R_FINITE(b)) {
      /* close to underflow when computing Tf(x) and overflow of tangent */
      /* we assume that the area is 0 */
      return (0.);
    }
    /* else: we have to split the interval */
    return (R_PosInf);
  }

  /* Transform b */
  z = s * b/a * (to-from);
  
  if (c == -0.5) {
    /* Case: T(x) = -1/sqrt(x) */

    if (fabs(z) > 0.5) {
      area = (-1./(a+b*(to-y)) + 1./(a+b*(from-y))) / b;
    } else {
      area = 1./(a*a) * (to-from) / (1. + z);
    }
    return (area);
  }

  if (c == -1.) {
    /* Case: T(x) = -1/x */

    if (fabs(z) > 1.e-6) {
      area = (-log(-a-b*(to-y)) + log(-a-b*(from-y))) / b;
    } else {
      /* Approximation by Taylor polynomial */
      area = -1./a * (to-from) * (1. - z/2. + z*z/3.);
    }
    return (area);
  }

  if (c == 1.) {
    /* Case: T(x) = x */

    area = 0.5 * a * (to-from) * (z+2.);
    return (area);
  }

  /* Case: T(x) = sgn(c) * x^c */
  /* For all other cases we only use a rough approximation in */
  /* case of numerical errors. */
  
  if (fabs(b)>1.e-10) {
    area = (FT(c, a+b*(to-y)) - FT(c, a+b*(from-y))) / b;
  } else {
    area = Tinv(c, a) * (to-from);
  }
  
  return (area);

} /* end of do_area() */

/*---------------------------------------------------------------------------*/

int
hat_iv (TINFLEX_IV *left, TINFLEX_IV *right, int link)
/*---------------------------------------------------------------------------*/
/* Compute hat and squeeze for a paricular interval.                         */
/*---------------------------------------------------------------------------*/
/*   left   ... data for boundary point to the left and for entire interval  */
/*   right  ... data for boundary point to the right                         */
/*   link   ... "pointer" to next interval to the right                      */
/*---------------------------------------------------------------------------*/
/* Return: vector with parameter of same kind as 'left'.                     */
/*---------------------------------------------------------------------------*/
{
  int type;       /* type of interval */
  double tl[3], tr[3], sc[3]; /* tangents and secants */
  double R, A_sq;

  /* macro for storing data in array */
#define setarray(a, x,y,z) {(a)[0]=(x); (a)[1]=(y); (a)[2]=(z);}
#define sethat(a) {left->ht_a = a[0]; left->ht_b = a[1]; left->ht_y = a[2];}
#define setsqueeze(a) {left->sq_a = a[0]; left->sq_b = a[1]; left->sq_y = a[2];}
#define setsqueezeNA() {left->sq_a = NA_REAL; left->sq_b = NA_REAL; left->sq_y = NA_REAL;}

  /* Insert pointer to next interval to the right */
  /* (next element in poor man's linked list). */
  left->next = link;

  /* Check for interval of length 0. */
  if ( left->x == right->x ) {
    setarray(tl, left->Tfx, 0., left->x);
    sethat(tl);
    setsqueeze(tl);
    left->A_ht = 0.;
    left->A_sq = 0.;
    left->type = 0;
    return (0);
  }
  
  /* Get type of distribution. */
  left->type = type = type_iv(left, right);

  /* Compute tangent at boundary points. */
  setarray(tl, left->Tfx, left->dTfx, left->x);
  setarray(tr, right->Tfx, right->dTfx, right->x);
  
  /* Compute secant. */
  R = (right->Tfx - left->Tfx) / (right->x - left->x);
  if (left->Tfx >= right->Tfx) { 
    setarray(sc, left->Tfx, R, left->x);
  }
  else {
    setarray(sc, right->Tfx, R, right->x);
  }
  
  /* Case: unbounded domains. */
  if (! R_FINITE(left->x) && type == type_IVa) {
    sethat(tr);
    setsqueezeNA();
  }
  else if (! R_FINITE(right->x) && type == type_IVa) {
    sethat(tl);
    setsqueezeNA();
  }

  /* Case: bounded domains. */
  else if (type == type_Ia) {
    sethat(tl);
    setsqueeze(tr);
  }
  else if (type == type_Ib) {
    sethat(tr);
    setsqueeze(tl);
  }
  else if (type == type_IIa) {
    sethat(tl);
    setsqueeze(sc);
  }
  else if (type == type_IIb) {
    sethat(tr);
    setsqueeze(sc);
  }
  else if (type == type_IIIa) {
    sethat(sc);
    setsqueeze(tr);
  }
  else if (type == type_IIIb) {
    sethat(sc);
    setsqueeze(tl);
  }
  else if (type == type_IVb) {
    sethat(sc);
    if (left->Tfx > right->Tfx) {
      setsqueeze(tr);
    }
    else {
      setsqueeze(tl);
    }
  }
   else if (type == type_IVa) {
     if (left->Tfx > right->Tfx) {
       sethat(tl);
     }
     else {
       sethat(tr);
     }
     setsqueeze(sc);
   }
    
  /* Compute area below hat. */
  left->A_ht = area(left->c, left->ht_a, left->ht_b, left->ht_y, left->x, right->x);

  /* Compute area below squeeze. */
  A_sq = area(left->c, left->sq_a, left->sq_b, left->sq_y, left->x, right->x);
  left->A_sq = (R_FINITE(A_sq)) ? A_sq : 0.;

  /* Return success */
  return (0);

#undef setarray
#undef sethat
#undef setsqueeze
#undef setsqueezeNA

} /* end of hat_iv() */

/*---------------------------------------------------------------------------*/

TINFLEX_GEN *
Tinflex_lib_setup (TINFLEX_FUNCT *lpdf, TINFLEX_FUNCT *dlpdf, TINFLEX_FUNCT *d2lpdf,
		   const void *params,
		   int n_ib, const double *ib, int n_c, const double *c,
		   double rho, int max_intervals)
/*---------------------------------------------------------------------------*/
/* Setup: compute hat and squeeze for density.                               */
/*---------------------------------------------------------------------------*/
/*   lpdf   ... log-density                                                  */
/*   dlpdf  ... derivative of log-density                                    */
/*   d2lpdf ... 2nd derivative of log-density                                */
/*   params ... parameters for log-density                                   */
/*   n_ib   ... length of ib                                                 */
/*   ib     ... interval boundaries of decomposition                         */
/*   n_c    ... length of c                                                  */
/*   c      ... parameter for transformation (global or for each interval)   */
/*   rho    ... performance parameter: requested upper bound for ratio       */
/*              between area below hat and area below squeeze                */
/*   max_intervals ... maximal numbers of intervals                          */
/*---------------------------------------------------------------------------*/
/* Return: Tfinflex generator object.                                        */
/*---------------------------------------------------------------------------*/
{
  int i;                       /* loop variable */
  double *cvec;                /* temporary vector of values for paramater 'c' */
  TINFLEX_IV *ivs;             /* pointer to list of all intervals */
  TINFLEX_IV *iv;              /* pointer to current and next interval */
  TINFLEX_IV *left, *right;    /* pointer to left (current) and right (next) interval */
  TINFLEX_GEN *gen;            /* generator object */
  int n_ivs;                   /* number of intervals  */
  int last_idx;                /* index of last interval used in linked list */ 
  int last_idx_save;           /* ... before starting current update loop */
  double A_ht_tot, A_sq_tot;   /* total area below hat / squeeze */
  double threshold;            /* threshold for splitting interval */
  double p;                    /* splitting point for interval */

  /* --- Check arguments */

  if (! (rho >= 1.0001)) {
    warning("Tinflex_setup(): argument 'rho' too small or invalid, using default");
    rho = 1.1;
  }

  if (! (max_intervals > 11)) {
    warning("Tinflex_setup(): argument 'max_intervals' too small or invalid, using default");
    max_intervals = 1001;
  }

  if (n_ib < 2 || n_ib > max_intervals / 2) {
    error("Tinflex_setup(): argument 'ib' invalid");
  }
  
  /* The boundaries must be sorted. */
  for (i=1; i<n_ib; i++) {
    if (ib[i-1] > ib[i]) {
      error("Tinflex_setup(): array 'ib' not sorted");
    }
  }

  /* Check parameters for transformation. */
  if (! (n_c == 1 || n_c == n_ib-1) ) {
    error("Tinflex_setup(): argument 'c' invalid: its length must equal either 1 or number of intervals");
  }

  if ((! R_FINITE(ib[0]) && ! (c[0] > -1.)) ||
      (! R_FINITE(ib[n_ib-1]) && ! (c[n_c-1] > -1.)) ) {
      error("Tinflex_setup(): (first and last) entry of argument 'c' must be greater than -1 for unbounded domains");
  }
  
  /* ........................................................................ */

  /* Get parameter 'c' for transformation. */
  cvec = R_Calloc(n_ib, double);
  if (n_c ==1 ) {
    for (i=0; i< n_ib; i++) {
      cvec[i] = c[0];
    }
  }
  else {
    memcpy(cvec, c, n_c * sizeof(double));
    /* duplicate last element. */
    cvec[n_ib-1] = c[n_c-1];
  }

  /* Create a table for hat and squeeze. */
  ivs = R_Calloc(max_intervals+1, TINFLEX_IV);

  /* Create generator object */
  gen = R_Calloc(1, TINFLEX_GEN);
  gen->c = cvec;
  gen->lpdf = lpdf;
  gen->dlpdf = dlpdf;
  gen->d2lpdf = d2lpdf;
  gen->params = params;
  gen->ivs = ivs;
  gen->n_ivs = 0;
  gen->Acum = NULL;
  gen->gt = NULL;

  /* We need transformed densities and their derivatives on */
  /* both boundaries of each interval. */
  /* However, the boundary to the right of one interval is the */
  /* boundary to the left of the consecutive interval. */
  /* If their parameters 'c' coincide we only need to store the */
  /* transformed density for the boundary to the left. */
  /* If these parameters differ we insert an interval of length 0 */
  /* in order to store the transformed density for the boundary */
  /* to the right.  */

  last_idx = -1L;
  for (i=0; i<n_ib; i++) {
    
    /* store data for interval */
    ++last_idx;
    iv = ivs + last_idx;
    iv->x = ib[i];
    iv->c = cvec[i];
    Tfdd (ib[i], gen, iv);

    /* Parameter 'c' changes for next interval. */
    /* Thus we insert an interval of length 0. */
    if (i+1 < n_ib && cvec[i] != cvec[i+1]) {
      ++last_idx;
      iv = ivs + last_idx;
      iv->x = ib[i+1];
      iv->c = cvec[i];
      Tfdd (ib[i+1], gen, iv);
    }
  }

  /* Terminate poor man's linked list. */
  (ivs+last_idx)->next = -1;

  /* remove temporary vector */
  R_Free(cvec);

  /* Compute parameters for hat and squeeze for initial intervals. */
  for (i=0; i<last_idx; i++) {
    hat_iv(/*left=*/ivs+i, /*right=*/ivs+i+1, /*link=*/i+1); 
#ifdef DEBUG
    debug_iv (ivs+i);
#endif
  }

#ifdef DEBUG
    debug_iv (ivs+last_idx);
#endif

  /* Compute total areas for initial hat and squeeze. */
  A_ht_tot = 0.;
  A_sq_tot = 0.;
  for (i=0; i<last_idx; i++) {
    A_ht_tot += (ivs+i)->A_ht;
    A_sq_tot += (ivs+i)->A_sq;
  }

  /* We have to split intervals where the area */
  /* between hat and squueze is too large. */
  while (! (A_ht_tot <= rho * A_sq_tot) && (last_idx < max_intervals)) {

    /* Compute average area. */
    threshold = 0.99 * (A_ht_tot - A_sq_tot) / last_idx;

    /* Split all intervals where area between hat and squeeze is too large. */
    last_idx_save = last_idx;
    for (i=0; i<=last_idx_save; i++) {

      /* scan through all intervals */
      iv = ivs+i;

      if ( iv->next < 0 || (iv->A_ht  - iv->A_sq) < threshold) {
	continue;
      }
      if ( last_idx >= max_intervals) {
	warning("Tinflex_setup(): maximum number of intervals reached");
	break;
      }
      
      /* split intervals where area between hat and squeeze is too large. */
      left = iv;
      right = ivs + iv->next;

      /* Splitting point for interval (use "arc mean"). */
      p = arcmean(iv->x, right->x);

      /* Create new interval. */
      ++last_idx;

      iv = ivs + last_idx;
      iv->x = p;
      iv->c = left->c;
      Tfdd (p, gen, iv);

      /* New interval on r.h.s. */
      hat_iv(/*left=*/iv, /*right=*/right, /*link=*/left->next); 
      /* Update interval on l.h.s. */
      hat_iv(/*left=*/left, /*right=*/iv, /*link=*/last_idx); 
    }

    /* Update total areas. */
    A_ht_tot = 0.;
    A_sq_tot = 0.;
    for (i=0; i<=last_idx; i++) {
      iv = ivs+i;
      if (iv->next >= 0L) {
	A_ht_tot += iv->A_ht;
	A_sq_tot += iv->A_sq;
      }
    }
  }

  /* Check result. */
  if (! R_FINITE(A_ht_tot)) { 
    Tinflex_lib_free(gen);
    error("Tinflex_setup(): Cannot create hat function. A_hat is not finite: A_hat=%g", A_ht_tot);
  }
    
  if (! (A_ht_tot > 0.)) {
    Tinflex_lib_free(gen);
    error("Tinflex_setup(): Cannot create hat function. A_hat is not > 0: A_hat=%g", A_ht_tot);
  }

  if (! R_FINITE(A_sq_tot)) {
    Tinflex_lib_free(gen);
    error("Tinflex_setup(): Cannot create squeeze. A_squeeze is not finite: A_squeeze=%g", A_sq_tot);
  }
  
  if (! R_FINITE(A_sq_tot >= 0.)) {
    Tinflex_lib_free(gen);
    error("Tinflex_setup(): Cannot create squeeze. A_squeeze is not >= 0: A_squeeze=%g", A_sq_tot);
  }

  if (! (A_ht_tot >= A_sq_tot)) {
    Tinflex_lib_free(gen);
    error("Tinflex_setup(): Invalid input: A_hat = %g < %g = A_squeeze!", A_ht_tot, A_sq_tot);
  }

  if (! (A_ht_tot <= rho * A_sq_tot)) {
    error("Tinflex_setup(): ratio A_hat / A_squeeze = %g is larger than rho = %g",
	  A_ht_tot / A_sq_tot, rho);
  }

  /* store areas */
  gen->A_ht_tot = A_ht_tot;
  gen->A_sq_tot = A_sq_tot;
  
  /* Number of intervals */
  n_ivs = last_idx;
  
  /* Truncate working array. */
  ivs = R_Realloc(ivs, n_ivs+1L, TINFLEX_IV);

  /* Sort list (so that the 'next' marker becomes obsolete). */
  qsort(ivs, (n_ivs+1L), sizeof(TINFLEX_IV), cmpiv);
  
  /* Store intervals in generator object */
  gen->ivs = ivs;
  gen->n_ivs = n_ivs;

  /* Compute guide table. */
  make_guide_table(gen);

  /* Return generator object.  */
  return (gen);
  
} /* end of Tinflex_lib_setup() */

/*...........................................................................*/

int cmpiv (const void *p1, const void *p2)
/* auxiliary function for sorting intervals w.r.t. is 'x' value */
{
  const TINFLEX_IV *iv1 = p1;
  const TINFLEX_IV *iv2 = p2;

  if (iv1->next < 0L)        /* iv1 is terminating interval */
    return (1L);
  else if (iv2->next < 0L)   /* iv2 is terminating interval */
    return (-1L);
  else if (iv1->x < iv2->x) 
    return (-1L);
  else if (iv1->x > iv2->x) 
    return (1L);
  else
    return (0L);
} /* end of cmpiv() */

/*---------------------------------------------------------------------------*/

void * Tinflex_lib_free (TINFLEX_GEN *gen)
/*---------------------------------------------------------------------------*/
/* Free allocated memory.                                                    */
/*---------------------------------------------------------------------------*/
/*   gen  ... Tinflex generator object                                       */
/*---------------------------------------------------------------------------*/
{
  R_Free(gen->ivs);
  R_Free(gen->Acum);
  R_Free(gen->gt);
  R_Free(gen);

  return (NULL);
} /* end of Tinflex_lib_free() */

/*---------------------------------------------------------------------------*/

int 
make_guide_table (TINFLEX_GEN *gen)
/*---------------------------------------------------------------------------*/
/* Create guide table for drawing interval at random.                        */
/* The result is stored in 'Acum' and 'gt'.                                  */
/*---------------------------------------------------------------------------*/
/*   gen  ... Tinflex generator object                                       */
/*---------------------------------------------------------------------------*/
{
  TINFLEX_IV *ivs;      /* data for hat and squeeze for each interval  */
  int n_ivs;            /* number of intervals                         */
  double *Acum;         /* array for storing cumulative areas          */
  int *gt;              /* array for storing guide table               */
  double sum, Astep;    /* variable for computing guide table          */
  int i,j;              /* loop variables                              */

  /* Get interval data. */
  ivs = gen->ivs;
  n_ivs = gen->n_ivs;

  /* allocate memory */
  gen->Acum = Acum = R_Calloc(n_ivs, double);
  gen->gt = gt = R_Calloc(n_ivs, int);

  /* Compute cumulated areas (probabilities). */
  for ( i=0, sum = 0.; i <n_ivs; i++ ) {
    Acum[i] = (sum += (ivs + i)->A_ht);
  }
  
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

  /* store total area */
  gen->A_ht_tot = Acum[n_ivs-1];

  /* success */
  return (0);
  
} /* end of make_guide_table() */

/*---------------------------------------------------------------------------*/

SEXP
Tinflex_lib_sample (TINFLEX_GEN *gen, int n)
/*---------------------------------------------------------------------------*/
/* Draw sample from Tinflex generator object.                                */
/*---------------------------------------------------------------------------*/
/*   gen  ... Tinflex generator object                                       */
/*   n    ... sample size                                                    */
/*---------------------------------------------------------------------------*/
{
  TINFLEX_IV *ivs;            /* data for hat and squeeze for all intervals  */
  int n_ivs;                  /* number of intervals                         */
  double *Acum;               /* array for storing cumulative areas          */
  double A_ht_tot;            /* total area below hat                        */
  int *gt;                    /* guide table                                 */
  SEXP sexp_res = R_NilValue; /* R object for storing random sample          */
  int i;                      /* auxiliary loop variable                     */

  /* Acceptance-rejection loop: */
  double U,V;                 /* uniform random numbers                      */
  int J;                      /* index of random interval                    */
  TINFLEX_IV *iv;             /* data for particular interval                */
  double a,b,y,x0,c;          /* parameters for hat (tangent)                */
  double z;                   /* auxiliary variable                          */
  double t;                   /* hat (tangent) at construction point         */
  double X;                   /* random point                                */
  double hx, sx;              /* hat and squeeze at X                        */

  /* check arguments */
  if (n<=0) {
    error("Tinflex_sample(): sample size 'n' must be positive integer"); 
  }

  /* extract data */
  ivs = gen->ivs;           /* parameters for hat and squeeze */
  n_ivs = gen->n_ivs;       /* number of intervals */
  Acum = gen->Acum;         /* accumulated areas */
  A_ht_tot = gen->A_ht_tot; /* total area below hat                        */
  gt = gen->gt;             /* guide table                                 */
  
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
       * 	sum += (ivs + J)->A_ht;
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
      iv = ivs + J;
      
      a  = iv->ht_a;
      b  = iv->ht_b;
      y  = iv->ht_y;
      x0 = iv->x;
      c  = iv->c;

      /* Value of transformed hat (tangent) at left boundary of interval. */
      t = a+b*(x0-y);
      
      /* Need a uniform random number in interval (0, iv->A_ht). */
      /*    U = iv->A_ht * unif_rand();                          */
      /* However, we can recycle U for this task.                */
      U = iv->A_ht + U - Acum[J];
      
      /* Generate from "hat distribution":                  */
      /*    X = y + ( FTinv(cT, FT(cT,t) + b*U) - a ) / b;  */
      
      /* For numerical reasons we have to distinguish */
      /* between different values of 'cT'.            */
      if (c == 0.) {
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
      
      else if (c == -0.5) {
      	/* Case: T(x) = -1/sqrt(x) */
	
      	z = U * b * t;
      	if (fabs(z) > 1.e-6)
      	  X = y + ( 1./(1./t - b*U) - a ) / b;
      	else
      	  X = x0 + U * t*t * (1 + z + z*z);
      }
      
      else if (c == 1.) {
      	/* Case: T(x) = x */
      	z = U * b / (t*t);
      	if (fabs(z) > 1.e-6)
      	  X = y + (FTinv(c, FT(c, t)+b*U)-a)/b;
	/* FIXME: why FT() and FTinv() ? */
	/* X = y + (sqrt(t*t+b*U)-a)/b; */
      	else
      	  X = x0 + U / t * (1 - z/2 + z*z/2);
      }
      
      else {
      	/* Case: T(x)=sign(c)*x^c */
      	/* For all other cases we only use a rough approximation in */
      	/* case of numerical errors. */
      
      	if (fabs(b)>1e-10) {
      	  X = y + ( FTinv(c, FT(c,t) + b*U) - a ) / b;
	}
	else { 
       	  U = U / iv->A_ht;
      	  X = (1.-U) * x0 + U * (iv+1)->x;
       	}
      }
      
      /* Compute hat and squeeze at X. */
      hx = Tinv(c, a + b * (X-y));
      if (iv->A_sq > 0)
      	sx = Tinv(c, iv->sq_a + iv->sq_b * (X-iv->sq_y));
      else
      	sx = 0.;
      
      /* Accept or reject. */
      V = hx * unif_rand();
      
      if ((V <= sx) ||
	  (V <= exp(gen->lpdf(X, gen->params))) )
	break;
    }
    
    /* Store random point. */
    NUMERIC_POINTER(sexp_res)[i] = X;
  }

  /* Update state for the R built-in URNG. */
  PutRNGstate();

  /* Return random sample to R. */
  UNPROTECT(1);
  return sexp_res;
} /* end of Tinflex_lib_sample() */

/*---------------------------------------------------------------------------*/

#ifdef DEBUG
void
debug_iv (TINFLEX_IV *iv)
{
  Rprintf("x = %g\tc = %g\n", iv->x, iv->c);
  Rprintf("ht_a = %g\tht_b = %g\tht_y = %g\n", iv->ht_a, iv->ht_b, iv->ht_y);
  Rprintf("sq_a = %g\tsq_b = %g\tsq_y = %g\n", iv->sq_a, iv->sq_b, iv->sq_y);
  Rprintf("A_ht = %g\tA_sq = %g\n", iv->A_ht, iv->A_sq);
  Rprintf("Tfx = %g\tdTfx = %g\td2Tfx = %g\n", iv->Tfx, iv->dTfx, iv->d2Tfx);
  Rprintf("type = %d\t next = %d\n", iv->type, iv->next);
  Rprintf("---------------------------------------\n");
} /* end of debug_iv() */

#endif

/*---------------------------------------------------------------------------*/
