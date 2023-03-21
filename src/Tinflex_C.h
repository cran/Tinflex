/*****************************************************************************/
/*                                                                           */
/*  Tinflex_C                                                                */
/*  Setup and sampling routine for generator                                 */
/*  (C only version)                                                         */
/*                                                                           */
/*****************************************************************************/

SEXP Tinflex_C_setup (SEXP sexp_obj, SEXP sexp_env,
		      SEXP sexp_lpdf, SEXP sexp_dlpdf, SEXP sexp_d2lpdf,
		      SEXP sexp_ib, SEXP sexp_c,
		      SEXP sexp_rho, SEXP sexp_max_intervals);

SEXP Tinflex_C_sample (SEXP sexp_gen, SEXP sexp_n);

SEXP Tinflex_C_2_R (SEXP sexp_gen);

/*---------------------------------------------------------------------------*/
