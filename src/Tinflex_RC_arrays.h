/*
 * -- Automatically created from 'setup.R'
 * -- by '../devel/make_Tinflex_RC_arrays_h.R'
 * Fri Jul  5 11:41:12 2019 
 */

/* Number of data for each interval */
static const int L_IVS = 15 ; /* = length(iv.struct) */

/* Parameters for each interval */
enum {
	IV_X = 0,
	IV_C = 1,
	IV_HT_A = 2,
	IV_HT_B = 3,
	IV_HT_Y = 4,
	IV_SQ_A = 5,
	IV_SQ_B = 6,
	IV_SQ_Y = 7,
	IV_A_HT = 8,
	IV_A_SQ = 9,
	IV_TYPE = 10,
	IV_TFX = 11,
	IV_DTFX = 12,
	IV_D2TFX = 13,
	IV_NEXT = 14,
};

/* Names of list elements for generator */
enum {
	GEN_IVS = 0,
	GEN_LPDF = 1,
	GEN_A_HT_TOT = 2,
	GEN_A_SQ_TOT = 3,
	GEN_ENV = 4,
	GEN_INIV = 5,
	GEN_ACUM = 6,
	GEN_GT = 7,
};

