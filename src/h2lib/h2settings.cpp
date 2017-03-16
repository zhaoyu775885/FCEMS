#include "h2settings.h"

const double r_zero = 0.0;
const double r_one = 1.0;
const double r_minusone = -1.0;

const field f_zero = 0.0;
const field f_one = 1.0;
const field f_minusone = -1.0;
const field f_i(0.0, 1.0);

const int i_zero = 0;
const int i_one = 1;

const char *_h2_ntrans = "Not transposed";
const char *_h2_trans = "Transposed";
const char *_h2_adj = "Conjugate transposed";
const char *_h2_left = "Left";
const char *_h2_right = "Right";
const char *_h2_lower = "Lower";
const char *_h2_upper = "Upper";
const char *_h2_unit = "Unit";
const char *_h2_nonunit = "Non-unit";
const char *_h2_vectors = "Vectors";
const char *_h2_skinnyvectors = "Skinny Vectors";
const char *_h2_novectors = "No vectors";

const CBLAS_LAYOUT h2_row_major = CblasRowMajor;
const CBLAS_LAYOUT h2_col_major = CblasColMajor;
const CBLAS_SIDE h2_left = CblasLeft;
const CBLAS_SIDE h2_right = CblasRight;
const CBLAS_UPLO h2_lower = CblasLower;
const CBLAS_UPLO h2_upper = CblasUpper;
const CBLAS_TRANSPOSE h2_ntrans = CblasNoTrans;
const CBLAS_TRANSPOSE h2_trans = CblasTrans;
const CBLAS_TRANSPOSE h2_ctrans = CblasConjTrans;
const CBLAS_DIAG h2_ndiag = CblasNonUnit;
const CBLAS_DIAG h2_udiag = CblasUnit;
