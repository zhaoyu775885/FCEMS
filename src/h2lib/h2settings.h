#ifndef H2SETTINGS_H
#define H2SETTINGS_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <complex>
#include <math.h>
#include <mkl.h>

using namespace std;

#define USE_BLAS 1
#define USE_COMPLEX 1

typedef complex<double> DComplex;
typedef DComplex field;
typedef DComplex *pfield;
typedef const DComplex *pcfield;
typedef unsigned uint;

extern const int i_zero;
extern const int i_one;

extern const double r_zero;
extern const double r_one;
extern const double r_minusone;

extern const field f_one;
extern const field f_zero;
extern const field f_minusone;

extern const field f_i;

extern const char *_h2_ntrans;
extern const char *_h2_trans;
extern const char *_h2_adj;
extern const char *_h2_left;
extern const char *_h2_right;
extern const char *_h2_lower;
extern const char *_h2_upper;
extern const char *_h2_unit;
extern const char *_h2_nonunit;
extern const char *_h2_vectors;
extern const char *_h2_skinnyvectors;
extern const char *_h2_novectors;

extern const CBLAS_LAYOUT h2_row_major;
extern const CBLAS_LAYOUT h2_col_major;
extern const CBLAS_SIDE h2_left;
extern const CBLAS_SIDE h2_right;
extern const CBLAS_UPLO h2_lower;
extern const CBLAS_UPLO h2_upper;
extern const CBLAS_TRANSPOSE h2_ntrans;
extern const CBLAS_TRANSPOSE h2_trans;
extern const CBLAS_TRANSPOSE h2_ctrans;
extern const CBLAS_DIAG h2_udiag;
extern const CBLAS_DIAG h2_ndiag;



typedef enum {
	AMATRIX = 0,
	HMATRIX = 1,
	H2MATRIX = 2
} matrixtype;

#endif
