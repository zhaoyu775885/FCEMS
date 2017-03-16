/* ------------------------------------------------------------
 This is the file "factorizations.c" of the H2Lib package.
 All rights reserved, Steffen Boerm 2010
 ------------------------------------------------------------ */

#include "factorizations.h"

#include "h2settings.h"
#include "h2basic.h"

#include <math.h>
#include <stdio.h>

/* ------------------------------------------------------------
 Diagonal matrices
 ------------------------------------------------------------ */

void
diagsolve_amatrix_avector(bool atrans, pcamatrix a, pavector x)
{
  pcfield   aa = a->a;
  int      lda = a->ld;
  pfield    xv = x->v;
  int      n = UINT_MIN(a->rows, a->cols);
  int      i;

  if (atrans) {
    for (i = 0; i < n; i++)
      xv[i] /= CONJ(aa[i + i * lda]);
  }
  else {
    for (i = 0; i < n; i++)
      xv[i] /= aa[i + i * lda];
  }
}

#ifdef USE_BLAS
void
diagsolve_amatrix(bool atrans, pcamatrix a, bool xtrans, pamatrix x)
{
  int      n = UINT_MIN(a->rows, a->cols);
  int      lda = a->ld;
  int      ldx = x->ld;
  field    *aa = a->a;
  field    *xa = x->a;
  field     alpha;
  int      i;

  if (xtrans) {
    for (i = 0; i < n; i++) {
      alpha = (atrans ? 1.0 / CONJ(aa[i + i * lda]) : 1.0 / aa[i + i * lda]);
      h2_scal(x->rows, &alpha, xa + i * ldx, i_one);
    }
  }
  else {
    for (i = 0; i < n; i++) {
      alpha = (atrans ? 1.0 / CONJ(aa[i + i * lda]) : 1.0 / aa[i + i * lda]);
      h2_scal(x->cols, &alpha, xa + i, ldx);
    }
  }
}
#else
void
diagsolve_amatrix(bool atrans, pcamatrix a, bool xtrans, pamatrix x)
{
  int      n = UINT_MIN(a->rows, a->cols);
  int      lda = a->ld;
  int      ldx = x->ld;
  field    *aa = a->a;
  field    *xa = x->a;
  field     alpha;
  int      i, j;

  if (xtrans) {
    for (i = 0; i < n; i++) {
      alpha = (atrans ? 1.0 / CONJ(aa[i + i * lda]) : 1.0 / aa[i + i * lda]);
      for (j = 0; j < x->rows; j++)
	xa[j + i * ldx] *= alpha;
    }
  }
  else {
    for (i = 0; i < n; i++) {
      alpha = (atrans ? 1.0 / CONJ(aa[i + i * lda]) : 1.0 / aa[i + i * lda]);
      for (j = 0; j < x->cols; j++)
	xa[i + j * ldx] *= alpha;
    }
  }
}
#endif

void
diageval_amatrix_avector(bool atrans, pcamatrix a, pavector x)
{
  int      n = UINT_MIN(a->rows, a->cols);
  int      lda = a->ld;
  field    *aa = a->a;
  field    *xv = x->v;
  field     alpha;
  int      i;

  assert(x->dim >= a->cols);

  for (i = 0; i < n; i++) {
    alpha = (atrans ? CONJ(aa[i + i * lda]) : aa[i + i * lda]);
    xv[i] *= alpha;
  }
}

#ifdef USE_BLAS
void
diageval_amatrix(bool atrans, pcamatrix a, bool xtrans, pamatrix x)
{
  int      n = UINT_MIN(a->rows, a->cols);
  int      lda = a->ld;
  int      ldx = x->ld;
  field    *aa = a->a;
  field    *xa = x->a;
  field     alpha;
  int      i;

  if (xtrans) {
    for (i = 0; i < n; i++) {
      alpha = (atrans ? CONJ(aa[i + i * lda]) : aa[i + i * lda]);
      h2_scal(x->rows, &alpha, xa + i * ldx, i_one);
    }
  }
  else {
    for (i = 0; i < n; i++) {
      alpha = (atrans ? CONJ(aa[i + i * lda]) : aa[i + i * lda]);
      h2_scal(x->cols, &alpha, xa + i, ldx);
    }
  }
}
#else
void
diageval_amatrix(bool atrans, pcamatrix a, bool xtrans, pamatrix x)
{
  int      n = UINT_MIN(a->rows, a->cols);
  int      lda = a->ld;
  int      ldx = x->ld;
  field    *aa = a->a;
  field    *xa = x->a;
  field     alpha;
  int      i, j;

  if (xtrans) {
    for (i = 0; i < n; i++) {
      alpha = (atrans ? CONJ(aa[i + i * lda]) : aa[i + i * lda]);
      for (j = 0; j < x->rows; j++)
	xa[j + i * ldx] *= alpha;
    }
  }
  else {
    for (i = 0; i < n; i++) {
      alpha = (atrans ? CONJ(aa[i + i * lda]) : aa[i + i * lda]);
      for (j = 0; j < x->cols; j++)
	xa[i + j * ldx] *= alpha;
    }
  }
}
#endif

#ifdef USE_BLAS
void
diageval_realavector_amatrix(field alpha, bool atrans, pcrealavector a,
			     bool xtrans, pamatrix x)
{
  double*     av = a->v;
  pfield    xa = x->a;
  int ldx = x->ld;
  field     beta;
  int  j;

  if (xtrans) {
    assert(x->cols <= a->dim);

    if (a->dim < 1 || x->rows < 1) {	/* Quick exit */
      return;
    }

    for (j = 0; j < x->cols; j++) {
      beta = (atrans ? CONJ(alpha) * av[j] : CONJ(alpha * av[j]));
      h2_scal(x->rows, &beta, xa + j * ldx, i_one);
    }
  }
  else {
    assert(x->rows <= a->dim);

    if (a->dim < 1 || x->cols < 1) {	/* Quick exit */
      return;
    }

    for (j = 0; j < x->rows; j++) {
      beta = (atrans ? alpha * av[j] : alpha * av[j]);
      h2_scal(x->cols, &beta, xa + j, ldx);
    }
  }
}
#else
void
diageval_realavector_amatrix(field alpha,
			     bool atrans, pcrealavector a, bool xtrans,
			     pamatrix x)
{
  double*     av = a->v;
  pfield    xa = x->a;
  longindex ldx = x->ld;
  field     beta;
  int  i, j;

  if (xtrans) {
    assert(x->cols <= a->dim);

    if (a->dim < 1 || x->rows < 1) {	/* Quick exit */
      return;
    }

    for (i = 0; i < x->cols; i++) {
      beta = (atrans ? CONJ(alpha) * av[i] : CONJ(alpha * av[i]));
      for (j = 0; j < x->rows; j++) {
	xa[j + i * ldx] *= beta;
      }
    }
  }
  else {
    assert(x->rows <= a->dim);

    if (a->dim < 1 || x->cols < 1) {	/* Quick exit */
      return;
    }

    for (i = 0; i < x->rows; i++) {
      beta = (atrans ? alpha * av[i] : alpha * av[i]);
      for (j = 0; j < x->cols; j++) {
	xa[i + j * ldx] *= beta;
      }
    }
  }
}
#endif

/* ------------------------------------------------------------
 Triangular matrices
 ------------------------------------------------------------ */

#ifdef USE_BLAS
static void
lowersolve_amatrix_avector(bool aunit, bool atrans, pcamatrix a, pavector x)
{
  int      n = UINT_MIN(a->rows, a->cols);
  // revised by Yu

  assert(x->dim >= a->rows);
  assert(x->dim >= a->cols);

  if (atrans) {
    h2_trsm(h2_col_major, h2_left, h2_lower, h2_ctrans, (aunit ? h2_udiag : h2_ndiag),
	    n, i_one, (DComplex*)&f_one, a->a, a->ld, x->v, x->dim);
  }
  else {
    h2_trsm(h2_col_major, h2_left, h2_lower, h2_ntrans, (aunit ? h2_udiag : h2_ndiag),
	    n, i_one, (DComplex*)&f_one, a->a, a->ld, x->v, x->dim);
  }
}

static void
uppersolve_amatrix_avector(bool aunit, bool atrans, pcamatrix a, pavector x)
{
  int      n = UINT_MIN(a->rows, a->cols);

  assert(x->dim >= a->rows);
  assert(x->dim >= a->cols);

  if (atrans) {
    h2_trsm(h2_col_major, h2_left, h2_upper, h2_ctrans, (aunit ? h2_udiag : h2_ndiag),
	    n, i_one, &f_one, a->a, a->ld, x->v, x->dim);
  }
  else {
    h2_trsm(h2_col_major, h2_left, h2_upper, h2_ntrans, (aunit ? h2_udiag : h2_ndiag),
	    n, i_one, &f_one, a->a, a->ld, x->v, x->dim);
  }
}
#else
static void
lowersolve_amatrix_avector(bool aunit, bool atrans, pcamatrix a, pavector x)
{
  pcfield   aa = a->a;
  int      lda = a->ld;
  pfield    xv = x->v;
  int      n = UINT_MIN(a->rows, a->cols);
  field     newval;
  int      i, j;

  assert(x->dim >= a->rows);
  assert(x->dim >= a->cols);

  if (atrans) {
    for (j = n; j-- > 0;) {
      newval = (aunit ? xv[j] : (xv[j] /= CONJ(aa[j + j * lda])));
      for (i = 0; i < j; i++)
	xv[i] -= CONJ(aa[j + i * lda]) * newval;
    }
  }
  else {
    for (j = 0; j < n; j++) {
      newval = (aunit ? xv[j] : (xv[j] /= aa[j + j * lda]));
      for (i = j + 1; i < n; i++)
	xv[i] -= aa[i + j * lda] * newval;
    }
  }
}

static void
uppersolve_amatrix_avector(bool aunit, bool atrans, pcamatrix a, pavector x)
{
  pcfield   aa = a->a;
  int      lda = a->ld;
  pfield    xv = x->v;
  int      n = UINT_MIN(a->rows, a->cols);
  field     newval;
  int      i, j;

  assert(x->dim >= a->rows);
  assert(x->dim >= a->cols);

  if (atrans) {
    for (j = 0; j < n; j++) {
      newval = (aunit ? xv[j] : (xv[j] /= CONJ(aa[j + j * lda])));
      for (i = j + 1; i < n; i++)
	xv[i] -= CONJ(aa[j + i * lda]) * newval;
    }
  }
  else {
    for (j = n; j-- > 0;) {
      newval = (aunit ? xv[j] : (xv[j] /= aa[j + j * lda]));
      for (i = 0; i < j; i++)
	xv[i] -= aa[i + j * lda] * newval;
    }
  }
}
#endif

void
triangularsolve_amatrix_avector(bool alower, bool aunit, bool atrans,
				pcamatrix a, pavector x)
{
  if (alower)
    lowersolve_amatrix_avector(aunit, atrans, a, x);
  else
    uppersolve_amatrix_avector(aunit, atrans, a, x);
}

#ifdef USE_BLAS
static void
lowersolve_amatrix(bool aunit, bool atrans, pcamatrix a,
		   bool xtrans, pamatrix x)
{
  int      n = UINT_MIN(a->rows, a->cols);
  field    *aa = a->a;
  int      lda = a->ld;
  field    *xa = x->a;
  int      ldx = x->ld;

  if (atrans) {
    if (xtrans) {
      assert(x->cols >= n);

      h2_trsm(h2_col_major, h2_right, h2_lower, h2_ntrans,
	      (aunit ? h2_udiag : h2_ndiag), x->rows, n, &f_one, aa, lda, xa, ldx);
    }
    else {
      assert(x->rows >= n);

      h2_trsm(h2_col_major, h2_left, h2_lower, h2_ctrans, (aunit ? h2_udiag : h2_ndiag), n, x->cols, &f_one, aa, lda, xa, ldx);
    }
  }
  else {
    if (xtrans) {
      assert(x->cols >= n);

      h2_trsm(h2_col_major, h2_right, h2_lower, h2_ctrans, (aunit ? h2_udiag : h2_ndiag), x->rows, n, &f_one, aa, lda, xa, ldx);
    }
    else {
      assert(x->rows >= n);

      h2_trsm(h2_col_major, h2_left, h2_lower, h2_ntrans, (aunit ? h2_udiag : h2_ndiag), n, x->cols, &f_one, aa, lda, xa, ldx);
    }
  }
}

static void
uppersolve_amatrix(bool aunit, bool atrans, pcamatrix a,
		   bool xtrans, pamatrix x)
{
  int      n = UINT_MIN(a->rows, a->cols);
  field    *aa = a->a;
  int      lda = a->ld;
  field    *xa = x->a;
  int      ldx = x->ld;

  if (atrans) {
    if (xtrans) {
      assert(x->cols >= n);

      h2_trsm(h2_col_major, h2_right, h2_upper, h2_ntrans, (aunit ? h2_udiag : h2_ndiag), x->rows, n, &f_one, aa, lda, xa, ldx);
    }
    else {
      assert(x->rows >= n);

      h2_trsm(h2_col_major, h2_left, h2_upper
, h2_ctrans, (aunit ? h2_udiag : h2_ndiag),
	      n, x->cols, &f_one, aa, lda, xa, ldx);
    }
  }
  else {
    if (xtrans) {
      assert(x->cols >= n);

      h2_trsm(h2_col_major, h2_right, h2_upper
, h2_ctrans, (aunit ? h2_udiag : h2_ndiag),
	      x->rows, n, &f_one, aa, lda, xa, ldx);
    }
    else {
      assert(x->rows >= n);

      h2_trsm(h2_col_major, h2_left, h2_upper
, h2_ntrans,
	      (aunit ? h2_udiag : h2_ndiag), n, x->cols, &f_one, aa,
	      lda, xa, ldx);
    }
  }
}
#else
static void
lowersolve_amatrix(bool aunit, bool atrans, pcamatrix a,
		   bool xtrans, pamatrix x)
{
  int      n = UINT_MIN(a->rows, a->cols);
  int      lda = a->ld;
  int      ldx = x->ld;
  pfield    aa = a->a;
  pfield    xa = x->a;
  int      i, j, k;
  field     alpha;

  if (atrans) {
    if (xtrans) {
      assert(x->cols >= n);

      for (k = n; k-- > 0;) {
	if (!aunit) {
	  alpha = 1.0 / aa[k + k * lda];
	  for (i = 0; i < x->rows; i++)
	    xa[i + k * ldx] *= alpha;
	}
	for (i = 0; i < x->rows; i++)
	  for (j = 0; j < k; j++)
	    xa[i + j * ldx] -= xa[i + k * ldx] * aa[k + j * lda];
      }
    }
    else {
      assert(x->rows >= n);

      for (k = n; k-- > 0;) {
	if (!aunit) {
	  alpha = 1.0 / CONJ(aa[k + k * lda]);
	  for (j = 0; j < x->cols; j++)
	    xa[k + j * ldx] *= alpha;
	}
	for (i = 0; i < k; i++)
	  for (j = 0; j < x->cols; j++)
	    xa[i + j * ldx] -= CONJ(aa[k + i * lda]) * xa[k + j * ldx];
      }
    }
  }
  else {
    if (xtrans) {
      assert(x->cols >= n);

      for (k = 0; k < n; k++) {
	if (!aunit) {
	  alpha = 1.0 / CONJ(aa[k + k * lda]);
	  for (i = 0; i < x->rows; i++)
	    xa[i + k * ldx] *= alpha;
	}
	for (i = 0; i < x->rows; i++)
	  for (j = k + 1; j < n; j++)
	    xa[i + j * ldx] -= xa[i + k * ldx] * CONJ(aa[j + k * lda]);
      }
    }
    else {
      assert(x->rows >= n);

      for (k = 0; k < n; k++) {
	if (!aunit) {
	  alpha = 1.0 / aa[k + k * lda];
	  for (j = 0; j < x->cols; j++)
	    xa[k + j * ldx] *= alpha;
	}
	for (i = k + 1; i < n; i++)
	  for (j = 0; j < x->cols; j++)
	    xa[i + j * ldx] -= aa[i + k * lda] * xa[k + j * ldx];
      }
    }
  }
}

static void
uppersolve_amatrix(bool aunit, bool atrans, pcamatrix a,
		   bool xtrans, pamatrix x)
{
  int      n = UINT_MIN(a->rows, a->cols);
  int      lda = a->ld;
  int      ldx = x->ld;
  pfield    aa = a->a;
  pfield    xa = x->a;
  int      i, j, k;
  field     alpha;

  if (atrans) {
    if (xtrans) {
      assert(x->cols >= n);

      for (k = 0; k < n; k++) {
	if (!aunit) {
	  alpha = 1.0 / aa[k + k * lda];
	  for (i = 0; i < x->rows; i++)
	    xa[i + k * ldx] *= alpha;
	}
	for (i = 0; i < x->rows; i++)
	  for (j = k + 1; j < n; j++)
	    xa[i + j * ldx] -= xa[i + k * ldx] * aa[k + j * lda];
      }
    }
    else {
      assert(x->rows >= n);

      for (k = 0; k < n; k++) {
	if (!aunit) {
	  alpha = 1.0 / CONJ(aa[k + k * lda]);
	  for (j = 0; j < x->cols; j++)
	    xa[k + j * ldx] *= alpha;
	}
	for (i = k + 1; i < n; i++)
	  for (j = 0; j < x->cols; j++)
	    xa[i + j * ldx] -= CONJ(aa[k + i * lda]) * xa[k + j * ldx];
      }
    }
  }
  else {
    if (xtrans) {
      assert(x->cols >= n);

      for (k = n; k-- > 0;) {
	if (!aunit) {
	  alpha = 1.0 / CONJ(aa[k + k * lda]);
	  for (i = 0; i < x->rows; i++)
	    xa[i + k * ldx] *= alpha;
	}
	for (i = 0; i < x->rows; i++)
	  for (j = 0; j < k; j++)
	    xa[i + j * ldx] -= xa[i + k * ldx] * CONJ(aa[j + k * lda]);
      }
    }
    else {
      assert(x->rows >= n);

      for (k = n; k-- > 0;) {
	if (!aunit) {
	  alpha = 1.0 / aa[k + k * lda];
	  for (j = 0; j < x->cols; j++)
	    xa[k + j * ldx] *= alpha;
	}
	for (i = 0; i < k; i++)
	  for (j = 0; j < x->cols; j++)
	    xa[i + j * ldx] -= aa[i + k * lda] * xa[k + j * ldx];
      }
    }
  }
}
#endif

void
triangularsolve_amatrix(bool alower, bool aunit, bool atrans, pcamatrix a,
			bool xtrans, pamatrix x)
{
  if (alower)
    lowersolve_amatrix(aunit, atrans, a, xtrans, x);
  else
    uppersolve_amatrix(aunit, atrans, a, xtrans, x);
}

#ifdef USE_BLAS
static void
lowereval_amatrix_avector(bool aunit, bool atrans, pcamatrix a, pavector x)
{
  field    *aa = a->a;
  field    *xv = x->v;
  int      lda = a->ld;
  int      n = UINT_MIN(a->rows, a->cols);
  int      n1, i;

  assert(x->dim >= a->rows);
  assert(x->dim >= a->cols);

  if (n == 0)			/* Quick exit */
    return;

  if (atrans) {
    /* Left upper part, upper triangular */
    h2_trmv(h2_col_major, h2_lower, h2_ctrans, (aunit ? h2_udiag : h2_ndiag), n, aa, lda, xv, i_one);

    /* Right part */
    if (n < a->rows) {
      n1 = a->rows - n;
      h2_gemv(h2_col_major, h2_ctrans, n1, n, &f_one, aa + n, lda, xv + n, i_one, &f_one, xv, i_one);
    }

    /* Lower part */
    if (n < a->cols)
      for (i = n; i < a->cols; i++)
	xv[i] = 0.0;
  }
  else {
    /* Lower part */
    if (n < a->rows) {
      for (i = n; i < a->rows; i++)
	xv[i] = 0.0;

      n1 = a->rows - n;
      h2_gemv(h2_col_major, h2_ntrans, n1, n, &f_one, aa + n, lda, xv, i_one, &f_one, xv + n, i_one);
    }

    /* Top part, lower triangular */
    h2_trmv(h2_col_major, h2_lower, h2_ntrans, (aunit ? h2_udiag : h2_ndiag), n, aa, lda, xv, i_one);
  }
}

static void
uppereval_amatrix_avector(bool aunit, bool atrans, pcamatrix a, pavector x)
{
  field    *aa = a->a;
  field    *xv = x->v;
  int      lda = a->ld;
  int      n = UINT_MIN(a->rows, a->cols);
  int      n1, i;

  assert(x->dim >= a->rows);
  assert(x->dim >= a->cols);

  if (n == 0)			/* Quick exit */
    return;

  if (atrans) {
    /* Lower part */
    if (n < a->cols) {
      for (i = n; i < a->cols; i++)
	xv[i] = 0.0;

      n1 = a->cols - n;
      h2_gemv(h2_col_major, h2_ctrans, n, n1, &f_one, aa + n * lda, lda, xv, i_one, &f_one, xv + n, i_one);
    }

    /* Top part, lower triangular */
    h2_trmv(h2_col_major, h2_upper, h2_ctrans, (aunit ? h2_udiag : h2_ndiag), n, aa, lda, xv, i_one);
  }
  else {
    /* Left upper part, upper triangular */
    h2_trmv(h2_col_major, h2_upper, h2_ntrans, (aunit ? h2_udiag : h2_ndiag), n, aa, lda, xv, i_one);

    /* Right part */
    if (n < a->cols) {
      n1 = a->cols - n;
      h2_gemv(h2_col_major, h2_ntrans, n, n1, &f_one, aa + n * lda, lda, xv + n, i_one, &f_one, xv, i_one);
    }

    /* Lower part */
    if (n < a->rows)
      for (i = n; i < a->rows; i++)
	xv[i] = 0.0;
  }
}
#else
static void
lowereval_amatrix_avector(bool aunit, bool atrans, pcamatrix a, pavector x)
{
  pcfield   aa = a->a;
  int      lda = a->ld;
  pfield    xv = x->v;
  int      n = UINT_MIN(a->rows, a->cols);
  field     newval;
  int      i, j;

  assert(x->dim >= a->rows);
  assert(x->dim >= a->cols);

  if (atrans) {
    /* Left upper part, upper triangular */
    for (i = 0; i < n; i++) {
      newval = (aunit ? xv[i] : CONJ(aa[i + i * lda]) * xv[i]);
      for (j = i + 1; j < n; j++)
	newval += CONJ(aa[j + i * lda]) * xv[j];
      xv[i] = newval;
    }

    /* Right part */
    if (n < a->rows)
      for (i = 0; i < n; i++) {
	newval = xv[i];
	for (j = n; j < a->rows; j++)
	  newval += CONJ(aa[j + i * lda]) * xv[j];
	xv[i] = newval;
      }

    /* Lower part */
    if (n < a->cols)
      for (i = n; i < a->cols; i++)
	xv[i] = 0.0;
  }
  else {
    /* Lower part */
    if (n < a->rows)
      for (i = n; i < a->rows; i++) {
	newval = 0.0;
	for (j = 0; j < a->cols; j++)
	  newval += aa[i + j * lda] * xv[j];
	xv[i] = newval;
      }

    /* Top part, lower triangular */
    for (i = n; i-- > 0;) {
      newval = (aunit ? xv[i] : aa[i + i * lda] * xv[i]);
      for (j = 0; j < i; j++)
	newval += aa[i + j * lda] * xv[j];
      xv[i] = newval;
    }
  }
}

static void
uppereval_amatrix_avector(bool aunit, bool atrans, pcamatrix a, pavector x)
{
  pcfield   aa = a->a;
  int      lda = a->ld;
  pfield    xv = x->v;
  int      n = UINT_MIN(a->rows, a->cols);
  field     newval;
  int      i, j;

  assert(x->dim >= a->rows);
  assert(x->dim >= a->cols);

  if (atrans) {
    /* Lower part */
    if (n < a->cols)
      for (i = n; i < a->cols; i++) {
	newval = 0.0;
	for (j = 0; j < n; j++)
	  newval += CONJ(aa[j + i * lda]) * xv[j];
	xv[i] = newval;
      }

    /* Top part, lower triangular */
    for (i = n; i-- > 0;) {
      newval = (aunit ? xv[i] : CONJ(aa[i + i * lda]) * xv[i]);
      for (j = 0; j < i; j++)
	newval += CONJ(aa[j + i * lda]) * xv[j];
      xv[i] = newval;
    }
  }
  else {
    /* Left upper part, upper triangular */
    for (i = 0; i < n; i++) {
      newval = (aunit ? xv[i] : xv[i] * aa[i + i * lda]);
      for (j = i + 1; j < n; j++)
	newval += aa[i + j * lda] * xv[j];
      xv[i] = newval;
    }

    /* Right part */
    if (n < a->cols)
      for (i = 0; i < n; i++) {
	newval = xv[i];
	for (j = n; j < a->cols; j++)
	  newval += aa[i + j * lda] * xv[j];
	xv[i] = newval;
      }

    /* Lower part */
    if (n < a->rows)
      for (i = n; i < a->rows; i++)
	xv[i] = 0.0;
  }
}
#endif

void
triangulareval_amatrix_avector(bool alower, bool aunit, bool atrans,
			       pcamatrix a, pavector x)
{
  if (alower)
    lowereval_amatrix_avector(aunit, atrans, a, x);
  else
    uppereval_amatrix_avector(aunit, atrans, a, x);
}

#ifdef USE_BLAS
static void
lowereval_amatrix(bool aunit, bool atrans, pcamatrix a, bool xtrans,
		  pamatrix x)
{
  field    *aa = a->a;
  int      lda = a->ld;
  field    *xa = x->a;
  int      ldx = x->ld;
  int      n = UINT_MIN(a->rows, a->cols);
  int      n1, i, j;

  if (n == 0)			/* Quick exit */
    return;

  if (xtrans) {
    assert(x->cols >= a->rows);
    assert(x->cols >= a->cols);

    if (atrans) {
      /* Left upper part, upper triangular */
      h2_trmm(h2_col_major, h2_right, h2_lower, h2_ntrans,
	      (aunit ? h2_udiag : h2_ndiag), x->rows, n, &f_one, aa,
	      lda, xa, ldx);

      /* Right part */
      if (n < a->rows) {
	n1 = a->rows - n;
	h2_gemm(h2_col_major, h2_ntrans, h2_ntrans, x->rows, n, n1, &f_one, xa + n * ldx, ldx, aa + n, lda, &f_one, xa, ldx);
      }

      /* Lower part */
      if (n < a->cols)
	for (j = 0; j < x->rows; j++)
	  for (i = n; i < a->cols; i++)
	    xa[j + i * ldx] = 0.0;
    }
    else {
      /* Lower part */
      if (n < a->rows) {
	for (j = 0; j < x->rows; j++)
	  for (i = n; i < a->rows; i++)
	    xa[j + i * ldx] = 0.0;

	n1 = a->rows - n;
	h2_gemm(h2_col_major, h2_ntrans, h2_ctrans, x->rows, n1, n, &f_one, xa, ldx, aa + n, lda, &f_one, xa + n * ldx, ldx);
      }

      /* Top part, lower triangular */
      h2_trmm(h2_col_major, h2_right, h2_lower, h2_ctrans, (aunit ? h2_udiag : h2_ndiag),
	      x->rows, n, &f_one, aa, lda, xa, ldx);
    }
  }
  else {
    assert(x->rows >= a->rows);
    assert(x->rows >= a->cols);

    if (atrans) {
      /* Left upper part, upper triangular */
      h2_trmm(h2_col_major, h2_left, h2_lower, h2_ctrans, (aunit ? h2_udiag : h2_ndiag), n, x->cols, &f_one, aa, lda, xa, ldx);

      /* Right part */
      if (n < a->rows) {
	n1 = a->rows - n;
	h2_gemm(h2_col_major, h2_ctrans, h2_ntrans, n, x->cols, n1, &f_one, aa + n, lda, xa + n, ldx, &f_one, xa, ldx);
      }

      /* Lower part */
      if (n < a->cols)
	for (j = 0; j < x->cols; j++)
	  for (i = n; i < a->cols; i++)
	    xa[i + j * ldx] = 0.0;
    }
    else {
      /* Lower part */
      if (n < a->rows) {
	for (j = 0; j < x->cols; j++)
	  for (i = n; i < a->rows; i++)
	    xa[i + j * ldx] = 0.0;

	n1 = a->rows - n;
	h2_gemm(h2_col_major, h2_ntrans, h2_ntrans, n1, x->cols, n, &f_one, aa + n, lda, xa, ldx, &f_one, xa + n, ldx);
      }

      /* Top part, lower triangular */
      h2_trmm(h2_col_major, h2_left, h2_lower, h2_ntrans,
	      (aunit ? h2_udiag : h2_ndiag), n, x->cols, &f_one, aa,
	      lda, xa, ldx);
    }
  }
}

static void
uppereval_amatrix(bool aunit, bool atrans, pcamatrix a, bool xtrans,
		  pamatrix x)
{
  field    *aa = a->a;
  int      lda = a->ld;
  field    *xa = x->a;
  int      ldx = x->ld;
  int      n = UINT_MIN(a->rows, a->cols);
  int      n1, i, j;

  if (n == 0)			/* Quick exit */
    return;

  if (xtrans) {
    assert(x->cols >= a->rows);
    assert(x->cols >= a->cols);

    if (atrans) {
      /* Lower part */
      if (n < a->cols) {
	for (j = 0; j < x->rows; j++)
	  for (i = n; i < a->cols; i++)
	    xa[j + i * ldx] = 0.0;

	n1 = a->cols - n;
	h2_gemm(h2_col_major, h2_ntrans, h2_ntrans, x->rows, n1, n, &f_one, xa, ldx, aa + n * lda, lda, &f_one, xa + n * ldx, ldx);
      }

      /* Top part, lower triangular */
      h2_trmm(h2_col_major, h2_right, h2_upper, h2_ntrans,
	      (aunit ? h2_udiag : h2_ndiag), x->rows, n, &f_one, aa,
	      lda, xa, ldx);
    }
    else {
      /* Left upper part, upper triangular */
      h2_trmm(h2_col_major, h2_right, h2_upper, h2_ctrans, (aunit ? h2_udiag : h2_ndiag), x->rows, n, &f_one, aa, lda, xa, ldx);

      /* Right part */
      if (n < a->cols) {
	n1 = a->cols - n;
	h2_gemm(h2_col_major, h2_ntrans, h2_ctrans, x->rows, n, n1, &f_one, xa + n * ldx, ldx, aa + n * lda, lda, &f_one, xa, ldx);
      }

      /* Lower part */
      if (n < a->rows)
	for (j = 0; j < x->rows; j++)
	  for (i = n; i < a->rows; i++)
	    xa[j + i * ldx] = 0.0;
    }
  }
  else {
    assert(x->rows >= a->rows);
    assert(x->rows >= a->cols);

    if (atrans) {
      /* Lower part */
      if (n < a->cols) {
	for (j = 0; j < x->cols; j++)
	  for (i = n; i < a->cols; i++)
	    xa[i + j * ldx] = 0.0;

	n1 = a->cols - n;
	h2_gemm(h2_col_major, h2_ctrans, h2_ntrans, n1, x->cols, n, &f_one, aa + n * lda,	lda, xa, ldx, &f_one, xa + n, ldx);
      }

      /* Top part, lower triangular */
      h2_trmm(h2_col_major, h2_left, h2_upper, h2_ctrans, (aunit ? h2_udiag : h2_ndiag), n, x->cols, &f_one, aa, lda, xa, ldx);
    }
    else {
      /* Left upper part, upper triangular */
      h2_trmm(h2_col_major, h2_left, h2_upper, h2_ntrans, (aunit ? h2_udiag : h2_ndiag), n, x->cols, &f_one, aa, lda, xa, ldx);

      /* Right part */
      if (n < a->cols) {
	n1 = a->cols - n;
	h2_gemm(h2_col_major, h2_ntrans, h2_ntrans, n, x->cols, n1, &f_one, aa + n * lda, lda, xa + n, ldx, &f_one, xa, ldx);
      }

      /* Lower part */
      if (n < a->rows)
	for (j = 0; j < x->cols; j++)
	  for (i = n; i < a->rows; i++)
	    xa[i + j * ldx] = 0.0;
    }
  }
}
#else
static void
lowereval_amatrix(bool aunit, bool atrans, pcamatrix a, bool xtrans,
		  pamatrix x)
{
  pcfield   aa = a->a;
  int      lda = a->ld;
  pfield    xa = x->a;
  int      ldx = x->ld;
  int      n = UINT_MIN(a->rows, a->cols);
  field     newval;
  int      i, j, k;

  if (xtrans) {
    assert(x->cols >= a->rows);
    assert(x->cols >= a->cols);

    if (atrans) {
      for (k = 0; k < x->rows; k++) {
	/* Left upper part, upper triangular */
	for (i = 0; i < n; i++) {
	  newval =
	    (aunit ? xa[k + i * ldx] : aa[i + i * lda] * xa[k + i * ldx]);
	  for (j = i + 1; j < n; j++)
	    newval += aa[j + i * lda] * xa[k + j * ldx];
	  xa[k + i * ldx] = newval;
	}

	/* Right part */
	if (n < a->rows)
	  for (i = 0; i < n; i++) {
	    newval = xa[k + i * ldx];
	    for (j = n; j < a->rows; j++)
	      newval += aa[j + i * lda] * xa[k + j * ldx];
	    xa[k + i * ldx] = newval;
	  }

	/* Lower part */
	if (n < a->cols)
	  for (i = n; i < a->cols; i++)
	    xa[k + i * ldx] = 0.0;
      }
    }
    else {
      for (k = 0; k < x->rows; k++) {
	/* Lower part */
	if (n < a->rows)
	  for (i = n; i < a->rows; i++) {
	    newval = 0.0;
	    for (j = 0; j < a->cols; j++)
	      newval += CONJ(aa[i + j * lda]) * xa[k + j * ldx];
	    xa[k + i * ldx] = newval;
	  }

	/* Top part, lower triangular */
	for (i = n; i-- > 0;) {
	  newval =
	    (aunit ? xa[k + i * ldx] : CONJ(aa[i + i * lda]) *
	     xa[k + i * ldx]);
	  for (j = 0; j < i; j++)
	    newval += CONJ(aa[i + j * lda]) * xa[k + j * ldx];
	  xa[k + i * ldx] = newval;
	}
      }
    }
  }
  else {
    assert(x->rows >= a->rows);
    assert(x->rows >= a->cols);

    if (atrans) {
      for (k = 0; k < x->cols; k++) {
	/* Left upper part, upper triangular */
	for (i = 0; i < n; i++) {
	  newval =
	    (aunit ? xa[i + k * ldx] : CONJ(aa[i + i * lda]) *
	     xa[i + k * ldx]);
	  for (j = i + 1; j < n; j++)
	    newval += CONJ(aa[j + i * lda]) * xa[j + k * ldx];
	  xa[i + k * ldx] = newval;
	}

	/* Right part */
	if (n < a->rows)
	  for (i = 0; i < n; i++) {
	    newval = xa[i + k * ldx];
	    for (j = n; j < a->rows; j++)
	      newval += CONJ(aa[j + i * lda]) * xa[j + k * ldx];
	    xa[i + k * ldx] = newval;
	  }

	/* Lower part */
	if (n < a->cols)
	  for (i = n; i < a->cols; i++)
	    xa[i + k * ldx] = 0.0;
      }
    }
    else {
      for (k = 0; k < x->cols; k++) {
	/* Lower part */
	if (n < a->rows)
	  for (i = n; i < a->rows; i++) {
	    newval = 0.0;
	    for (j = 0; j < a->cols; j++)
	      newval += aa[i + j * lda] * xa[j + k * ldx];
	    xa[i + k * ldx] = newval;
	  }

	/* Top part, lower triangular */
	for (i = n; i-- > 0;) {
	  newval =
	    (aunit ? xa[i + k * ldx] : aa[i + i * lda] * xa[i + k * ldx]);
	  for (j = 0; j < i; j++)
	    newval += aa[i + j * lda] * xa[j + k * ldx];
	  xa[i + k * ldx] = newval;
	}
      }
    }
  }
}

static void
uppereval_amatrix(bool aunit, bool atrans, pcamatrix a, bool xtrans,
		  pamatrix x)
{
  pcfield   aa = a->a;
  int      lda = a->ld;
  pfield    xa = x->a;
  int      ldx = x->ld;
  int      n = UINT_MIN(a->rows, a->cols);
  field     newval;
  int      i, j, k;

  if (xtrans) {
    assert(x->cols >= a->rows);
    assert(x->cols >= a->cols);

    if (atrans) {
      for (k = 0; k < x->rows; k++) {
	/* Lower part */
	if (n < a->cols)
	  for (i = n; i < a->cols; i++) {
	    newval = 0.0;
	    for (j = 0; j < a->rows; j++)
	      newval += aa[j + i * lda] * xa[k + j * ldx];
	    xa[k + i * ldx] = newval;
	  }

	/* Top part, lower triangular */
	for (i = n; i-- > 0;) {
	  newval =
	    (aunit ? xa[k + i * ldx] : aa[i + i * lda] * xa[k + i * ldx]);
	  for (j = 0; j < i; j++)
	    newval += aa[j + i * lda] * xa[k + j * ldx];
	  xa[k + i * ldx] = newval;
	}
      }
    }
    else {
      for (k = 0; k < x->rows; k++) {
	/* Left upper part, upper triangular */
	for (i = 0; i < n; i++) {
	  newval =
	    (aunit ? xa[k + i * ldx] : CONJ(aa[i + i * lda]) *
	     xa[k + i * ldx]);
	  for (j = i + 1; j < n; j++)
	    newval += CONJ(aa[i + j * lda]) * xa[k + j * ldx];
	  xa[k + i * ldx] = newval;
	}

	/* Right part */
	if (n < a->cols)
	  for (i = 0; i < n; i++) {
	    newval = xa[k + i * ldx];
	    for (j = n; j < a->cols; j++)
	      newval += CONJ(aa[i + j * lda]) * xa[k + j * ldx];
	    xa[k + i * ldx] = newval;
	  }

	/* Lower part */
	if (n < a->rows)
	  for (i = n; i < a->rows; i++)
	    xa[k + i * ldx] = 0.0;
      }
    }
  }
  else {
    assert(x->rows >= a->rows);
    assert(x->rows >= a->cols);

    if (atrans) {
      for (k = 0; k < x->cols; k++) {
	/* Lower part */
	if (n < a->cols)
	  for (i = n; i < a->cols; i++) {
	    newval = 0.0;
	    for (j = 0; j < a->rows; j++)
	      newval += CONJ(aa[j + i * lda]) * xa[j + k * ldx];
	    xa[i + k * ldx] = newval;
	  }

	/* Top part, lower triangular */
	for (i = n; i-- > 0;) {
	  newval =
	    (aunit ? xa[i + k * ldx] : CONJ(aa[i + i * lda]) *
	     xa[i + k * ldx]);
	  for (j = 0; j < i; j++)
	    newval += CONJ(aa[j + i * lda]) * xa[j + k * ldx];
	  xa[i + k * ldx] = newval;
	}
      }
    }
    else {
      for (k = 0; k < x->cols; k++) {
	/* Left upper part, upper triangular */
	for (i = 0; i < n; i++) {
	  newval =
	    (aunit ? xa[i + k * ldx] : aa[i + i * lda] * xa[i + k * ldx]);
	  for (j = i + 1; j < n; j++)
	    newval += aa[i + j * lda] * xa[j + k * ldx];
	  xa[i + k * ldx] = newval;
	}

	/* Right part */
	if (n < a->cols)
	  for (i = 0; i < n; i++) {
	    newval = xa[i + k * ldx];
	    for (j = n; j < a->cols; j++)
	      newval += aa[i + j * lda] * xa[j + k * ldx];
	    xa[i + k * ldx] = newval;
	  }

	/* Lower part */
	if (n < a->rows)
	  for (i = n; i < a->rows; i++)
	    xa[i + k * ldx] = 0.0;
      }
    }
  }
}
#endif

void
triangulareval_amatrix(bool alower, bool aunit, bool atrans, pcamatrix a,
		       bool xtrans, pamatrix x)
{
  if (alower)
    lowereval_amatrix(aunit, atrans, a, xtrans, x);
  else
    uppereval_amatrix(aunit, atrans, a, xtrans, x);
}

void
triangularaddmul_amatrix(field alpha, bool alower, bool atrans,
			 pcamatrix a, bool blower, bool btrans, pcamatrix b,
			 pamatrix c)
{
  pcfield   aa = a->a;
  int      lda = a->ld;
  pcfield   ba = b->a;
  int      ldb = b->ld;
  pfield    ca = c->a;
  int      ldc = c->ld;
  int      aoff, adim, ainc, boff, bdim, binc;
  int      j;
  //#ifndef USE_BLAS
  int      i, k;
  //#endif

  if (atrans) {
    assert(c->rows == a->cols);

    ainc = lda;
    lda = 1;

    if (btrans) {		/* C += alpha A^* B^* */
      assert(a->rows == b->cols);
      assert(c->cols == b->rows);

      binc = 1;

      for (j = 0; j < a->rows; j++) {
	if (alower) {		/* A^* upper triangular */
	  aoff = 0;
	  adim = UINT_MIN(j + 1, a->cols);
	}
	else {			/* A^* lower triangular */
	  aoff = j;
	  adim = a->cols - UINT_MIN(j, a->cols);
	}

	if (blower) {		/* B^* upper triangular */
	  boff = j;
	  bdim = b->rows - UINT_MIN(j, b->rows);
	}
	else {			/* B^* lower triangular */
	  boff = 0;
	  bdim = UINT_MIN(j + 1, b->rows);
	}

	//TODO Need a 'ger' equivalent for A \gets A + alpha * CONJ(x) * y ** H
	//#ifdef USE_BLAS
	//        h2_geru(&adim, &bdim, &alpha, aa + aoff * ainc + j * lda, &ainc,
	//            ba + boff * binc + j * ldb, &binc, ca + aoff + boff * ldc, &ldc);
	//#else

	for (k = 0; k < bdim; k++) {
	  for (i = 0; i < adim; i++) {
	    ca[(aoff + i) + (boff + k) * ldc] += alpha
	      * CONJ(aa[(aoff + i) * ainc + j * lda]) *
	      CONJ(ba[(boff + k) * binc + j * ldb]);
	  }
	}
	//#endif
      }
    }
    else {			/* C += alpha A^* B */
      assert(a->rows == b->rows);
      assert(c->cols == b->cols);

      binc = ldb;
      ldb = 1;

      for (j = 0; j < a->rows; j++) {
	if (alower) {		/* A^* upper triangular */
	  aoff = 0;
	  adim = UINT_MIN(j + 1, a->cols);
	}
	else {			/* A^* lower triangular */
	  aoff = j;
	  adim = a->cols - UINT_MIN(j, a->cols);
	}

	if (blower) {		/* B lower triangular */
	  boff = 0;
	  bdim = UINT_MIN(j + 1, b->cols);
	}
	else {			/* B upper triangular */
	  boff = j;
	  bdim = b->cols - UINT_MIN(j, b->cols);
	}

	//TODO Need a 'ger' equivalent for A \gets A + alpha * CONJ(x) * y ** T
	//#ifdef USE_BLAS
	//        h2_gerc(&adim, &bdim, &alpha, aa + aoff * ainc + j * lda, &ainc,
	//            ba + boff * binc + j * ldb, &binc, ca + aoff + boff * ldc, &ldc);
	//#else
	for (k = 0; k < bdim; k++) {
	  for (i = 0; i < adim; i++) {
	    ca[(aoff + i) + (boff + k) * ldc] += alpha
	      * CONJ(aa[(aoff + i) * ainc + j * lda])
	      * ba[(boff + k) * binc + j * ldb];
	  }
	}
	//#endif
      }
    }
  }
  else {
    assert(c->rows == a->rows);

    ainc = 1;

    if (btrans) {		/* C += alpha A B^* */
      assert(a->cols == b->cols);
      assert(c->cols == b->rows);

      binc = 1;

      for (j = 0; j < a->cols; j++) {
	if (alower) {		/* A lower triangular */
	  aoff = j;
	  adim = a->rows - UINT_MIN(j, a->rows);
	}
	else {			/* A upper triangular */
	  aoff = 0;
	  adim = UINT_MIN(j + 1, a->rows);
	}

	if (blower) {		/* B^* upper triangular */
	  boff = j;
	  bdim = b->rows - UINT_MIN(j, b->rows);
	}
	else {			/* B^* lower triangular */
	  boff = 0;
	  bdim = UINT_MIN(j + 1, b->rows);
	}

#ifdef USE_BLAS
	h2_gerc(h2_col_major, adim, bdim, &alpha, aa + aoff * ainc + j * lda, ainc, ba + boff * binc + j * ldb, binc, ca + aoff + boff * ldc, ldc);
#else
	for (k = 0; k < bdim; k++) {
	  for (i = 0; i < adim; i++) {
	    ca[(aoff + i) + (boff + k) * ldc] += alpha
	      * aa[(aoff + i) * ainc + j * lda]
	      * CONJ(ba[(boff + k) * binc + j * ldb]);
	  }
	}
#endif
      }
    }
    else {			/* C += alpha A B */
      assert(a->cols == b->rows);
      assert(c->cols == b->cols);

      binc = ldb;
      ldb = 1;

      for (j = 0; j < a->cols; j++) {
	if (alower) {		/* A lower triangular */
	  aoff = j;
	  adim = a->rows - UINT_MIN(j, a->rows);
	}
	else {			/* A upper triangular */
	  aoff = 0;
	  adim = UINT_MIN(j + 1, a->rows);
	}

	if (blower) {		/* B lower triangular */
	  boff = 0;
	  bdim = UINT_MIN(j + 1, b->cols);
	}
	else {			/* B upper triangular */
	  boff = j;
	  bdim = b->cols - UINT_MIN(j, b->cols);
	}

#ifdef USE_BLAS
	h2_geru(h2_col_major, adim, bdim, &alpha, aa + aoff * ainc + j * lda, ainc, ba + boff * binc + j * ldb, binc, ca + aoff + boff * ldc, ldc);
#else

	for (k = 0; k < bdim; k++) {
	  for (i = 0; i < adim; i++) {
	    ca[(aoff + i) + (boff + k) * ldc] += alpha
	      * aa[(aoff + i) * ainc + j * lda]
	      * ba[(boff + k) * binc + j * ldb];
	  }
	}
#endif
      }
    }
  }
}

void
copy_lower_amatrix(pcamatrix a, bool aunit, pamatrix b)
{
  pfield    aa = a->a;
  int      lda = a->ld;
  pfield    ba = b->a;
  int      ldb = b->ld;
  int      rows;
  int      cols;
  int      i, j;

  rows = UINT_MIN(a->rows, b->rows);
  cols = UINT_MIN(a->cols, b->cols);

  for (j = 0; j < cols; j++) {
    if (aunit) {
      for (i = 0; i < rows && i < j; i++) {
	ba[i + j * ldb] = 0.0;
      }

      if (i < rows && i == j) {
	ba[i + j * ldb] = 1.0;
	i++;
      }
    }
    else {
      for (i = 0; i < rows && i < j; i++) {
	ba[i + j * ldb] = 0.0;
      }
    }

    for (; i < rows; i++) {
      ba[i + j * ldb] = aa[i + j * lda];
    }
  }
  for (; j < b->cols; j++) {
    for (i = 0; i < rows; i++) {
      ba[i + j * ldb] = 0.0;
    }
  }
}

void
copy_upper_amatrix(pcamatrix a, bool aunit, pamatrix b)
{
  pfield    aa = a->a;
  int      lda = a->ld;
  pfield    ba = b->a;
  int      ldb = b->ld;
  int      rows;
  int      cols;
  int      i, j;

  rows = UINT_MIN(a->rows, b->rows);
  cols = UINT_MIN(a->cols, b->cols);

  for (j = 0; j < cols; j++) {
    if (aunit) {
      for (i = 0; i < rows && i < j; i++) {
	ba[i + j * ldb] = aa[i + j * lda];
      }

      if (i < rows && i == j) {
	ba[i + j * ldb] = 1.0;
	i++;
      }
    }
    else {
      for (i = 0; i < rows && i <= j; i++) {
	ba[i + j * ldb] = aa[i + j * lda];
      }
    }

    for (; i < b->rows; i++) {
      ba[i + j * ldb] = 0.0;
    }
  }
}

/* ------------------------------------------------------------
 LR decomposition
 ------------------------------------------------------------ */

#ifdef USE_BLAS
int
lrdecomp_amatrix(pamatrix a)
{
  field    *aa = a->a;
  int      lda = a->ld;
  int      n = a->rows;
  field     alpha;
  int      i, n1;

  assert(n == a->cols);

  for (i = 0; i < n - 1; i++) {
	  if (aa[i + i * lda] == 0.0)
		  return i + 1;

	  alpha = 1.0 / aa[i + i * lda];

	  n1 = n - i - 1;

  	  h2_scal(n1, &alpha, aa + (i + 1) + i * lda, i_one);

  	  h2_geru(h2_col_major, n1, n1, &f_minusone, aa + (i + 1) + i * lda, i_one, aa + i + (i + 1) * lda, lda, aa + (i + 1) + (i + 1) * lda, lda);

  }

  if (aa[i + i * lda] == 0.0)
    return i + 1;

  return 0;
}
#else
int
lrdecomp_amatrix(pamatrix a)
{
  pfield    aa = a->a;
  int      lda = a->ld;
  int      n = a->rows;
  field     alpha;
  int      i, j, k;

  assert(n == a->cols);

  for (i = 0; i < n - 1; i++) {
    if (aa[i + i * lda] == 0.0)
      return i + 1;

    alpha = 1.0 / aa[i + i * lda];

    for (j = i + 1; j < n; j++) {
      aa[j + i * lda] *= alpha;
      for (k = i + 1; k < n; k++)
	aa[j + k * lda] -= aa[j + i * lda] * aa[i + k * lda];
    }
  }

  if (aa[i + i * lda] == 0.0)
    return i + 1;

  return 0;
}
#endif

void
lrsolve_amatrix_avector(pcamatrix a, pavector x)
{
  triangularsolve_amatrix_avector(true, true, false, a, x);
  triangularsolve_amatrix_avector(false, false, false, a, x);
}

void
lrsolve_amatrix(pcamatrix a, pamatrix x)
{
  triangularsolve_amatrix(true, true, false, a, false, x);
  triangularsolve_amatrix(false, false, false, a, false, x);
}

/* ------------------------------------------------------------
 Cholesky decomposition
 ------------------------------------------------------------ */

#ifdef USE_BLAS
int
choldecomp_amatrix(pamatrix a)
{
  field    *aa = a->a;
  int      lda = a->ld;
  int      n = a->rows;

  int       info;

  assert(n == a->cols);

  info = h2_potrf(h2_col_major, h2_lower, n, (MKL_Complex16 *)aa, lda);

  return info;
}
#else
int
choldecomp_amatrix(pamatrix a)
{
  pfield    aa = a->a;
  int      lda = a->ld;
  int      n = a->rows;
  double      diag;
  field		alpha;
  int      i, j, k;

  assert(n == a->cols);

  for (i = 0; i < n - 1; i++) {
    diag = REAL(aa[i + i * lda]);

    if (ABS(aa[i + i * lda] - diag) > 1e-12 || diag <= 0.0)
      return i + 1;

    aa[i + i * lda] = REAL_SQRT(diag);
    alpha = 1.0 / aa[i + i * lda];
    for (j = i + 1; j < n; j++)
      aa[j + i * lda] *= alpha;

    for (j = i + 1; j < n; j++)
      for (k = i + 1; k <= j; k++)
	aa[j + k * lda] -= aa[j + i * lda] * CONJ(aa[k + i * lda]);
  }

  diag = REAL(aa[i + i * lda]);
  if (ABS(aa[i + i * lda] - diag) > 1e-12 || diag <= 0.0)
    return i + 1;

  aa[i + i * lda] = REAL_SQRT(diag);

  return 0;
}
#endif

void
cholsolve_amatrix_avector(pcamatrix a, pavector x)
{
  triangularsolve_amatrix_avector(true, false, false, a, x);
  triangularsolve_amatrix_avector(true, false, true, a, x);
}

void
cholsolve_amatrix(pcamatrix a, pamatrix x)
{
  assert(a->cols == x->rows);

  triangularsolve_amatrix(true, false, false, a, false, x);
  triangularsolve_amatrix(true, false, true, a, false, x);
}

/* ------------------------------------------------------------
 LDL^T decomposition
 ------------------------------------------------------------ */

#ifdef USE_BLAS
int
ldltdecomp_amatrix(pamatrix a)
{
  field    *aa = a->a;
  int      lda = a->ld;
  int      n = a->rows;
  double      diag, alpha;
  int      i, n1;

  assert(n == a->cols);

  for (i = 0; i < n - 1; i++) {
    diag = REAL(aa[i + i * lda]);

    if (ABS(aa[i + i * lda] - diag) > 1e-12 || diag == 0.0)
      return i + 1;

    alpha = 1.0 / diag;
    n1 = n - i - 1;
    h2_rscal(n1, alpha, aa + (i + 1) + i * lda, i_one);

    alpha = -diag;
    h2_syr(h2_col_major, h2_lower, n1, alpha, aa + (i + 1) + i * lda, i_one, aa + (i + 1) + (i + 1) * lda, lda);
  }

  diag = REAL(aa[i + i * lda]);
  if (ABS(aa[i + i * lda] - diag) > 1e-12 || diag == 0.0)
    return i + 1;

  return 0;
}
#else
int
ldltdecomp_amatrix(pamatrix a)
{
  pfield    aa = a->a;
  int      lda = a->ld;
  int      n = a->rows;
  double      diag, alpha;
  int      i, j, k;

  assert(n == a->cols);

  for (i = 0; i < n - 1; i++) {
    diag = REAL(aa[i + i * lda]);

    if (ABS(aa[i + i * lda] - diag) > 1e-12 || diag == 0.0)
      return i + 1;

    alpha = 1.0 / diag;
    for (j = i + 1; j < n; j++)
      aa[j + i * lda] *= alpha;

    for (j = i + 1; j < n; j++)
      for (k = i + 1; k <= j; k++)
	aa[j + k * lda] -= diag * aa[j + i * lda] * CONJ(aa[k + i * lda]);
  }

  diag = REAL(aa[i + i * lda]);
  if (ABS(aa[i + i * lda] - diag) > 1e-12 || diag == 0.0)
    return i + 1;

  return 0;
}
#endif

void
ldltsolve_amatrix_avector(pcamatrix a, pavector x)
{
  triangularsolve_amatrix_avector(true, true, false, a, x);
  diagsolve_amatrix_avector(false, a, x);
  triangularsolve_amatrix_avector(true, true, true, a, x);
}

void
ldltsolve_amatrix(pcamatrix a, pamatrix x)
{
  assert(a->cols == x->rows);

  triangularsolve_amatrix(true, true, false, a, false, x);
  diagsolve_amatrix(false, a, false, x);
  triangularsolve_amatrix(true, true, true, a, false, x);
}

/* ------------------------------------------------------------
 Orthogonal decompositions
 ------------------------------------------------------------ */

#ifdef USE_BLAS
void
qrdecomp_amatrix(pamatrix a, pavector tau)
{
  int      rows = a->rows;
  int      cols = a->cols;
  int      refl = min<int>(rows, cols);
  //field    *work;
  //int       lwork;
  int		info;

  assert(a->ld >= rows);
  /* Quick exit if no reflections used */
  if (refl == 0)
    return;

  //lwork = 4 * cols;
  //work = allocfield(lwork);

  if (tau->dim < refl)
    resize_avector(tau, refl);

  //h2_geqrf(h2_col_major, rows, cols, a->a, a->ld, tau->v, work, &lwork, &info);
  info = h2_geqrf(h2_col_major, rows, cols, (MKL_Complex16 *)(a->a), a->ld, (MKL_Complex16 *)tau->v);
  assert(info == 0);

  //freemem(work);
}
#else
void
qrdecomp_amatrix(pamatrix a, pavector tau)
{
  pfield    aa = a->a;
  int      lda = a->ld;
  pfield    tauv;
  int      rows = a->rows;
  int      cols = a->cols;
  int      refl = UINT_MIN(rows, cols);
  field     alpha, beta, gamma, diag;
  double      norm2, norm;
  int      i, j, k;

  /* Provide enough storage for scaling factors */
  if (tau->dim < refl)
    resize_avector(tau, refl);
  tauv = tau->v;

  for (k = 0; k < refl; k++) {
    /* Compute norm of k-th column */
    norm2 = 0.0;
    for (i = k; i < rows; i++)
      norm2 += ABSSQR(aa[i + k * lda]);
    norm = REAL_SQRT(norm2);

    if (norm2 == 0.0)
      tauv[k] = 0.0;
    else {
      /* Determine reflection vector v */
      diag = aa[k + k * lda];
      alpha = -SIGN1(diag) * norm;

      /* Compute norm of v */
      beta = 1.0 / (norm2 - CONJ(alpha) * diag);

      /* Rescale to ensure v_1 = 1 */
      beta *= ABSSQR(diag - alpha);
      gamma = 1.0 / (diag - alpha);
      for (i = k + 1; i < rows; i++)
	aa[i + k * lda] *= gamma;
      tauv[k] = beta;

      /* Compute k-th column */
      aa[k + k * lda] = alpha;

      /* Update remaining columns */
      for (j = k + 1; j < cols; j++) {
	gamma = aa[k + j * lda];
	for (i = k + 1; i < rows; i++)
	  gamma += CONJ(aa[i + k * lda]) * aa[i + j * lda];

	gamma *= beta;

	aa[k + j * lda] -= gamma;
	for (i = k + 1; i < rows; i++)
	  aa[i + j * lda] -= gamma * aa[i + k * lda];
      }
    }
  }
}
#endif

#ifdef USE_BLAS
int
qrdecomp_pivot_amatrix(pamatrix a, pavector tau, int * colpiv)
{
  pfield    aa = a->a;
  int      lda = a->ld;
  pfield    tauv;
  int      rows = a->rows;
  int      cols = a->cols;
  int      refl = UINT_MIN(rows, cols);
  field     alpha, beta, gamma, diag;
  double      norm2, norm, maxnorm2;
  int      j, k, jmax, n1;

  /* Provide enough storage for scaling factors */
  if (tau->dim < refl)
    resize_avector(tau, refl);
  tauv = tau->v;

  /* Initialize column pivots */
  if (colpiv) {
    for (j = 0; j < cols; j++)
      colpiv[j] = j;
  }

  for (k = 0; k < refl; k++) {
    /* Compute norm of k-th column */
    n1 = rows - k;
    norm2 = REAL_SQR(h2_nrm2(n1, aa + k + k * lda, i_one));
    maxnorm2 = norm2;
    jmax = k;

    /* Find maximal norm */
    for (j = k + 1; j < cols; j++) {
      norm2 = REAL_SQR(h2_nrm2(n1, aa + k + j * lda, i_one));
      if (norm2 > maxnorm2) {
	maxnorm2 = norm2;
	jmax = j;
      }
    }

    /* Swap columns */
    if (jmax != k) {
      h2_swap(rows, aa + k * lda, i_one, aa + jmax * lda, i_one);

      if (colpiv) {
	j = colpiv[k];
	colpiv[k] = jmax;
	colpiv[jmax] = j;
      }
    }

    /* Remember pivot */
    if (colpiv)
      colpiv[k] = jmax;

    /* Prepare norm */
    norm2 = maxnorm2;
    norm = REAL_SQRT(norm2);

    if (norm2 == 0.0)
      tauv[k] = 0.0;
    else {
      /* Determine reflection vector v */
      diag = aa[k + k * lda];
      alpha = -SIGN1(diag) * norm;

      /* Compute norm of v */
      beta = 1.0 / (norm2 - CONJ(alpha) * diag);

      /* Rescale to ensure v_1 = 1 */
      beta *= ABSSQR(diag - alpha);
      gamma = 1.0 / (diag - alpha);
      n1 = rows - k - 1;
      h2_scal(n1, &gamma, aa + k + 1 + k * lda, i_one);
      tauv[k] = beta;

      /* Compute k-th column */
      aa[k + k * lda] = alpha;

      /* Update remaining columns */
      for (j = k + 1; j < cols; j++) {
	n1 = rows - k - 1;
	field tmp;
	h2_dot(n1, aa + k + 1 + k * lda, i_one, aa + k + 1 + j * lda, i_one, &tmp);
	gamma = aa[k + j * lda] + tmp;

	gamma *= -beta;

	aa[k + j * lda] += gamma;
	h2_axpy(n1, &gamma, aa + k + 1 + k * lda, i_one,
		aa + k + 1 + j * lda, i_one);
      }
    }
  }

  return refl;
}

int
qrdecomp_rank_amatrix(pamatrix a, pavector tau, pctruncmode tm, double eps,
		      int * colpiv)
{
  pfield    aa = a->a;
  int      lda = a->ld;
  pfield    tauv;
  int      rows = a->rows;
  int      cols = a->cols;
  int      refl = UINT_MIN(rows, cols);
  field     alpha, beta, gamma, diag;
  double      norm2, norm, maxnorm2, frobnorm2, firstnorm2;
  int      j, k, jmax, n1;

  /* Provide enough storage for scaling factors */
  if (tau->dim < refl)
    resize_avector(tau, refl);
  tauv = tau->v;

  /* Initialize column pivots */
  if (colpiv) {
    for (j = 0; j < cols; j++)
      colpiv[j] = j;
  }

  /* Initialize norm estimate for relative error criterion */
  firstnorm2 = 0.0;

  for (k = 0; k < refl; k++) {
    /* Compute norm of k-th column */
    n1 = rows - k;
    norm2 = REAL_SQR(h2_nrm2(n1, aa + k + k * lda, i_one));
    maxnorm2 = norm2;
    frobnorm2 = norm2;
    jmax = k;

    /* Find maximal norm */
    for (j = k + 1; j < cols; j++) {
      norm2 = REAL_SQR(h2_nrm2(n1, aa + k + j * lda, i_one));
      if (norm2 > maxnorm2) {
	maxnorm2 = norm2;
	jmax = j;
      }
      frobnorm2 += norm2;
    }

    /* Swap columns */
    if (jmax != k) {
      h2_swap(rows, aa + k * lda, i_one, aa + jmax * lda, i_one);

      if (colpiv) {
	j = colpiv[k];
	colpiv[k] = jmax;
	colpiv[jmax] = j;
      }
    }

    /* Prepare norm */
    norm2 = maxnorm2;
    norm = REAL_SQRT(norm2);

    /* Exit if the norm is small enough */
    if (tm && tm->absolute) {
      /* Check the Frobenius norm */
      if (frobnorm2 <= eps * eps)
	break;
    }
    else {
      if (k == 0)
	/* Compute Frobenius norm or estimate spectral norm */
	firstnorm2 = ((tm && tm->frobenius) ? frobnorm2 : norm2);
      else
	/* Compare Frobenius norm with estimate for the entire matrix */
      if (frobnorm2 <= eps * eps * firstnorm2)
	break;
    }

    if (norm2 == 0.0)
      tauv[k] = 0.0;
    else {
      /* Determine reflection vector v */
      diag = aa[k + k * lda];
      alpha = -SIGN1(diag) * norm;

      /* Compute norm of v */
      beta = 1.0 / (norm2 - CONJ(alpha) * diag);

      /* Rescale to ensure v_1 = 1 */
      beta *= ABSSQR(diag - alpha);
      gamma = 1.0 / (diag - alpha);
      n1 = rows - k - 1;
      h2_scal(n1, &gamma, aa + k + 1 + k * lda, i_one);
      tauv[k] = beta;

      /* Compute k-th column */
      aa[k + k * lda] = alpha;

      /* Update remaining columns */
      for (j = k + 1; j < cols; j++) {
	n1 = rows - k - 1;
	field tmp;
	h2_dot(n1, aa + k + 1 + k * lda, i_one, aa + k + 1 + j * lda, i_one, &tmp);
	gamma = aa[k + j * lda] + tmp;

	gamma *= -beta;

	aa[k + j * lda] += gamma;
	h2_axpy(n1, &gamma, aa + k + 1 + k * lda, i_one, aa + k + 1 + j * lda, i_one);
      }
    }
  }

  shrink_avector(tau, k);

  return k;
}
#else
int
qrdecomp_pivot_amatrix(pamatrix a, pavector tau, int * colpiv)
{
  pfield    aa = a->a;
  int      lda = a->ld;
  pfield    tauv;
  int      rows = a->rows;
  int      cols = a->cols;
  int      refl = UINT_MIN(rows, cols);
  field     alpha, beta, gamma, diag;
  double      norm2, norm, maxnorm2;
  int      i, j, k, jmax;

  /* Provide enough storage for scaling factors */
  if (tau->dim < refl)
    resize_avector(tau, refl);
  tauv = tau->v;

  /* Initialize column pivots */
  if (colpiv) {
    for (j = 0; j < cols; j++)
      colpiv[j] = j;
  }

  for (k = 0; k < refl; k++) {
    /* Compute norm of k-th column */
    norm2 = 0.0;
    for (i = k; i < rows; i++)
      norm2 += ABSSQR(aa[i + k * lda]);
    maxnorm2 = norm2;
    jmax = k;

    /* Find maximal norm */
    for (j = k + 1; j < cols; j++) {
      norm2 = 0.0;
      for (i = k; i < rows; i++)
	norm2 += ABSSQR(aa[i + j * lda]);
      if (norm2 > maxnorm2) {
	maxnorm2 = norm2;
	jmax = j;
      }
    }

    /* Swap columns */
    if (jmax != k) {
      for (i = 0; i < rows; i++) {
	alpha = aa[i + k * lda];
	aa[i + k * lda] = aa[i + jmax * lda];
	aa[i + jmax * lda] = alpha;
      }

      if (colpiv) {
	j = colpiv[k];
	colpiv[k] = jmax;
	colpiv[jmax] = j;
      }
    }

    /* Remember pivot */
    if (colpiv)
      colpiv[k] = jmax;

    /* Prepare norm */
    norm2 = maxnorm2;
    norm = REAL_SQRT(norm2);

    if (norm2 == 0.0)
      tauv[k] = 0.0;
    else {
      /* Determine reflection vector v */
      diag = aa[k + k * lda];
      alpha = -SIGN1(diag) * norm;

      /* Compute norm of v */
      beta = 1.0 / (norm2 - CONJ(alpha) * diag);

      /* Rescale to ensure v_1 = 1 */
      beta *= ABSSQR(diag - alpha);
      gamma = 1.0 / (diag - alpha);
      for (i = k + 1; i < rows; i++)
	aa[i + k * lda] *= gamma;
      tauv[k] = beta;

      /* Compute k-th column */
      aa[k + k * lda] = alpha;

      /* Update remaining columns */
      for (j = k + 1; j < cols; j++) {
	gamma = aa[k + j * lda];
	for (i = k + 1; i < rows; i++)
	  gamma += CONJ(aa[i + k * lda]) * aa[i + j * lda];

	gamma *= beta;

	aa[k + j * lda] -= gamma;
	for (i = k + 1; i < rows; i++)
	  aa[i + j * lda] -= gamma * aa[i + k * lda];
      }
    }
  }

  return refl;
}

int
qrdecomp_rank_amatrix(pamatrix a, pavector tau, pctruncmode tm, double eps,
		      int * colpiv)
{
  pfield    aa = a->a;
  int      lda = a->ld;
  pfield    tauv;
  int      rows = a->rows;
  int      cols = a->cols;
  int      refl = UINT_MIN(rows, cols);
  field     alpha, beta, gamma, diag;
  double      norm2, norm, maxnorm2, frobnorm2, firstnorm2;
  int      i, j, k, jmax;

  /* Provide enough storage for scaling factors */
  if (tau->dim < refl)
    resize_avector(tau, refl);
  tauv = tau->v;

  /* Initialize column pivots */
  if (colpiv) {
    for (j = 0; j < cols; j++)
      colpiv[j] = j;
  }

  /* Initialize norm estimate for relative error criterion */
  firstnorm2 = 0.0;

  for (k = 0; k < refl; k++) {
    /* Compute norm of k-th column */
    norm2 = 0.0;
    for (i = k; i < rows; i++)
      norm2 += ABSSQR(aa[i + k * lda]);
    maxnorm2 = norm2;
    frobnorm2 = norm2;
    jmax = k;

    /* Find maximal norm */
    for (j = k + 1; j < cols; j++) {
      norm2 = 0.0;
      for (i = k; i < rows; i++)
	norm2 += ABSSQR(aa[i + j * lda]);
      if (norm2 > maxnorm2) {
	maxnorm2 = norm2;
	jmax = j;
      }
      frobnorm2 += norm2;
    }

    /* Swap columns */
    if (jmax != k) {
      for (i = 0; i < rows; i++) {
	alpha = aa[i + k * lda];
	aa[i + k * lda] = aa[i + jmax * lda];
	aa[i + jmax * lda] = alpha;
      }

      if (colpiv) {
	j = colpiv[k];
	colpiv[k] = jmax;
	colpiv[jmax] = j;
      }
    }

    /* Prepare norm */
    norm2 = maxnorm2;
    norm = REAL_SQRT(norm2);

    /* Exit if the norm is small enough */
    if (tm && tm->absolute) {
      /* Check the Frobenius norm */
      if (frobnorm2 <= eps * eps)
	break;
    }
    else {
      if (k == 0)
	/* Compute Frobenius norm or estimate spectral norm */
	firstnorm2 = ((tm && tm->frobenius) ? frobnorm2 : norm2);
      else
	/* Compare Frobenius norm with estimate for the entire matrix */
      if (frobnorm2 <= eps * eps * firstnorm2)
	break;
    }

    if (norm2 == 0.0)
      tauv[k] = 0.0;
    else {
      /* Determine reflection vector v */
      diag = aa[k + k * lda];
      alpha = -SIGN1(diag) * norm;

      /* Compute norm of v */
      beta = 1.0 / (norm2 - CONJ(alpha) * diag);

      /* Rescale to ensure v_1 = 1 */
      beta *= ABSSQR(diag - alpha);
      gamma = 1.0 / (diag - alpha);
      for (i = k + 1; i < rows; i++)
	aa[i + k * lda] *= gamma;
      tauv[k] = beta;

      /* Compute k-th column */
      aa[k + k * lda] = alpha;

      /* Update remaining columns */
      for (j = k + 1; j < cols; j++) {
	gamma = aa[k + j * lda];
	for (i = k + 1; i < rows; i++)
	  gamma += CONJ(aa[i + k * lda]) * aa[i + j * lda];

	gamma *= beta;

	aa[k + j * lda] -= gamma;
	for (i = k + 1; i < rows; i++)
	  aa[i + j * lda] -= gamma * aa[i + k * lda];
      }
    }
  }

  shrink_avector(tau, k);

  return k;
}
#endif

#ifdef USE_BLAS
/* Remark: if compiled the wrong way, DORMQR, and by extension DORMBR
 * and DGESVD, are currently not thread-safe.
 * gfortran does the right thing if called with "-frecursive", but this
 * appears not to be the standard in, e.g., OpenSUSE Linux. */
#if defined(THREADSAFE_LAPACK) || !defined(USE_OPENMP)
void
qreval_amatrix_avector(bool qtrans, pcamatrix a, pcavector tau, pavector x)
{
  int      rows = a->rows;
  int      cols = a->cols;
  int      refl;
  int		info;

  refl = UINT_MIN3(rows, cols, tau->dim);

  if (refl < 1)
    return;

  assert(x->dim >= rows);
  char left = 'L';
  char trans;

  if (qtrans) {
      trans = 'C';
	  info = h2_ormqr(h2_col_major, left, trans, rows, i_one, refl, (MKL_Complex16 *)(a->a), a->ld,
                   (MKL_Complex16 *)(tau->v), (MKL_Complex16 *)(x->v), x->dim);
    assert(info == 0);
  }
  else {
      trans = 'N';
    info = h2_ormqr(h2_col_major, left, trans, rows, i_one, refl,
	     (MKL_Complex16 *)(a->a), a->ld, (MKL_Complex16 *)(tau->v), (MKL_Complex16 *)(x->v), x->dim);
    assert(info == 0);
  }
}

void
qreval_amatrix(bool qtrans, pcamatrix a, pcavector tau, pamatrix x)
{
  int      rows = a->rows;
  int      cols = a->cols;
  int      refl;
  int      info;

  refl = UINT_MIN3(rows, cols, tau->dim);

  if (refl < 1 || x->cols < 1)
    return;

  assert(x->rows >= rows);


  if (qtrans) {
    info = h2_ormqr(h2_col_major, h2_left, h2_ctrans, rows, x->cols, refl, (MKL_Complex16 *)(a->a), a->ld, (MKL_Complex16 *)(tau->v), (MKL_Complex16 *)(x->a), x->ld);
    assert(info == 0);
  }
  else {
    info = h2_ormqr(h2_col_major, 'L', 'N', rows, x->cols, refl, (MKL_Complex16 *)(a->a), a->ld, (MKL_Complex16 *)(tau->v), (MKL_Complex16 *)(x->a), x->ld);
	//cout << info << ": wrong!" << endl;
    assert(info == 0);
  }
}
#else
void
qreval_amatrix_avector(bool qtrans, pcamatrix a, pcavector tau, pavector x)
{
  pfield    aa = a->a;
  field     work[1];
  int      lda = a->ld;
  pcfield   tauv = tau->v;
  pfield    xv = x->v;
  int      rows = a->rows;
  int      cols = a->cols;
  int      refl, rows1;
  field     beta, gamma;
  int      k;
  char      side='L';

  refl = UINT_MIN3(rows, cols, tau->dim);

  assert(x->dim >= rows);

  if (qtrans) {
    for (k = 0; k < refl; k++) {
      rows1 = rows - k;
      beta = CONJ(tauv[k]);
      gamma = aa[k + k * lda];
      aa[k + k * lda] = 1.0;
      zlarf(&side, &rows1, &i_one, (MKL_Complex16 *)aa + k + k * lda, &i_one,
            (MKL_Complex16 *)&beta, (MKL_Complex16 *)xv + k, &rows, (MKL_Complex16 *)work);
      aa[k + k * lda] = gamma;
    }
  }
  else {
    for (k = refl; k-- > 0;) {
      rows1 = rows - k;
      beta = tauv[k];
      gamma = aa[k + k * lda];
      aa[k + k * lda] = 1.0;
      h2_larf(&side, &rows1, &i_one, (MKL_Complex16 *)aa + k + k * lda, &i_one,
      (MKL_Complex16 *)&beta, (MKL_Complex16 *)xv + k, &rows, (MKL_Complex16 *)work);
      aa[k + k * lda] = gamma;
    }
  }
}

void
qreval_amatrix(bool qtrans, pcamatrix a, pcavector tau, pamatrix x)
{
  pfield    aa = a->a;
  pfield    work;
  int      lda = a->ld;
  pcfield   tauv = tau->v;
  pfield    xa = x->a;
  int      ldx = x->ld;
  int      rows = a->rows;
  int      cols = a->cols;
  int      refl, rows1;
  field     beta, gamma;
  int      k;
  char      side='L';

  refl = UINT_MIN3(rows, cols, tau->dim);

  assert(x->rows >= rows);

  work = allocfield(UINT_MAX(x->rows, x->cols));

  if (qtrans) {
    for (k = 0; k < refl; k++) {
      rows1 = rows - k;
      beta = CONJ(tauv[k]);
      gamma = aa[k + k * lda];
      aa[k + k * lda] = 1.0;
      h2_larf(&side, &rows1, &x->cols, (MKL_Complex16 *)aa + k + k * lda, &i_one, (MKL_Complex16 *)&beta,
	      (MKL_Complex16 *)xa + k, &ldx, (MKL_Complex16 *)work);
      aa[k + k * lda] = gamma;
    }
  }
  else {
    for (k = refl; k-- > 0;) {
      rows1 = rows - k;
      beta = tauv[k];
      gamma = aa[k + k * lda];
      aa[k + k * lda] = 1.0;
      h2_larf(&side, &rows1, &x->cols, (MKL_Complex16 *)aa + k + k * lda, &i_one, (MKL_Complex16 *)&beta,
	      (MKL_Complex16 *)xa + k, &ldx, (MKL_Complex16 *)work);
      aa[k + k * lda] = gamma;
    }
  }

  freemem(work);
}
#endif
#else
void
qreval_amatrix_avector(bool qtrans, pcamatrix a, pcavector tau, pavector x)
{
  pcfield   aa = a->a;
  int      lda = a->ld;
  pcfield   tauv = tau->v;
  pfield    xv = x->v;
  int      rows = a->rows;
  int      cols = a->cols;
  int      refl;
  field     beta, gamma;
  int      i, k;

  refl = UINT_MIN3(rows, cols, tau->dim);

  assert(x->dim >= rows);

  if (qtrans) {
    for (k = 0; k < refl; k++) {
      beta = CONJ(tauv[k]);

      if (beta != 0.0) {
	gamma = xv[k];
	for (i = k + 1; i < rows; i++)
	  gamma += CONJ(aa[i + k * lda]) * xv[i];

	gamma *= beta;

	xv[k] -= gamma;
	for (i = k + 1; i < rows; i++)
	  xv[i] -= gamma * aa[i + k * lda];
      }
    }
  }
  else {
    for (k = refl; k-- > 0;) {
      beta = tauv[k];

      if (beta != 0.0) {
	gamma = xv[k];
	for (i = k + 1; i < rows; i++)
	  gamma += CONJ(aa[i + k * lda]) * xv[i];

	gamma *= beta;

	xv[k] -= gamma;
	for (i = k + 1; i < rows; i++)
	  xv[i] -= gamma * aa[i + k * lda];
      }
    }
  }
}

void
qreval_amatrix(bool qtrans, pcamatrix a, pcavector tau, pamatrix x)
{
  pcfield   aa = a->a;
  int      lda = a->ld;
  pcfield   tauv = tau->v;
  pfield    xa = x->a;
  int      ldx = x->ld;
  int      rows = a->rows;
  int      cols = a->cols;
  int      refl;
  field     beta, gamma;
  int      i, j, k;

  refl = UINT_MIN3(rows, cols, tau->dim);

  assert(x->rows >= rows);

  if (qtrans) {
    for (k = 0; k < refl; k++) {
      beta = CONJ(tauv[k]);

      if (beta != 0.0)
	for (j = 0; j < x->cols; j++) {
	  gamma = xa[k + j * ldx];
	  for (i = k + 1; i < rows; i++)
	    gamma += CONJ(aa[i + k * lda]) * xa[i + j * ldx];

	  gamma *= beta;

	  xa[k + j * ldx] -= gamma;
	  for (i = k + 1; i < rows; i++)
	    xa[i + j * ldx] -= gamma * aa[i + k * lda];
	}
    }
  }
  else {
    for (k = refl; k-- > 0;) {
      beta = tauv[k];

      if (beta != 0.0)
	for (j = 0; j < x->cols; j++) {
	  gamma = xa[k + j * ldx];
	  for (i = k + 1; i < rows; i++)
	    gamma += CONJ(aa[i + k * lda]) * xa[i + j * ldx];

	  gamma *= beta;

	  xa[k + j * ldx] -= gamma;
	  for (i = k + 1; i < rows; i++)
	    xa[i + j * ldx] -= gamma * aa[i + k * lda];
	}
    }
  }
}
#endif

void
qrsolve_amatrix_avector(pcamatrix a, pcavector tau, pavector x)
{
  qreval_amatrix_avector(true, a, tau, x);
  triangularsolve_amatrix_avector(false, false, false, a, x);
}

void
qrinvert_amatrix(pamatrix a)
{
  pamatrix  acopy;
  pavector  tau, b;
  amatrix   atmp;
  avector   ttmp, btmp;
  int      n = a->rows;
  int      i;

  assert(n == a->cols);

  acopy = init_amatrix(&atmp, n, n);
  tau = init_avector(&ttmp, n);

  copy_amatrix(false, a, acopy);

  qrdecomp_amatrix(acopy, tau);

  identity_amatrix(a);

  for (i = 0; i < n; i++) {
    b = init_row_avector(&btmp, a, i);
    qrsolve_amatrix_avector(acopy, tau, b);
    uninit_avector(b);
  }

  uninit_avector(tau);
  uninit_amatrix(acopy);
}

#ifdef USE_BLAS
void
qrexpand_amatrix(pcamatrix a, pcavector tau, pamatrix q)
{
  field    *work;
  int      refl, lwork, info;

  refl = UINT_MIN(UINT_MIN(q->cols, tau->dim), UINT_MIN(a->rows, a->cols));

  /* Quick exit if no reflections used */
  if (refl == 0) {
    identity_amatrix(q);
    return;
  }

  copy_sub_amatrix(false, a, q);

  lwork = 4 * a->rows;
  work = allocfield(lwork);

  info = h2_orgqr(h2_col_major, q->rows, q->cols, refl, (MKL_Complex16 *)(q->a), q->ld, (MKL_Complex16 *)(tau->v));
  assert(info == 0);

  freemem(work);
}
#else
void
qrexpand_amatrix(pcamatrix a, pcavector tau, pamatrix q)
{
  pcfield   aa = a->a;
  int      lda = a->ld;
  pcfield   tauv = tau->v;
  pfield    qa = q->a;
  int      ldq = q->ld;
  int      rows = a->rows;
  int      refl;
  field     beta, gamma;
  int      i, j, k;

  /* Determine number of relevant elementary reflections */
  refl = UINT_MIN(UINT_MIN(q->cols, tau->dim), UINT_MIN(a->rows, a->cols));

  assert(q->rows >= rows);
  assert(q->cols <= rows);

  /* Create identity matrix */
  for (j = 0; j < q->cols; j++) {
    for (i = 0; i < rows; i++)
      qa[i + j * ldq] = 0.0;
    qa[j + j * ldq] = 1.0;
  }

  /* Apply reflections in reversed order */
  for (k = refl; k-- > 0;) {
    beta = tauv[k];
    if (beta != 0.0)
      for (j = k; j < q->cols; j++) {
	gamma = qa[k + j * ldq];
	for (i = k + 1; i < rows; i++)
	  gamma += CONJ(aa[i + k * lda]) * qa[i + j * ldq];

	gamma *= beta;

	qa[k + j * ldq] -= gamma;
	for (i = k + 1; i < rows; i++)
	  qa[i + j * ldq] -= gamma * aa[i + k * lda];
      }
  }
}
#endif
