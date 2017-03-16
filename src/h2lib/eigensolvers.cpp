
/* ------------------------------------------------------------
 * This is the file "eigensolvers.c" of the H2Lib package.
 * All rights reserved, Steffen Boerm 2009
 * ------------------------------------------------------------ */

#include "eigensolvers.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "h2basic.h"
#include "factorizations.h"

/** Relative accuracy used in the stopping criterion of the QR iteration. */
#define H2_QR_EPS 1e-14

/** Bound for determining when a number is essentially zero. */
#define H2_ALMOST_ZERO 1e-300

/** Run-time checks for self-made solvers. */
/* #define RUNTIME_CHECK_EIGENSOLVERS */

/** Relative tolerance for run-time checks. */
#ifdef USE_FLOAT
#define H2_CHECK_TOLERANCE 1.0e-6
#else
#define H2_CHECK_TOLERANCE 1.0e-12
#endif

/* ------------------------------------------------------------
 * Constructors and destructors
 * ------------------------------------------------------------ */

ptridiag
init_tridiag(ptridiag T, int size)
{
  T->d = allocreal(size);
  T->u = (size > 1 ? allocreal(size - 1) : NULL);
  T->l = (size > 1 ? allocreal(size - 1) : NULL);
  T->size = size;
  T->owner = NULL;

  return T;
}

ptridiag
init_sub_tridiag(ptridiag T, ptridiag src, int size, int off)
{
  assert(off + size <= src->size);

  T->d = src->d + off;
  T->u = (size > 1 ? src->u + off : NULL);
  T->l = (size > 1 ? src->l + off : NULL);
  T->size = size;
  T->owner = src;

  return T;
}

ptridiag
init_vec_tridiag(ptridiag T, prealavector src, int size)
{
  assert(3 * size - 2 <= src->dim);

  T->d = src->v;
  T->u = (size > 1 ? src->v + size : NULL);
  T->l = (size > 1 ? src->v + 2 * size - 1 : NULL);
  T->size = size;
  T->owner = (ptridiag) src;

  return T;
}

void
uninit_tridiag(ptridiag T)
{
  if (T->owner == NULL) {
    freemem(T->d);
    if (T->size > 1) {
      freemem(T->l);
      freemem(T->u);
    }
    else {
      assert(T->l == NULL);
      assert(T->u == NULL);
    }
  }
}

ptridiag
new_tridiag(int size)
{
  ptridiag  T;

  T = (ptridiag) allocmem(sizeof(tridiag));

  init_tridiag(T, size);

  return T;
}

void
del_tridiag(ptridiag T)
{
  uninit_tridiag(T);

  freemem(T);
}

void
copy_tridiag(pctridiag T, ptridiag Tcopy)
{
  int      size = T->size;
  double*     d = T->d;
  double*     u = T->u;
  double*     l = T->l;
  int      i;

  assert(Tcopy->size >= T->size);

  for (i = 0; i < size; i++)
    Tcopy->d[i] = d[i];
  for (i = 0; i < size - 1; i++) {
    Tcopy->u[i] = u[i];
    Tcopy->l[i] = l[i];
  }
}

ptridiag
clone_tridiag(pctridiag T)
{
  ptridiag  Tcopy;

  Tcopy = new_tridiag(T->size);
  copy_tridiag(T, Tcopy);

  return Tcopy;
}

double
check_tridiag(pctridiag T, pcamatrix Ts)
{
  int      n = T->size;
  const double *    d = T->d;
  const double *    l = T->l;
  const double *    u = T->u;
  pcfield   a = Ts->a;
  int      lda = Ts->ld;
  int      i, j;
  field     val;
  double      norm;

  assert(Ts->rows >= n);
  assert(Ts->cols >= n);

  norm = 0.0;
  for (j = 0; j < n; j++)
    for (i = 0; i < n; i++) {
      val = a[i + j * lda];
      if (i == j)
	val -= d[i];
      else if (i == j + 1)
	val -= l[j];
      else if (i + 1 == j)
	val -= u[i];
      norm += ABSSQR(val);
    }

  return REAL_SQRT(norm);
}

double
check_lower_tridiag(pctridiag T, pcamatrix Ts)
{
  int      n = T->size;
  const double*    d = T->d;
  const double*    l = T->l;
  pcfield   a = Ts->a;
  int      lda = Ts->ld;
  int      i, j;
  field     val;
  double      norm;

  assert(Ts->rows >= n);
  assert(Ts->cols >= n);

  norm = 0.0;
  for (j = 0; j < n; j++)
    for (i = 0; i < n; i++) {
      val = a[i + j * lda];
      if (i == j)
	val -= d[i];
      else if (i == j + 1)
	val -= l[j];

      if (ABS(val) > 1e-10)
	printf("%u %u %g\n", i, j, ABS(val));
      norm += ABSSQR(val);
    }

  return REAL_SQRT(norm);
}

#ifdef USE_BLAS
void
diageval_tridiag_amatrix(field alpha, bool atrans, pctridiag a,
			 bool xtrans, pamatrix x)
{
  const double*     d = a->d;
  int      n = a->size;
  pfield    xa = x->a;
  int      ldx = x->ld;
  field     beta;
  int  i;

  if (n < 1)			/* Quick exit */
    return;

  if (xtrans) {
    assert(x->cols == a->size);

    for (i = 0; i < n; i++) {
      beta = (atrans ? CONJ(alpha * d[i]) : CONJ(alpha) * d[i]);

      h2_scal(x->rows, &beta, xa + i * ldx, i_one);
    }
  }
  else {
    assert(x->rows == a->size);

    for (i = 0; i < n; i++) {
      beta = (atrans ? alpha * d[i] : alpha * d[i]);

      h2_scal(x->cols, &beta, xa + i, ldx);
    }
  }
}
#else
void
diageval_tridiag_amatrix(field alpha,
			 bool atrans, pctridiag a, bool xtrans, pamatrix x)
{
  const double*     d = a->d;
  int      n = a->size;
  pfield    xa = x->a;
  int      ldx = x->ld;
  field     beta;
  int  i, j;

  if (n < 1)			/* Quick exit */
    return;

  if (xtrans) {
    assert(x->cols == a->size);

    for (i = 0; i < n; i++) {
      beta = (atrans ? CONJ(alpha * d[i]) : CONJ(alpha) * d[i]);

      for (j = 0; j < x->rows; j++)
	xa[j + i * ldx] *= beta;
    }
  }
  else {
    assert(x->rows == a->size);

    for (i = 0; i < n; i++) {
      beta = (atrans ? alpha * d[i] : alpha * d[i]);

      for (j = 0; j < x->cols; j++)
	xa[i + j * ldx] *= beta;
    }
  }
}
#endif

#ifdef USE_BLAS
void
lowereval_tridiag_amatrix(field alpha, bool atrans, pctridiag a,
			  bool xtrans, pamatrix x)
{
  const double*     d = a->d;
  const double*     l = a->l;
  int      n = a->size;
  pfield    xa = x->a;
  int      ldx = x->ld;
  field     beta, gamma;
  int      i;

  if (n < 1)			/* Quick exit */
    return;

  if (atrans) {
    if (xtrans) {
      assert(x->cols == a->size);

      for (i = 0; i < n - 1; i++) {
	beta = CONJ(alpha) * d[i];
	gamma = CONJ(alpha) * l[i];

	h2_scal(x->rows, &beta, xa + i * ldx, i_one);
	h2_axpy(x->rows, &gamma, xa + (i + 1) * ldx, i_one, xa + i * ldx,
		i_one);
      }
      beta = CONJ(alpha) * d[n - 1];

      h2_scal(x->rows, &beta, xa + (n - 1) * ldx, i_one);
    }
    else {
      assert(x->rows == a->size);

      for (i = 0; i < n - 1; i++) {
	beta = alpha * d[i];
	gamma = alpha * l[i];

	h2_scal(x->cols, &beta, xa + i, ldx);
	h2_axpy(x->cols, &gamma, xa + (i + 1), ldx, xa + i, ldx);
      }
      beta = alpha * d[n - 1];

      h2_scal(x->cols, &beta, xa + (n - 1), ldx);
    }
  }
  else {
    if (xtrans) {
      assert(x->cols == a->size);

      for (i = n; i-- > 1;) {
	beta = CONJ(alpha * d[i]);
	gamma = CONJ(alpha * l[i - 1]);

	h2_scal(x->rows, &beta, xa + i * ldx, i_one);
	h2_axpy(x->rows, &gamma, xa + (i - 1) * ldx, i_one, xa + i * ldx,
		i_one);
      }
      beta = CONJ(alpha * d[0]);

      h2_scal(x->rows, &beta, xa, i_one);
    }
    else {
      assert(x->rows == a->size);

      for (i = n; i-- > 1;) {
	beta = alpha * d[i];
	gamma = alpha * l[i - 1];

	h2_scal(x->cols, &beta, xa + i, ldx);
	h2_axpy(x->cols, &gamma, xa + (i - 1), ldx, xa + i, ldx);
      }
      beta = alpha * d[0];

      h2_scal(x->cols, &beta, xa, ldx);
    }
  }
}
#else
void
lowereval_tridiag_amatrix(field alpha,
			  bool atrans, pctridiag a, bool xtrans, pamatrix x)
{
  const double*     d = a->d;
  const double*     l = a->l;
  int      n = a->size;
  pfield    xa = x->a;
  int      ldx = x->ld;
  field    beta, gamma;
  int  i, j;

  if (n < 1)			/* Quick exit */
    return;

  if (atrans) {
    if (xtrans) {
      assert(x->cols == a->size);

      for (i = 0; i < n - 1; i++) {
	beta = CONJ(alpha) * d[i];
	gamma = CONJ(alpha) * l[i];

	for (j = 0; j < x->rows; j++)
	  xa[j + i * ldx] =
	    beta * xa[j + i * ldx] + gamma * xa[j + (i + 1) * ldx];
      }
      beta = CONJ(alpha) * d[n - 1];

      for (j = 0; j < x->rows; j++)
	xa[j + (n - 1) * ldx] *= beta;
    }
    else {
      assert(x->rows == a->size);

      for (i = 0; i < n - 1; i++) {
	beta = alpha * (d[i]);
	gamma = alpha * (l[i]);

	for (j = 0; j < x->cols; j++)
	  xa[i + j * ldx] =
	    beta * xa[i + j * ldx] + gamma * xa[(i + 1) + j * ldx];
      }
      beta = alpha * (d[n - 1]);

      for (j = 0; j < x->cols; j++)
	xa[(n + 1) + j * ldx] *= beta;
    }
  }
  else {
    if (xtrans) {
      assert(x->cols == a->size);

      for (i = n; i-- > 1;) {
	beta = CONJ(alpha * d[i]);
	gamma = CONJ(alpha * l[i - 1]);

	for (j = 0; j < x->rows; j++)
	  xa[j + i * ldx] =
	    beta * xa[j + i * ldx] + gamma * xa[j + (i - 1) * ldx];
      }
      beta = CONJ(alpha * d[0]);

      for (j = 0; j < x->rows; j++)
	xa[j] *= beta;
    }
    else {
      assert(x->rows == a->size);

      for (i = n; i-- > 1;) {
	beta = alpha * d[i];
	gamma = alpha * l[i - 1];

	for (j = 0; j < x->cols; j++)
	  xa[i + j * ldx] =
	    beta * xa[i + j * ldx] + gamma * xa[(i - 1) + j * ldx];
      }
      beta = alpha * d[0];

      for (j = 0; j < x->cols; j++)
	xa[j * ldx] *= beta;
    }
  }
}
#endif

/* ------------------------------------------------------------
 * Givens rotation
 * ------------------------------------------------------------ */

static void
givens(double a, double b, double* c, double* s)
{
  double      t;

  if (a == 0.0 && b == 0.0) {
    *c = 1.0;
    *s = 0.0;
  }
  else if (REAL_ABS(a) > REAL_ABS(b)) {
    t = b / a;
    *c = REAL_RSQRT(1.0 + REAL_SQR(t));
    *s = t * (*c);
  }
  else {
    t = a / b;
    *s = REAL_RSQRT(1.0 + REAL_SQR(t));
    *c = t * (*s);
  }
}

/* ------------------------------------------------------------
 * One QR step for the symmetric tridiagonal matrix
 * ------------------------------------------------------------ */

void
qrstep_tridiag(ptridiag T, field shift, pamatrix Q)
{
  double      a0, a1, b0, b0t, carry;
  double      c, s;
  field     x, y;
  double*     a = T->d;
  double*     b = T->l;
  pfield    qa;
  int      ldq;
  int      n = T->size;
  int      i, j;

  assert(Q == NULL || Q->cols == n);

  qa = (Q ? Q->a : NULL);
  ldq = (Q ? Q->ld : 0);

  /* Quick exit if trivial matrix */
  if (n < 2)
    return;

  /* Determine Givens rotation */
  givens((a[0] - shift).real(), b[0], &c, &s);

  /* Apply to first two rows */
  a0 = c * a[0] + s * b[0];
  b0 = -s * a[0] + c * b[0];
  b0t = c * b[0] + s * a[1];
  a1 = -s * b[0] + c * a[1];

  /* Apply to first two columns */
  a[0] = c * a0 + s * b0t;
  b[0] = c * b0 + s * a1;
  a[1] = -s * b0 + c * a1;

  /* Apply to first two columns of Q */
  if (Q)
    for (j = 0; j < Q->rows; j++) {
      x = qa[j];
      y = qa[j + ldq];
      qa[j] = c * x + s * y;
      qa[j + ldq] = -s * x + c * y;
    }

  /* Chase away coefficient at (2,0) */
  for (i = 1; i < n - 1; i++) {
    /* Apply Givens rotation to rows i, i+1 */
    carry = s * b[i];
    b[i] = c * b[i];

    /* Determine Givens rotation to eliminate carry... */
    givens(b[i - 1], carry, &c, &s);

    /* ... and apply it */
    b[i - 1] = c * b[i - 1] + s * carry;

    /* Apply it also to remainder of rows i, i+1, */
    a0 = c * a[i] + s * b[i];
    b0 = -s * a[i] + c * b[i];
    b0t = c * b[i] + s * a[i + 1];
    a1 = -s * b[i] + c * a[i + 1];

    /* apply it to the columns i, i+1, */
    a[i] = c * a0 + s * b0t;
    b[i] = c * b0 + s * a1;
    a[i + 1] = -s * b0 + c * a1;

    /* and if necessary to the columns i, i+1 of Q */
    if (Q)
      for (j = 0; j < Q->rows; j++) {
	x = qa[j + i * ldq];
	y = qa[j + (i + 1) * ldq];
	qa[j + i * ldq] = c * x + s * y;
	qa[j + (i + 1) * ldq] = -s * x + c * y;
      }
  }
}

/* ------------------------------------------------------------
 * QR iteration for self-adjoint tridiagonal matrices
 * ------------------------------------------------------------ */

typedef struct {
  ptridiag  T;
  pamatrix  U, Vt;
} evpsortdata;

static bool
evp_leq(int i, int j, void *data)
{
  evpsortdata *esd = (evpsortdata *) data;
  ptridiag  T = esd->T;

  assert(i < T->size);
  assert(j < T->size);

  return (T->d[i] <= T->d[j]);
}

static bool
evp_geq(int i, int j, void *data)
{
  evpsortdata *esd = (evpsortdata *) data;
  ptridiag  T = esd->T;

  assert(i < T->size);
  assert(j < T->size);

  return (T->d[i] >= T->d[j]);
}

static void
evp_swap(int i, int j, void *data)
{
  evpsortdata *esd = (evpsortdata *) data;
  ptridiag  T = esd->T;
  pamatrix  U = esd->U;
  pamatrix  Vt = esd->Vt;
  int      ldu, ldv;
  double      h;
  field     hf;
  int      k;

  assert(i < T->size);
  assert(j < T->size);
  assert(U == NULL || U->cols >= T->size);
  assert(Vt == NULL || Vt->rows >= T->size);

  h = T->d[i];
  T->d[i] = T->d[j];
  T->d[j] = h;

  if (U) {
    ldu = U->ld;
    for (k = 0; k < U->rows; k++) {
      hf = U->a[k + i * ldu];
      U->a[k + i * ldu] = U->a[k + j * ldu];
      U->a[k + j * ldu] = hf;
    }
  }

  if (Vt) {
    ldv = Vt->ld;
    for (k = 0; k < Vt->cols; k++) {
      hf = Vt->a[i + k * ldv];
      Vt->a[i + k * ldv] = Vt->a[j + k * ldv];
      Vt->a[j + k * ldv] = hf;
    }
  }
}

int
sb_muleig_tridiag(ptridiag T, pamatrix Q, int maxiter)
{
  tridiag   tmp;
  amatrix   tmp2;
  evpsortdata esd;
  ptridiag  Tsub;
  pamatrix  Qsub;
  double*     a = T->d;
  double*     b = T->l;
  int      n = T->size;
  double      d, m, a0, a1, b0, db, lambda, lambda1;
  int      off, size, iter;

  if (n < 1)
    return 0;

  /* Find first non-zero subdiagonal coefficient */
  off = 0;
  while (off < n - 1 &&
	 REAL_ABS(b[off]) <
	 H2_QR_EPS * (REAL_ABS(a[off]) + REAL_ABS(a[off + 1]))) {
    b[off] = 0.0;
    off++;
  }

  iter = 0;
  while (off < n - 1 && (maxiter == 0 || iter < maxiter)) {
    iter++;

    /* Find dimension of unreduced submatrix */
    size = 2;
    while (size < n - off
	   &&
	   REAL_ABS(b[off + size - 1]) >= H2_QR_EPS
	   * (REAL_ABS(a[off + size - 1]) + REAL_ABS(a[off + size])))
      size++;

    /* Determine Wilkinson shift */
    a0 = a[off + size - 2];
    a1 = a[off + size - 1];
    b0 = b[off + size - 2];
    m = (a0 + a1) * 0.5;
    d = (a0 - a1) * 0.5;
    db = (REAL_ABS(d) > REAL_ABS(b0) ?
	  REAL_ABS(d) * REAL_SQRT(1.0 + REAL_SQR(b0 / d)) :
	  REAL_ABS(b0) * REAL_SQRT(1.0 + REAL_SQR(d / b0)));
    lambda1 = (m > 0.0 ? m + db : m - db);
    lambda = (m * d < 0.0 ? lambda1 : (a0 * a1 - ABSSQR(b0)) / lambda1);

    /* Set up submatrices */
    Tsub = init_sub_tridiag(&tmp, T, size, off);
    Qsub = (Q ? init_sub_amatrix(&tmp2, Q, Q->rows, 0, size, off) : NULL);

    /* Perform QR step on submatrices */
    qrstep_tridiag(Tsub, lambda, Qsub);

    /* Clean up submatrices */
    if (Qsub)
      uninit_amatrix(Qsub);
    uninit_tridiag(Tsub);

    /* Find non-zero subdiagonal coefficient */
    while (off < n - 1 &&
	   REAL_ABS(b[off]) <
	   H2_QR_EPS * (REAL_ABS(a[off]) + REAL_ABS(a[off + 1]))) {
      b[off] = 0.0;
      off++;
    }
  }

  if (maxiter > 0 && iter == maxiter) {
    (void)
      printf("  ### Eigensolver warning: did not converge after %u steps\n",
	     iter);
  }

  /* Sort eigenvalues */
  esd.T = T;
  esd.U = Q;
  esd.Vt = 0;
  heapsort(n, evp_leq, evp_swap, &esd);

  return iter;
}

#ifdef USE_BLAS
int
muleig_tridiag(ptridiag T, pamatrix Q)
{
  double*     work;
  int       info = 0;

  if (T->size > 1) {
    if (Q) {
      work = allocreal(2 * T->size + 2);
      info = h2_steqr(h2_col_major, 'N', T->size, T->d, T->l, (MKL_Complex16 *)(Q->a), Q->ld);
      assert(info >= 0);
      freemem(work);
    }
    else {
      info = h2_stev(h2_col_major, 'N', T->size, T->d, T->l, NULL, i_one);
      assert(info >= 0);
    }
  }
  return (info != 0);
}
#else
int
muleig_tridiag(ptridiag T, pamatrix Q)
{
  int      iter, maxiter;

  maxiter = 32 * T->size;
  iter = sb_muleig_tridiag(T, Q, maxiter);

  return !(iter < maxiter);
}
#endif

int
eig_tridiag(ptridiag T, pamatrix Q)
{
  if (Q)
    identity_amatrix(Q);

  return muleig_tridiag(T, Q);
}

/* ------------------------------------------------------------
 * Hessenberg tridiagonalization
 * ------------------------------------------------------------ */

void
sb_tridiagonalize_amatrix(pamatrix A, ptridiag T, pamatrix Q)
{
  pfield    aa;
  int      lda;
  pfield    qa;
  int      ldq;
  double     *d, *l, *u;
  int      n;
  int      i, j, k;
  double      norm2, norm;
  field     first, alpha, beta, gamma;
  field     nsign, psign;

  n = A->rows;

  assert(A->cols == n);
  assert(Q->rows == n);
  assert(Q->cols == n);
  assert(T->size == n);

  aa = A->a;
  lda = A->ld;

  d = T->d;
  l = T->l;
  u = T->u;

  qa = (Q ? Q->a : 0);
  ldq = (Q ? Q->ld : 0);

  if (Q)
    identity_amatrix(Q);

  /* Quick exit if trivial matrix */
  if (n < 1)
    return;
  else if (n < 2) {
    d[0] = aa[0].real();
    return;
  }

  /* Initialize sign */
  psign = 1.0;

  for (k = 0; k < n - 1; k++) {
    /* Compute norm of k-th column */
    norm2 = ABSSQR(aa[(k + 1) + k * lda]);
    for (i = k + 2; i < n; i++)
      norm2 += ABSSQR(aa[i + k * lda]);
    norm = REAL_SQRT(norm2);

    if (norm2 <= H2_ALMOST_ZERO) {
      d[k] = REAL(aa[k + k * lda]);
      l[k] = 0.0;
      u[k] = 0.0;
    }
    else {
      /* Determine Householder reflection vector */
      first = aa[(k + 1) + k * lda];
      alpha = -SIGN1(first) * norm;
      gamma = first - alpha;

      /* Compute 2 / |v|^2 */
      beta = 1.0 / (norm2 + ABS(first) * norm);

      /* Normalize reflection vector */
      for (i = k + 2; i < n; i++)
	aa[i + k * lda] /= gamma;
      beta *= ABSSQR(gamma);

      /* Compute k-th column */
      nsign = SIGN1(alpha);
      d[k] = REAL(aa[k + k * lda]);
      l[k] = u[k] = ABS(alpha);

      /* Update remaining columns */
      for (j = k + 1; j < n; j++) {
	gamma = aa[(k + 1) + j * lda];
	for (i = k + 2; i < n; i++)
	  gamma += CONJ(aa[i + k * lda]) * aa[i + j * lda];

	gamma *= beta;

	aa[(k + 1) + j * lda] -= gamma;
	for (i = k + 2; i < n; i++)
	  aa[i + j * lda] -= gamma * aa[i + k * lda];
      }

      /* Update remaining rows */
      for (j = k + 1; j < n; j++) {
	gamma = aa[j + (k + 1) * lda];
	for (i = k + 2; i < n; i++)
	  gamma += aa[i + k * lda] * aa[j + i * lda];

	gamma *= beta;

	aa[j + (k + 1) * lda] -= gamma;
	for (i = k + 2; i < n; i++)
	  aa[j + i * lda] -= gamma * CONJ(aa[i + k * lda]);
      }

      if (Q) {
	/* Apply reflection to all rows of Q */
	for (j = 0; j < n; j++) {
	  gamma = qa[j + (k + 1) * ldq];
	  for (i = k + 2; i < n; i++)
	    gamma += aa[i + k * lda] * qa[j + i * ldq];

	  gamma *= beta;

	  qa[j + (k + 1) * ldq] -= gamma;
	  for (i = k + 2; i < n; i++)
	    qa[j + i * ldq] -= gamma * CONJ(aa[i + k * lda]);
	}

	/* Adjust sign of column k+1 */
	for (j = 0; j < n; j++)
	  qa[j + (k + 1) * ldq] *= nsign * CONJ(psign);
      }

      /* Update sign */
      psign *= CONJ(nsign);
    }
  }
  T->d[k] = REAL(aa[k + k * lda]);
}

#if 0
void
tridiagonalize_amatrix(pamatrix A, ptridiag T, pamatrix Q)
{
  pfield    aa;
  int      lda;
  pfield    qa;
  int      ldq;
  double    *d,*l, *u;
  pfield    work;
  int      n, n1;
  int      k;
  field     alpha, beta;
  field     nsign, psign;

  n = A->rows;

  assert(A->cols == n);
  assert(Q->rows == n);
  assert(Q->cols == n);
  assert(T->size == n);

  aa = A->a;
  lda = A->ld;

  d = T->d;
  l = T->l;
  u = T->u;

  qa = (Q ? Q->a : 0);
  ldq = (Q ? Q->ld : 0);

  work = allocfield(UINT_MAX(A->rows, A->cols));

  if (Q)
    identity_amatrix(Q);

  /* Quick exit if trivial matrix */
  if (n < 1)
    return;
  else if (n < 2) {
    d[0] = REAL(aa[0]);
    return;
  }

  /* Initialize sign */
  psign = 1.0;

  for (k = 0; k < n - 1; k++) {
    /* Determine Householder vector */
    n1 = n - k - 1;
    alpha = aa[(k + 1) + k * lda];
    h2_larfg(n1, (MKL_Complex16 *)&alpha, (MKL_Complex16 *)( aa + (k + 2) + k * lda), i_one, (MKL_Complex16 *)&beta);

    /* Compute k-th column */
    nsign = SIGN1(alpha);
    d[k] = REAL(aa[k + k * lda]);
    l[k] = u[k] = ABS(alpha);

    /* Update remaining columns */
    aa[(k + 1) + k * lda] = 1.0;
    beta = CONJ(beta);
    h2_larf(_h2_left, &n1, &n1, aa + (k + 1) + k * lda, &i_one, &beta,
	    aa + (k + 1) + (k + 1) * lda, &lda, work);

    /* Update remaining rows */
    beta = CONJ(beta);
    h2_larf(_h2_right, &n1, &n1, aa + (k + 1) + k * lda, &i_one, &beta,
	    aa + (k + 1) + (k + 1) * lda, &lda, work);

    if (Q) {
      /* Update Q */
      h2_larf(_h2_right, &n, &n1, aa + (k + 1) + k * lda, &i_one, &beta,
	      qa + (k + 1) * ldq, &ldq, work);

      /* Adjust sign of column k+1 */
      beta = nsign * CONJ(psign);
      h2_scal(n, &beta, qa + (k + 1) * ldq, i_one);

      /* Update sign */
      psign *= CONJ(nsign);
    }
  }
  T->d[k] = REAL(aa[k + k * lda]);

  freemem(work);
}
#else
void
tridiagonalize_amatrix(pamatrix A, ptridiag T, pamatrix Q)
{
  sb_tridiagonalize_amatrix(A, T, Q);
}
#endif

/* ------------------------------------------------------------
 * Self-adjoint eigenvalue problem
 * ------------------------------------------------------------ */

int
sb_eig_amatrix(pamatrix A, prealavector lambda, pamatrix Q, int maxiter)
{
  tridiag   tmp;
  ptridiag  T;
  int      n, iter;
  int      i;

  n = A->rows;

  assert(A->cols == n);
  assert(Q == 0 || Q->rows == n);
  assert(Q == 0 || Q->cols == n);

  /* Set up auxiliary matrix */
  T = init_tridiag(&tmp, n);

  /* Tridiagonalize A */
  sb_tridiagonalize_amatrix(A, T, Q);

  /* Solve tridiagonal eigenproblem */
  iter = sb_muleig_tridiag(T, Q, maxiter);

  /* Copy eigenvalues */
  if (lambda->dim < n) {
    resize_realavector(lambda, n);
  }
  for (i = 0; i < n; i++) {
    lambda->v[i] = T->d[i];
  }

  /* Clean up */
  uninit_tridiag(T);

  return iter;
}

#ifdef USE_BLAS
/* Remark: if compiled the wrong way, DORMQR, and by extension DORMBR
 * and DGESVD, are currently not thread-safe.
 * gfortran does the right thing if called with "-frecursive", but this
 * appears not to be the standard in, e.g., OpenSUSE Linux. */
#if defined(THREADSAFE_LAPACK) || !defined(USE_OPENMP)
int
eig_amatrix(pamatrix A, prealavector lambda, pamatrix Q)
{
  pfield    work;

  double*     rwork;

  int  lwork;
  int      n = A->rows;
  int       info = 0;

  /* Quick exit if trivial matrix */
  if (n < 2) {
    if (Q)
      identity_amatrix(Q);
    return 0;
  }

  lwork = 12 * n;
  work = allocfield(lwork);

  rwork = allocreal(3 * n);

  if (Q) {
    h2_heev(h2_col_major, 'V', h2_lower, n, (MKL_Complex16 *)(A->a), A->ld, lambda->v);
    copy_amatrix(false, A, Q);
  }
  else {
    h2_heev(h2_col_major, 'N', h2_lower, n, (MKL_Complex16 *)(A->a), A->ld, lambda->v);
  }

#ifdef USE_COMPLEX
  freemem(rwork);
#endif

  freemem(work);

  return (info != 0);
}
#else
int
eig_amatrix(pamatrix A, prealavector lambda, pamatrix Q)
{
  tridiag   tmp;
  ptridiag  T;
  int      n, info;
  int      i;

  n = A->rows;

  assert(A->cols == n);
  assert(Q == 0 || Q->rows == n);
  assert(Q == 0 || Q->cols == n);

  /* Set up auxiliary matrix */
  T = init_tridiag(&tmp, n);

  /* Tridiagonalize A */
  tridiagonalize_amatrix(A, T, Q);

  /* Solve tridiagonal eigenproblem */
  info = muleig_tridiag(T, Q);

  /* Copy eigenvalues */
  if (lambda->dim < n) {
    resize_realavector(lambda, n);
  }
  for (i = 0; i < n; i++) {
    lambda->v[i] = T->d[i];
  }

  /* Clean up */
  uninit_tridiag(T);

  return info;
}
#endif
#else
int
eig_amatrix(pamatrix A, prealavector lambda, pamatrix Q)
{
  int      iter, maxiter;

  assert(A->rows == A->cols);

  maxiter = 32 * A->rows;

  iter = sb_eig_amatrix(A, lambda, Q, maxiter);

  return !(iter < maxiter);
}
#endif

int
geig_amatrix(pamatrix A, pamatrix M, prealavector lambda, pamatrix Q)
{
  int      info;

  /* Compute Cholesky factorization M = L L^* */
  choldecomp_amatrix(M);

  /* Compute B = L^{-1} A L^{-*} */
  triangularsolve_amatrix(true, false, false, M, true, A);
  triangularsolve_amatrix(true, false, false, M, false, A);

  /* Solve eigenvalue problem for B */
  info = eig_amatrix(A, lambda, Q);

  /* Compute eigenvector matrix L^{-*} Q */
  triangularsolve_amatrix(true, false, true, M, false, Q);

  return info;
}

/* ------------------------------------------------------------
 * One SVD step for a sub-bidiagonal matrix
 * ------------------------------------------------------------ */

void
svdstep_tridiag(ptridiag T, field shift, pamatrix U, pamatrix Vt)
{
  double      a0, a1, b0, b0t, carry;
  double      c, s, xr, yr;
  pfield    ua, va;
  field     x, y;
  int      ldu, ldv;
  double*     a = T->d;
  double*     b = T->l;
  int      n = T->size;
  int      i, j;

  assert(U == NULL || U->cols == n);
  assert(Vt == NULL || Vt->rows == n);

  ua = (U ? U->a : NULL);
  ldu = (U ? U->ld : 0);

  va = (Vt ? Vt->a : NULL);
  ldv = (Vt ? Vt->ld : 0);

  /* Quick exit if trivial matrix */
  if (n < 2)
    return;

  /* Determine Givens rotation */
  givens((REAL_SQR(a[0]) - shift).real(), b[0] * a[0], &c, &s);

  /* Apply to first two rows */
  a0 = c * a[0] + s * b[0];
  b0 = -s * a[0] + c * b[0];
  b0t = s * a[1];
  a1 = c * a[1];

  /* Apply to first two columns of U */
  if (U)
    for (j = 0; j < U->rows; j++) {
      x = ua[j];
      y = ua[j + ldu];
      ua[j] = c * x + s * y;
      ua[j + ldu] = -s * x + c * y;
    }

  /* Determine Givens rotation to eliminate (0,1) */
  givens(a0, b0t, &c, &s);

  /* Apply to first two columns of rows 0 and 1 */
  a[0] = c * a0 + s * b0t;
  b[0] = c * b0 + s * a1;
  a[1] = -s * b0 + c * a1;

  /* Apply to first two rows of Vt */
  if (Vt)
    for (j = 0; j < Vt->cols; j++) {
      x = va[j * ldv];
      y = va[1 + j * ldv];
      va[j * ldv] = c * x + s * y;
      va[1 + j * ldv] = -s * x + c * y;
    }

  /* Chase away coefficient at (2,0) */
  for (i = 1; i < n - 1; i++) {
    /* Apply Givens rotation to columns i-1, i of row i+1 */
    carry = s * b[i];
    b[i] = c * b[i];

    /* Determine Givens rotation to eliminate carry... */
    givens(b[i - 1], carry, &c, &s);

    /* ... and apply it to rows i,i+1 of column i-1, ... */
    b[i - 1] = c * b[i - 1] + s * carry;

    /* ... rows i,i+1 of column i, ... */
    xr = a[i];
    yr = b[i];
    a[i] = c * xr + s * yr;
    b[i] = -s * xr + c * yr;

    /* ... and rows i,i+1 of column i+1 */
    carry = s * a[i + 1];
    a[i + 1] = c * a[i + 1];

    /* Apply to columns i,i+1 of U */
    if (U)
      for (j = 0; j < U->rows; j++) {
	x = ua[j + i * ldu];
	y = ua[j + (i + 1) * ldu];
	ua[j + i * ldu] = c * x + s * y;
	ua[j + (i + 1) * ldu] = -s * x + c * y;
      }

    /* Determine Givens rotation to eliminate carry... */
    givens(a[i], carry, &c, &s);

    /* ... and apply it to columns i,i+1 of row i, ... */
    a[i] = c * a[i] + s * carry;

    /* ... and columns i,i+1 of row i+1 */
    xr = b[i];
    yr = a[i + 1];
    b[i] = c * xr + s * yr;
    a[i + 1] = -s * xr + c * yr;

    /* Apply to rows i,i+1 of Vt */
    if (Vt)
      for (j = 0; j < Vt->cols; j++) {
	x = va[i + j * ldv];
	y = va[(i + 1) + j * ldv];
	va[i + j * ldv] = c * x + s * y;
	va[(i + 1) + j * ldv] = -s * x + c * y;
      }
  }
}

/* ------------------------------------------------------------
 * Eliminate subdiagonal element for zero diagonal element
 * ------------------------------------------------------------ */

static void
elimsubdiag(ptridiag T, pamatrix Vt)
{
  double*     a = T->d;
  double*     b = T->l;
  int      n = T->size;
  pfield    va;
  int      ldv;
  double      c, s, carry;
  field     x, y;
  double      yr;
  int      i, j;

  assert(Vt == NULL || Vt->rows == n);

  va = (Vt ? Vt->a : NULL);
  ldv = (Vt ? Vt->ld : 0);

  /* Quick exit if trivial matrix */
  if (n < 2)
    return;

  /* Determine Givens rotation to eliminate subdiagonal... */
  givens(a[1], b[0], &c, &s);

  /* ... and apply it column-wise to the second row */
  a[1] = c * a[1] + s * b[0];
  b[0] = 0.0;

  /* ... and its adjoint row-wise to Vt */
  if (Vt) {
    for (j = 0; j < Vt->cols; j++) {
      x = va[j * ldv];
      y = va[1 + j * ldv];
      va[j * ldv] = c * x - s * y;
      va[1 + j * ldv] = s * x + c * y;
    }
  }

  /* Chase away coefficient at (2,0) */
  for (i = 1; i < n - 1; i++) {
    /* Apply Givens rotation to row i+1 */
    yr = b[i];
    b[i] = c * yr;
    carry = -s * yr;

    /* Determine Givens rotation to eliminate (i+1,0)... */
    givens(a[i + 1], carry, &c, &s);

    /* ... and apply it column-wise to the (i+1)th row */
    a[i + 1] = c * a[i + 1] + s * carry;

    /* ... and its adjoint row-wise to Vt */
    if (Vt) {
      for (j = 0; j < Vt->cols; j++) {
	x = va[j * ldv];
	y = va[(i + 1) + j * ldv];
	va[j * ldv] = c * x - s * y;
	va[(i + 1) + j * ldv] = s * x + c * y;
      }
    }
  }
}

/* ------------------------------------------------------------
 * Singular value decomposition of a sub-bidiagonal matrix
 * ------------------------------------------------------------ */

int
sb_mulsvd_tridiag(ptridiag T, pamatrix U, pamatrix Vt, int maxiter)
{
  tridiag   tmp;
  amatrix   tmp2, tmp3;
  evpsortdata esd;
  ptridiag  Tsub;
  pamatrix  Usub, Vsub;
  double*     a = T->d;
  double*     b = T->l;
  int      n = T->size;
  pfield    ua;
  int      ldu;
  double      a0, b0, a1, lambda, lambda1, d, m, db, Tnorm;
  int      i, j, off, size, iter;

  ua = (U ? U->a : NULL);
  ldu = (U ? U->ld : 0);

  /* Simple case: 1-by-1 matrix */
  if (n < 2) {
    /* Ensure positive sign */
    if (a[0] < 0.0) {
      a[0] = -a[0];

      if (U)
	for (j = 0; j < U->rows; j++)
	  ua[j] = -ua[j];
    }

    return 0;
  }

  /* Compute maximum norm of T */
  Tnorm = REAL_ABS(a[0]);
  for (i = 1; i < n; i++) {
    d = REAL_ABS(b[i - 1]) + REAL_ABS(a[i]);
    if (d > Tnorm)
      Tnorm = d;
  }

  /* Find first non-zero subdiagonal coefficient */
  off = 0;
  while (off < n - 1 &&
	 REAL_ABS(b[off]) <=
	 H2_QR_EPS * (REAL_ABS(a[off]) + REAL_ABS(a[off + 1]))) {
    b[off] = 0.0;
    off++;
  }

  iter = 0;
  while (off < n - 1 && (maxiter == 0 || iter < maxiter)) {
    iter++;

    /* Find dimension of unreduced submatrix */
    size = 2;
    while (size < n - off
	   &&
	   REAL_ABS(b[off + size - 1]) > H2_QR_EPS
	   * (REAL_ABS(a[off + size - 1]) + REAL_ABS(a[off + size])))
      size++;

    /* Check for zero diagonal entries */
    i = 0;
    while (i < size && REAL_ABS(a[off + i]) > H2_QR_EPS * Tnorm)
      i++;

    if (i + 1 < size) {
      /* Eliminate subdiagonal entry */
      a[off + i] = 0.0;
      Tsub = init_sub_tridiag(&tmp, T, size - i, off + i);
      Vsub =
	(Vt ? init_sub_amatrix(&tmp3, Vt, size - i, off + i, Vt->cols, 0) :
	 NULL);

      elimsubdiag(Tsub, Vsub);

      if (Vsub)
	uninit_amatrix(Vsub);
      uninit_tridiag(Tsub);
    }
    else {
      /* Set up submatrices */
      Tsub = init_sub_tridiag(&tmp, T, size, off);
      Usub = (U ? init_sub_amatrix(&tmp2, U, U->rows, 0, size, off) : NULL);
      Vsub =
	(Vt ? init_sub_amatrix(&tmp3, Vt, size, off, Vt->cols, 0) : NULL);

      /* Determine Wilkinson shift for T T^* */
      a0 = REAL_SQR(a[off + size - 2]);
      if (off + size > 2)
	a0 += REAL_SQR(b[off + size - 3]);
      b0 = b[off + size - 2] * a[off + size - 2];
      a1 = REAL_SQR(b[off + size - 2]) + REAL_SQR(a[off + size - 1]);

      /* Find the eigenvalue with larger absolute value */
      m = (a0 + a1) * 0.5;
      d = (a0 - a1) * 0.5;
      db = (REAL_ABS(d) > REAL_ABS(b0) ?
	    REAL_ABS(d) * REAL_SQRT(1.0 + REAL_SQR(b0 / d)) :
	    REAL_ABS(b0) * REAL_SQRT(1.0 + REAL_SQR(d / b0)));
      lambda1 = (m > 0 ? m + db : m - db);

      /* Find the eigenvalue closest to a1 */
      lambda = (m * d < 0.0 ? lambda1 : (a0 * a1 - REAL_SQR(b0)) / lambda1);

      /* Perform SVD step on submatrices */
      svdstep_tridiag(Tsub, lambda, Usub, Vsub);

      /* Clean up submatrices */
      if (Vsub)
	uninit_amatrix(Vsub);
      if (Usub)
	uninit_amatrix(Usub);
      uninit_tridiag(Tsub);
    }

    /* Find non-zero subdiagonal coefficient */
    while (off < n - 1 &&
	   REAL_ABS(b[off]) <=
	   H2_QR_EPS * (REAL_ABS(a[off]) + REAL_ABS(a[off + 1]))) {
      b[off] = 0.0;
      off++;
    }
  }

  if (maxiter > 0 && iter >= maxiter) {
    (void) printf("  ### SVD warning: did not converge after %u steps\n",
		  iter);
  }

  /* Ensure positive signs */
  for (i = 0; i < n; i++)
    if (a[i] < 0.0) {
      a[i] = -a[i];

      if (U)
	for (j = 0; j < U->rows; j++)
	  ua[j + i * ldu] = -ua[j + i * ldu];
    }

  /* Sort singular values */
  esd.T = T;
  esd.U = U;
  esd.Vt = Vt;
  heapsort(n, evp_geq, evp_swap, &esd);

  return iter;
}

#ifdef USE_BLAS
int
mulsvd_tridiag(ptridiag T, pamatrix U, pamatrix Vt)
{
  int      lwork;
  double     *work;
  int       info;

  if (T->size < 1)
    return 0;

  lwork = 4 * T->size;
  work = allocreal(lwork);

  info = h2_bdsqr(LAPACK_COL_MAJOR, 'L', T->size, (Vt ? Vt->cols : i_zero), (U ? U->rows : i_zero), i_zero, T->d, T->l, (MKL_Complex16 *)(Vt ? Vt->a : 0), (Vt ? Vt->ld : i_one), (MKL_Complex16 *)(U ? U->a : 0), (U ? U->ld : i_one), 0, i_one);

  freemem(work);

  return (info != 0);
}
#else
int
mulsvd_tridiag(ptridiag T, pamatrix U, pamatrix Vt)
{
  int      iter, maxiter;

  if (T->size < 1)
    return 0;

  maxiter = 32 * T->size;

  iter = sb_mulsvd_tridiag(T, U, Vt, maxiter);

  return !(iter < maxiter);
}
#endif

int
svd_tridiag(ptridiag T, pamatrix U, pamatrix Vt)
{
  if (U)
    identity_amatrix(U);

  if (Vt)
    identity_amatrix(Vt);

  return mulsvd_tridiag(T, U, Vt);
}

/* ------------------------------------------------------------
 * Bidiagonalize a matrix
 * ------------------------------------------------------------ */

void
sb_bidiagonalize_amatrix(pamatrix A, ptridiag T, pamatrix U, pamatrix Vt)
{
  pfield    a, ua, va;
  double    *d, *l, *tau;
  field     alpha, gamma, diag, nsign;
  double      norm, norm2, beta;
  int      rows, cols, lda, ldu, ldv;
  int      size;
  int      i, j, k;

  rows = A->rows;
  cols = A->cols;
  size = UINT_MIN(rows, cols);

  assert(T->size == size);
  assert(U == NULL || U->rows >= A->rows);
  assert(U == NULL || U->cols >= min<int>(A->rows, A->cols));
  assert(Vt == NULL || Vt->cols >= A->cols);
  assert(Vt == NULL || Vt->rows >= min<int>(A->rows, A->cols));

  a = A->a;
  lda = A->ld;

  d = T->d;
  l = T->l;
  tau = d;

  ua = (U ? U->a : NULL);
  ldu = (U ? U->ld : 0);

  va = (Vt ? Vt->a : NULL);
  ldv = (Vt ? Vt->ld : 0);

  if (U)
    identity_amatrix(U);
  if (Vt)
    identity_amatrix(Vt);

  if (rows > cols) {		/* Make A upper triangular */
    for (k = 0; k < cols; k++) {
      norm2 = ABSSQR(a[k + k * lda]);
      for (i = k + 1; i < rows; i++)
	norm2 += ABSSQR(a[i + k * lda]);
      norm = REAL_SQRT(norm2);

      if (norm2 <= H2_ALMOST_ZERO) {
	tau[k] = 0.0;
	a[k + k * lda] = 0.0;
      }
      else {
	/* Determine Householder reflection vector v */
	diag = a[k + k * lda];
	alpha = -SIGN1(diag) * norm;
	gamma = diag - alpha;

	/* Compute norm of v */
	beta = 1.0 / (norm2 + ABS(diag) * norm);

	/* Normalize reflection vector */
	for (i = k + 1; i < rows; i++)
	  a[i + k * lda] /= gamma;
	beta *= ABSSQR(gamma);

	/* Store scaling factor */
	tau[k] = beta;

	/* Store diagonal element */
	a[k + k * lda] = alpha;

	/* Update columns k+1,...,cols */
	for (j = k + 1; j < cols; j++) {
	  gamma = a[k + j * lda];
	  for (i = k + 1; i < rows; i++)
	    gamma += CONJ(a[i + k * lda]) * a[i + j * lda];

	  gamma *= beta;

	  a[k + j * lda] -= gamma;
	  for (i = k + 1; i < rows; i++)
	    a[i + j * lda] -= gamma * a[i + k * lda];
	}
      }
    }

    /* Apply reflections in reversed order to U */
    if (U)
      for (k = cols; k-- > 0;) {
	beta = tau[k];
	if (beta != 0.0) {
	  for (j = 0; j < U->cols; j++) {
	    gamma = ua[k + j * ldu];
	    for (i = k + 1; i < rows; i++)
	      gamma += CONJ(a[i + k * lda]) * ua[i + j * ldu];

	    gamma *= beta;

	    ua[k + j * ldu] -= gamma;
	    for (i = k + 1; i < rows; i++)
	      ua[i + j * ldu] -= gamma * a[i + k * lda];
	  }
	}
      }

    /* Clear subdiagonal entries */
    for (k = 0; k < cols; k++)
      for (i = k + 1; i < cols; i++)
	a[i + k * lda] = 0.0;

    /* Now A is quadratic */
    rows = cols;
  }
  else if (cols > rows) {	/* Make A lower triangular */
    for (k = 0; k < rows; k++) {
      norm2 = ABSSQR(a[k + k * lda]);
      for (i = k + 1; i < cols; i++)
	norm2 += ABSSQR(a[k + i * lda]);
      norm = REAL_SQRT(norm2);

      if (norm2 <= H2_ALMOST_ZERO) {
	tau[k] = 0.0;
	a[k + k * lda] = 0.0;
      }
      else {
	/* Determine Householder reflection vector v */
	diag = a[k + k * lda];
	alpha = -SIGN1(diag) * norm;
	gamma = diag - alpha;

	/* Compute norm of v */
	beta = 1.0 / (norm2 + ABS(diag) * norm);

	/* Normalize reflection vector */
	for (i = k + 1; i < cols; i++)
	  a[k + i * lda] /= gamma;
	beta *= ABSSQR(gamma);

	/* Store scaling factor */
	tau[k] = beta;

	/* Store diagonal element */
	a[k + k * lda] = alpha;

	/* Update rows k+1,...,rows */
	for (j = k + 1; j < rows; j++) {
	  gamma = a[j + k * lda];
	  for (i = k + 1; i < cols; i++)
	    gamma += CONJ(a[k + i * lda]) * a[j + i * lda];

	  gamma *= beta;

	  a[j + k * lda] -= gamma;
	  for (i = k + 1; i < cols; i++)
	    a[j + i * lda] -= gamma * a[k + i * lda];
	}
      }
    }

    /* Apply reflections in reversed order to Vt */
    if (Vt)
      for (k = rows; k-- > 0;) {
	beta = tau[k];
	if (beta != 0.0) {
	  for (j = 0; j < Vt->rows; j++) {
	    gamma = va[j + k * ldv];
	    for (i = k + 1; i < cols; i++)
	      gamma += CONJ(a[k + i * lda]) * va[j + i * ldv];

	    gamma *= beta;

	    va[j + k * ldv] -= gamma;
	    for (i = k + 1; i < cols; i++)
	      va[j + i * ldv] -= gamma * a[k + i * lda];
	  }
	}
      }

    /* Clear superdiagonal entries */
    for (k = 0; k < rows; k++)
      for (i = k + 1; i < rows; i++)
	a[k + i * lda] = 0.0;

    /* Now A is quadratic */
    cols = rows;
  }

  assert(size == rows);
  assert(size == cols);

  /* Golub-Kahan bidiagonalization */
  for (k = 0; k < size; k++) {
    /* Eliminate (k,k+1) to (k,cols) by column reflections */
    norm2 = ABSSQR(a[k + k * lda]);
    for (i = k + 1; i < cols; i++)
      norm2 += ABSSQR(a[k + i * lda]);
    norm = REAL_SQRT(norm2);

    if (norm2 <= H2_ALMOST_ZERO)
      d[k] = 0.0;
    else {
      /* Determine Householder reflection vector v */
      diag = a[k + k * lda];
      nsign = -SIGN1(diag);
      alpha = nsign * norm;
      gamma = diag - alpha;

      /* Compute 2 / |v|^2 */
      beta = 1.0 / (norm2 + ABS(diag) * norm);

      /* Store diagonal element */
      d[k] = norm;

      /* Normalize reflection vector */
      for (i = k + 1; i < cols; i++)
	a[k + i * lda] /= gamma;
      beta *= ABSSQR(gamma);

      /* Update rows k+1,...,rows */
      for (j = k + 1; j < rows; j++) {
	gamma = a[j + k * lda];
	for (i = k + 1; i < cols; i++)
	  gamma += CONJ(a[k + i * lda]) * a[j + i * lda];

	gamma *= beta;

	a[j + k * lda] -= gamma;
	for (i = k + 1; i < cols; i++)
	  a[j + i * lda] -= gamma * a[k + i * lda];

	a[j + k * lda] *= CONJ(nsign);
      }

      /* Update columns of Vt */
      if (Vt)
	for (j = 0; j < Vt->cols; j++) {
	  gamma = va[k + j * ldv];
	  for (i = k + 1; i < cols; i++)
	    gamma += a[k + i * lda] * va[i + j * ldv];

	  gamma *= beta;

	  va[k + j * ldv] -= gamma;
	  for (i = k + 1; i < cols; i++)
	    va[i + j * ldv] -= gamma * CONJ(a[k + i * lda]);

	  va[k + j * ldv] *= nsign;
	}
    }

    /* Eliminate (k+2,k) to (rows,k) by row reflections */
    if (k + 1 < rows) {
      norm2 = ABSSQR(a[(k + 1) + k * lda]);
      for (i = k + 2; i < rows; i++)
	norm2 += ABSSQR(a[i + k * lda]);
      norm = REAL_SQRT(norm2);

      if (norm2 <= H2_ALMOST_ZERO)
	l[k] = 0.0;
      else {
	/* Determine Householder reflection vector v */
	diag = a[(k + 1) + k * lda];
	nsign = -SIGN1(diag);
	alpha = nsign * norm;
	gamma = diag - alpha;

	/* Compute 2 / |v|^2 */
	beta = 1.0 / (norm2 + ABS(diag) * norm);

	/* Store subdiagonal element */
	l[k] = norm;

	/* Normalize reflection vector */
	for (i = k + 2; i < rows; i++)
	  a[i + k * lda] /= gamma;
	beta *= ABSSQR(gamma);

	/* Update columns k+1,...,cols */
	for (j = k + 1; j < cols; j++) {
	  gamma = a[(k + 1) + j * lda];
	  for (i = k + 2; i < rows; i++)
	    gamma += CONJ(a[i + k * lda]) * a[i + j * lda];

	  gamma *= beta;

	  a[(k + 1) + j * lda] -= gamma;
	  for (i = k + 2; i < rows; i++)
	    a[i + j * lda] -= gamma * a[i + k * lda];

	  a[(k + 1) + j * lda] *= CONJ(nsign);
	}

	/* Update rows of U */
	if (U)
	  for (j = 0; j < U->rows; j++) {
	    gamma = ua[j + (k + 1) * ldu];
	    for (i = k + 2; i < rows; i++)
	      gamma += a[i + k * lda] * ua[j + i * ldu];

	    gamma *= beta;

	    ua[j + (k + 1) * ldu] -= gamma;
	    for (i = k + 2; i < rows; i++)
	      ua[j + i * ldu] -= gamma * CONJ(a[i + k * lda]);

	    ua[j + (k + 1) * ldu] *= nsign;
	  }
      }
    }
  }
}

#if 0
void
bidiagonalize_amatrix(pamatrix A, ptridiag T, pamatrix U, pamatrix Vt)
{
  pfield    a, ua, va, work, tau;
  double     *d, *l;
  field     alpha, beta, nsign;
  int      rows, cols, lda, ldu, ldv;
  int      rows1, cols1;
  int      size, lwork;
  int      i, j, k;

  rows = A->rows;
  cols = A->cols;
  size = UINT_MIN(rows, cols);

  assert(T->size == size);
  assert(U == NULL || U->rows >= A->rows);
  assert(U == NULL || U->cols >= UINT_MIN(A->rows, A->cols));
  assert(Vt == NULL || Vt->cols >= A->cols);
  assert(Vt == NULL || Vt->rows >= UINT_MIN(A->rows, A->cols));

  lwork = UINT_MAX(A->rows, A->cols);
  work = allocfield(lwork);

  a = A->a;
  lda = A->ld;

  d = T->d;
  l = T->l;

  ua = (U ? U->a : NULL);
  ldu = (U ? U->ld : 0);

  va = (Vt ? Vt->a : NULL);
  ldv = (Vt ? Vt->ld : 0);

  if (U)
    identity_amatrix(U);
  if (Vt)
    identity_amatrix(Vt);

  if (rows > cols) {		/* Make A upper triangular */
    tau = allocfield(cols);

    for (k = 0; k < cols; k++) {
      rows1 = rows - k;

      /* Determine Householder reflection vector v */
      alpha = a[k + k * lda];
      h2_larfg(rows1, (MKL_Complex16 *)&alpha, (MKL_Complex16 *)(a + (k + 1) + k * lda), i_one, (MKL_Complex16 *)&beta);

      /* Store scaling factor */
      tau[k] = beta;
      beta = CONJ(beta);

      /* dlarf/zlarf expects the full Householder vector */
      a[k + k * lda] = 1.0;

      /* Update columns k+1,...,cols */
      cols1 = cols - (k + 1);
      h2_larf(_h2_left, &rows1, &cols1, a + k + k * lda, &i_one, &beta,
	      a + k + (k + 1) * lda, &lda, work);

      /* Store diagonal element */
      a[k + k * lda] = alpha;
    }

    /* Apply reflections in reversed order to U */
    if (U)
      for (k = cols; k-- > 0;) {
	beta = tau[k];

	if (beta != 0.0) {
	  alpha = a[k + k * lda];
	  a[k + k * lda] = 1.0;

	  rows1 = rows - k;
	  cols1 = U->cols;
	  h2_larf(_h2_left, &rows1, &cols1, a + k + k * lda, &i_one, &beta,
		  ua + k, &ldu, work);

	  a[k + k * lda] = alpha;
	}
      }

    /* Clear subdiagonal entries */
    for (k = 0; k < cols; k++)
      for (i = k + 1; i < cols; i++)
	a[i + k * lda] = 0.0;

    /* Clean up */
    freemem(tau);

    /* Now A is quadratic */
    rows = cols;
  }
  else if (cols > rows) {	/* Make A lower triangular */
    tau = allocfield(rows);

    for (k = 0; k < rows; k++) {
      cols1 = cols - k;

      /* Determine Householder reflection vector v */
      alpha = CONJ(a[k + k * lda]);
#ifdef USE_COMPLEX
      for (i = k + 1; i < cols; i++)
	a[k + i * lda] = CONJ(a[k + i * lda]);
#endif
      h2_larfg(&cols1, &alpha, a + k + (k + 1) * lda, &lda, &beta);

      /* Store scaling factor */
      tau[k] = beta;

      /* dlarf/zlarf expects the full Householder vector */
      a[k + k * lda] = 1.0;

      /* Update rows k+1,...,rows */
      rows1 = rows - (k + 1);
      h2_larf(_h2_right, &rows1, &cols1, a + k + k * lda, &lda, &beta,
	      a + (k + 1) + k * lda, &lda, work);

      /* Store diagonal element */
      a[k + k * lda] = CONJ(alpha);
    }

    /* Apply reflections in reversed order to Vt */
    if (Vt)
      for (k = rows; k-- > 0;) {
	beta = CONJ(tau[k]);

	if (beta != 0.0) {
	  alpha = a[k + k * lda];
	  a[k + k * lda] = 1.0;

	  rows1 = Vt->rows;
	  cols1 = cols - k;
	  h2_larf(_h2_right, &rows1, &cols1, a + k + k * lda, &lda, &beta,
		  va + k * ldv, &ldv, work);

	  a[k + k * lda] = alpha;
	}
      }

    /* Clear superdiagonal entries */
    for (k = 0; k < rows; k++)
      for (i = k + 1; i < rows; i++)
	a[k + i * lda] = 0.0;

    /* Clean up */
    freemem(tau);

    /* Now A is quadratic */
    cols = rows;
  }

  assert(size == rows);
  assert(size == cols);

  /* Golub-Kahan bidiagonalization */
  for (k = 0; k < size; k++) {
    /* Eliminate (k,k+1) to (k,cols) by column reflections */
    cols1 = size - k;

    /* Determine Householder reflection vector v */
    alpha = CONJ(a[k + k * lda]);

    for (i = k + 1; i < size; i++)
      a[k + i * lda] = CONJ(a[k + i * lda]);

    h2_larfg(&cols1, &alpha, a + k + (k + 1) * lda, &lda, &beta);

    /* Store diagonal element */
    nsign = SIGN1(CONJ(alpha));
    d[k] = ABS(alpha);

    /* dlarf/zlarf expects the full Householder vector */
    a[k + k * lda] = 1.0;

    /* Update rows k+1,...,rows */
    rows1 = size - (k + 1);
    h2_larf(_h2_right, &rows1, &cols1, a + k + k * lda, &lda, &beta,
	    a + (k + 1) + k * lda, &lda, work);

    /* Scale to preserve sign */
    for (j = k + 1; j < rows; j++)
      a[j + k * lda] *= nsign;

    /* Update columns of Vt */
    if (Vt) {
      beta = CONJ(beta);
      rows1 = size - k;
      cols1 = Vt->cols;
      h2_larf(_h2_left, &rows1, &cols1, a + k + k * lda, &lda, &beta, va + k,
	      &ldv, work);

      for (j = 0; j < Vt->cols; j++)
	va[k + j * ldv] *= CONJ(nsign);
    }

    /* Eliminate (k+2,k) to (rows,k) by row reflections */
    if (k + 1 < size) {
      rows1 = size - (k + 1);

      /* Determine Householder reflection vector v */
      alpha = a[(k + 1) + k * lda];
      h2_larfg(&rows1, &alpha, a + (k + 2) + k * lda, &i_one, &beta);

      /* Store subdiagonal element */
      nsign = SIGN1(alpha);
      l[k] = ABS(alpha);

      /* dlarf/zlarf expects the full Householder vector */
      a[(k + 1) + k * lda] = 1.0;

      /* Update columns k+1,...,cols */
      beta = CONJ(beta);
      cols1 = size - k - 1;
      h2_larf(_h2_left, &rows1, &cols1, a + (k + 1) + k * lda, &i_one, &beta,
	      a + (k + 1) + (k + 1) * lda, &lda, work);

      for (j = k + 1; j < cols; j++)
	a[(k + 1) + j * lda] *= CONJ(nsign);

      /* Update rows of U */
      if (U) {
	beta = CONJ(beta);
	rows1 = U->rows;
	cols1 = rows - (k + 1);
	h2_larf(_h2_right, &rows1, &cols1, a + (k + 1) + k * lda, &i_one,
		&beta, ua + (k + 1) * ldu, &ldu, work);

	for (j = 0; j < U->rows; j++)
	  ua[j + (k + 1) * ldu] *= nsign;
      }
    }
  }

  /* Clean up */
  freemem(work);
}
#else
void
bidiagonalize_amatrix(pamatrix A, ptridiag T, pamatrix U, pamatrix Vt)
{
  sb_bidiagonalize_amatrix(A, T, U, Vt);
}
#endif

void
bidiagonalize_verified_amatrix(pamatrix A, ptridiag T, pamatrix U,
			       pamatrix Vt)
{
  pamatrix  Acopy, Utmp, Vttmp;
  amatrix   tmp1, tmp2, tmp3;
  int      k;
  double      norm, error;

  k = UINT_MIN(A->rows, A->cols);

  Acopy = init_amatrix(&tmp1, A->rows, A->cols);
  copy_amatrix(false, A, Acopy);

  Utmp = init_amatrix(&tmp2, A->rows, k);
  Vttmp = init_amatrix(&tmp3, k, A->cols);

  bidiagonalize_amatrix(A, T, Utmp, Vttmp);

  if (U)
    copy_amatrix(false, Utmp, U);

  if (Vt)
    copy_amatrix(false, Vttmp, Vt);

  if (T->size > 0) {
    norm = normfrob_amatrix(Acopy);
    lowereval_tridiag_amatrix(1.0, true, T, true, Utmp);
    addmul_amatrix(-1.0, false, Utmp, false, Vttmp, Acopy);
    error = normfrob_amatrix(Acopy);
    if (error > 1e-12 * norm) {
      printf("  bidiag: %.4g\n", error);
      abort();
    }
  }

  uninit_amatrix(Vttmp);
  uninit_amatrix(Utmp);
  uninit_amatrix(Acopy);
}

/* ------------------------------------------------------------
 * Singular value decomposition of an arbitrary matrix
 * ------------------------------------------------------------ */

int
sb_svd_amatrix(pamatrix A, prealavector sigma, pamatrix U, pamatrix Vt,
	       int maxiter)
{
#ifdef RUNTIME_CHECK_EIGENSOLVERS
  tridiag   tmp1;
  amatrix   tmp2, tmp3, tmp4, tmp5;
  ptridiag  T;
  pamatrix  Acopy, Utmp, Vttmp, UT;
  double      norm, error;
  int      k;
  int      iter;
  int      i;

  k = UINT_MIN(A->rows, A->cols);

  if (k < 1)			/* Quick exit */
    return 0;

  Acopy = init_amatrix(&tmp2, A->rows, A->cols);
  Utmp = init_amatrix(&tmp3, A->rows, k);
  Vttmp = init_amatrix(&tmp4, k, A->cols);
  UT = init_amatrix(&tmp5, A->rows, k);

  copy_amatrix(false, A, Acopy);

  /* Set up auxiliary matrix */
  T = init_tridiag(&tmp1, k);

  /* Bidiagonalize A */
  sb_bidiagonalize_amatrix(A, T, Utmp, Vttmp);

  if (U)
    copy_amatrix(false, Utmp, U);

  if (Vt)
    copy_amatrix(false, Vttmp, Vt);

  /* Check bidiagonalization */
  norm = normfrob_amatrix(Acopy);

  copy_amatrix(false, Utmp, UT);
  lowereval_tridiag_amatrix(1.0, true, T, true, UT);

  addmul_amatrix(-1.0, false, UT, false, Vttmp, Acopy);

  error = normfrob_amatrix(Acopy);

  if (error > H2_CHECK_TOLERANCE * norm)
    (void) printf("  ### Poor bidiagonalization accuracy, %g\n",
		  error / norm);

  /* Compute SVD of bidiagonal matrix */
  iter = sb_mulsvd_tridiag(T, U, Vt, maxiter);

  if (maxiter > 0 && iter >= maxiter)
    (void) printf("  ## SVD iteration did not converge after %u steps\n",
		  iter);

  /* Copy singular values */
  for (i = 0; i < k; i++)
    sigma->v[i] = T->d[i];

  /* Clean up */
  uninit_tridiag(T);
  uninit_amatrix(UT);
  uninit_amatrix(Vttmp);
  uninit_amatrix(Utmp);
  uninit_amatrix(Acopy);
#else
  tridiag   tmp;
  ptridiag  T;
  int      k;
  int      iter;
  int      i;

  cout << "in sb svd" << endl;

  k = UINT_MIN(A->rows, A->cols);

  if (k < 1)			/* Quick exit */
    return 0;

  /* Set up auxiliary matrix */
  T = init_tridiag(&tmp, k);

  /* Bidiagonalize A */
  sb_bidiagonalize_amatrix(A, T, U, Vt);

  /* Compute SVD of bidiagonal matrix */
  iter = sb_mulsvd_tridiag(T, U, Vt, maxiter);

  /* Copy singular values */
  for (i = 0; i < k; i++)
    sigma->v[i] = T->d[i];

  /* Clean up */
  uninit_tridiag(T);
#endif

  return iter;
}

#ifdef USE_BLAS
#if defined(THREADSAFE_LAPACK) || !defined(USE_OPENMP)
int
svd_amatrix(pamatrix A, prealavector sigma, pamatrix U, pamatrix Vt)
{
  pfield    work;
  int  lwork;
  double*     rwork;

  int       info = 0;

  assert(sigma->dim >= min<int>(A->rows, A->cols));

  if (A->rows > 0 && A->cols > 0) {
    lwork = 10 * UINT_MAX(A->rows, A->cols);
    work = allocfield(lwork);

    rwork = allocreal(5 * UINT_MIN(A->rows, A->cols));

    h2_gesvd(h2_col_major, (U ? 'S' : 'N'),
	     (Vt ? 'S' : 'N'),
	     A->rows, A->cols,
	     (MKL_Complex16 *)(A->a), A->ld,
	     sigma->v,
	     (MKL_Complex16 *)(U ? U->a : NULL), (U ? U->ld : i_one),
	     (MKL_Complex16 *)(Vt ? Vt->a : NULL), (Vt ? Vt->ld : i_one),
	     rwork);

    freemem(rwork);
    freemem(work);
  }

  return (info != 0);
}
#else
int
svd_amatrix(pamatrix A, prealavector sigma, pamatrix U, pamatrix Vt)
{
  tridiag   tmp;
  ptridiag  T;
  int      size;
  int      info;
  int      i;

  size = UINT_MIN(A->rows, A->cols);

  /* Set up auxiliary matrix */
  T = init_tridiag(&tmp, size);

  /* Bidiagonalize A */
  bidiagonalize_amatrix(A, T, U, Vt);

  /* Compute SVD of bidiagonal matrix */
  info = mulsvd_tridiag(T, U, Vt);

  /* Copy singular values */
  for (i = 0; i < size; i++)
    sigma->v[i] = T->d[i];

  /* Clean up */
  uninit_tridiag(T);

  return info;
}
#endif
#else
int
svd_amatrix(pamatrix A, prealavector sigma, pamatrix U, pamatrix Vt)
{
#ifdef RUNTIME_CHECK_EIGENSOLVERS
  int      iter, maxiter;
  pamatrix  Acopy, Utmp, Vttmp;
  amatrix   tmp1, tmp2, tmp3;
  double      norm, error;
  int      k;

  k = UINT_MIN(A->rows, A->cols);


  if (k < 1)			/* Quick exit */
    return 0;

  Acopy = init_amatrix(&tmp1, A->rows, A->cols);
  copy_amatrix(false, A, Acopy);

  Utmp = init_amatrix(&tmp2, A->rows, k);
  Vttmp = init_amatrix(&tmp3, k, A->cols);

  maxiter = 32 * k;

  iter = sb_svd_amatrix(A, sigma, Utmp, Vttmp, maxiter);

  if (U)
    copy_amatrix(false, Utmp, U);
  if (Vt)
    copy_amatrix(false, Vttmp, Vt);

  norm = normfrob_amatrix(Acopy);

  diageval_realavector_amatrix(1.0, true, sigma, true, Utmp);
  addmul_amatrix(-1.0, false, Utmp, false, Vttmp, Acopy);

  error = normfrob_amatrix(Acopy);

  if (error > H2_CHECK_TOLERANCE * norm)
    (void) printf("  ### Poor SVD accuracy, %g\n", error / norm);

  uninit_amatrix(Vttmp);
  uninit_amatrix(Utmp);
  uninit_amatrix(Acopy);
#else

  int      iter, maxiter;
  int      k;

  k = UINT_MIN(A->rows, A->cols);

  if (k < 1)			/* Quick exit */
    return 0;

  maxiter = 32 * k;

  iter = sb_svd_amatrix(A, sigma, U, Vt, maxiter);
#endif

  return !(iter < maxiter);
}
#endif
