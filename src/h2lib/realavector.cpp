
#include <stdio.h>
#include "realavector.h"

static int active_realavector = 0;

prealavector
init_realavector(prealavector v, int dim)
{
  assert(v != NULL);

  v->v = (dim > 0 ? allocreal(dim) : NULL);
  v->dim = dim;
  v->owner = NULL;

#ifdef USE_OPENMP
#pragma omp atomic
#endif
  active_realavector++;

  return v;
}

prealavector
init_sub_realavector(prealavector v, prealavector src, int dim, int off)
{
  assert(v != NULL);
  assert(src != NULL);
  assert(off + dim <= src->dim);

  v->v = src->v + off;
  v->dim = dim;
  v->owner = src;

#ifdef USE_OPENMP
#pragma omp atomic
#endif
  active_realavector++;

  return v;
}

prealavector
init_pointer_realavector(prealavector v, double* src, int dim)
{
  assert(v != NULL);
  assert(dim == 0 || src != NULL);

  v->v = src;
  v->dim = dim;
  v->owner = src;

#ifdef USE_OPENMP
#pragma omp atomic
#endif
  active_realavector++;

  return v;
}

void
uninit_realavector(prealavector v)
{
  if (!v->owner)
    freemem(v->v);

  assert(active_realavector > 0);

#ifdef USE_OPENMP
#pragma omp atomic
#endif
  active_realavector--;
}

prealavector
new_realavector(int dim)
{
  prealavector v;

  v = (prealavector) allocmem(sizeof(realavector));

  init_realavector(v, dim);

  return v;
}

prealavector
new_sub_realavector(prealavector src, int dim, int off)
{
  prealavector v;

  v = (prealavector) allocmem(sizeof(realavector));

  init_sub_realavector(v, src, dim, off);

  return v;
}

prealavector
new_pointer_realavector(double* src, int dim)
{
  prealavector v;

  v = (prealavector) allocmem(sizeof(realavector));

  init_pointer_realavector(v, src, dim);

  return v;
}

void
del_realavector(prealavector v)
{
  uninit_realavector(v);
  freemem(v);
}

void
resize_realavector(prealavector v, int dim)
{
  assert(v->owner == NULL);

  if (dim != v->dim) {
    freemem(v->v);
    v->v = allocreal(dim);
    v->dim = dim;
  }
}

void
shrink_realavector(prealavector v, int dim)
{
  assert(dim <= v->dim);

  v->dim = dim;
}

/* ------------------------------------------------------------
 Statistics
 ------------------------------------------------------------ */

int
getactives_realavector()
{
  return active_realavector;
}

size_t
getsize_realavector(pcrealavector v)
{
  size_t    sz;

  sz = sizeof(realavector);
  if (v->owner == NULL)
    sz += (size_t) sizeof(double) * v->dim;

  return sz;
}

size_t
getsize_heap_realavector(pcrealavector v)
{
  size_t    sz;

  sz = 0;
  if (v->owner == NULL)
    sz += (size_t) sizeof(double) * v->dim;

  return sz;
}

/* ------------------------------------------------------------
 Simple utility functions
 ------------------------------------------------------------ */

void
clear_realavector(prealavector v)
{
  int      i;

  for (i = 0; i < v->dim; i++)
    v->v[i] = 0.0;
}

void
fill_realavector(prealavector v, double x)
{
  int      i;

  for (i = 0; i < v->dim; i++)
    v->v[i] = x;
}

void
random_realavector(prealavector v)
{
  int      i;

  for (i = 0; i < v->dim; i++)
    v->v[i] = REAL_RAND();
}

void
copy_realavector(pcrealavector v, prealavector w)
{
  int      i;

  assert(v->dim == w->dim);

  for (i = 0; i < v->dim; i++)
    w->v[i] = v->v[i];
}

void
copy_sub_realavector(pcrealavector v, prealavector w)
{
  int      i, n;

  n = UINT_MIN(v->dim, w->dim);

  for (i = 0; i < n; i++)
    w->v[i] = v->v[i];
}

void
print_realavector(pcrealavector v)
{
  int      dim = v->dim;
  int      i;

  (void) printf("realavector(%u)\n", dim);
  if (dim == 0)
    return;

  (void) printf("  (%.5e", v->v[0]);
  for (i = 1; i < dim; i++)
    (void) printf(" %.5e", v->v[i]);
  (void) printf(")\n");
}

/* ------------------------------------------------------------
 Very basic linear algebra
 ------------------------------------------------------------ */

#ifdef USE_BLAS
void
scale_realavector(double alpha, prealavector v)
{
    int      dim = v->dim / 2;
    h2_rscal(dim, alpha, (field *) v->v, i_one);
}
#else
void
scale_realavector(double alpha, prealavector v)
{
  int      i;

  for (i = 0; i < v->dim; i++)
    v->v[i] *= alpha;
}
#endif

#ifdef USE_BLAS
double norm2_realavector(pcrealavector v)
{
  int      dim = v->dim / 2;

  return h2_nrm2(dim, (field *) v->v, i_one);
}
#else
double
norm2_realavector(pcrealavector v)
{
  double      sum;
  int      i;

  sum = 0.0;
  for (i = 0; i < v->dim; i++)
    sum += REAL_SQR(v->v[i]);

  return REAL_SQRT(sum);
}
#endif

#ifdef USE_BLAS
double
dotprod_realavector(pcrealavector x, pcrealavector y)
{
  assert(x->dim == y->dim);

  return h2_rdot(x->dim, x->v, i_one, y->v, i_one);
}
#else
double
dotprod_realavector(pcrealavector x, pcrealavector y)
{
  double alpha;
  int      i;

  assert(x->dim == y->dim);

  alpha = 0.0;
  for (i = 0; i < x->dim; i++)
    alpha += x->v[i] * y->v[i];

  return alpha;
}
#endif

#ifdef USE_BLAS
void
add_realavector(double alpha, pcrealavector x, prealavector y)
{
    assert(x->dim == y->dim);
    h2_raxpy(x->dim, &alpha, x->v, i_one, y->v, i_one);
}
#else
void
add_realavector(double alpha, pcrealavector x, prealavector y)
{
  int      i;

  for (i = 0; i < x->dim; i++)
    y->v[i] += alpha * x->v[i];
}
#endif
