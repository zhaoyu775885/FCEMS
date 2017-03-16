#include "amatrix.h"
#include "block.h"
#include <math.h>
#include <assert.h>
#include "blas.h"

/// The Layout option should be left for special ctrl
/// however, the trans option must be specified clearly

pamatrix
init_amatrix(pamatrix a, int rows, int cols)
{
	assert(a != NULL);
	a->a = (rows>0 && cols>0 ? new field [rows*cols]() : NULL);
    a->ld = rows;// 2 Row Major now!
    a->rows = rows;
    a->cols = cols;
    a->owner = NULL;

	return a;
}

void
uninit_amatrix(pamatrix a)
{
    if (!a->owner) {
		delete [] a->a;
    }
}

void
print_amatrix(pamatrix a)
{
		for (int i=0;i<a->rows;++i) {
	for (int j=0;j<a->cols;++j) {
			cout << a->a[i+j*a->rows] << " ";
		}
		cout << endl;
	}
	cout << endl;
}

pamatrix
new_amatrix(int rows, int cols)
{
    pamatrix am = new amatrix;
    init_amatrix(am, rows, cols);
    return am;
}

void
del_amatrix(pamatrix a)
{
    uninit_amatrix(a);
    free(a);
}

pamatrix
new_zero_amatrix(int rows, int cols)
{
    pamatrix am = new amatrix;
    init_amatrix(am, rows, cols);
    for (int j=0;j<cols;++j)
        for (int i=0;i<rows;++i)
            am->a[i+j*rows] = 0;
    return am;
}

pamatrix
new_ones_amatrix(int rows, int cols)
{
    pamatrix am = new amatrix;
    init_amatrix(am, rows, cols);
    for (int j=0;j<cols;++j)
        for (int i=0;i<rows;++i)
            am->a[i+j*rows] = 1;
    return am;
}

void
addeval_amatrix_avector(field alpha, pamatrix am, char layout, char trans, pcavector x, field beta, pavector y)
{
    switch (trans) {
        case CblasNoTrans:      y->dim = am->rows; break;
        case CblasTrans:        y->dim = am->cols; break;
        case CblasConjTrans:    y->dim = am->cols; break;
        default:    cout << "Error! Clarify the Trans Type\n" << endl;
    }
    cblas_zgemv((CBLAS_ORDER)layout, (CBLAS_TRANSPOSE)trans, am->rows, am->cols, &alpha, am->a, am->ld, x->v, 1, &beta, y->v, 1);
}


void
init_eye_amatrix(pamatrix am, int n)
{
        for (int j=0;j<n;j++) {
    for (int i=0;i<n;i++) {
            if (i == j) {
                am->a[i+j*n] = 1;
            }
            else {
                am->a[i*n+j] = 0;
            }
        }
    }
}

void
init_ones_amatrix(pamatrix am, int n)
{
        for (int j=0;j<n;j++) {
    for (int i=0;i<n;i++) {
            am->a[i+j*n] = 1;
        }
    }
}

void
init_zeros_amatrix(pamatrix am, int n)
{
        for (int j=0;j<n;j++) {
    for (int i=0;i<n;i++) {
            am->a[i+j*n] = 0;
        }
    }
}
/*
void
init_copy_amatrix(pamatrix am, field *b)
{
    for (int i=0;i<am->rows;i++) {
        for (int j=0;j<am->cols;j++) {
            am->a[i*am->cols + j] = b[i*am->cols + j];
        }
    }
}
*/

pavector
init_avec_from_amat(pamatrix am, int basenum, int nu)
{
    pavector pv = new avector;

    pv->dim = basenum;

    pv->v = am->a + nu*basenum;

    return pv;
}

void
uninit_avec_from_amat(pavector pv)
{
    delete [] pv->v;
}

void
copy_sub_amatrix(bool atrans, pcamatrix a, pamatrix b)
{
	int lda = a->ld;
	int ldb = b->ld;

	int rows, cols;

	if (atrans) {
		rows = min<int>(a->cols, b->rows);
		cols = min<int>(a->rows, b->cols);
		for (int j=0;j<cols;j++) {
			for (int i=0;i<rows;i++) {
				b->a[i+j*ldb] = conj(a->a[j+i*lda]);
			}
		}
	}
	else {
		rows = min<int>(a->rows, b->rows);
		cols = min<int>(a->cols, b->cols);

		for (int j=0;j<cols;j++) {
			for (int i=0;i<rows;i++) {
				b->a[j*ldb+i] = a->a[j*lda+i];
			}
		}
	}
}

void
copy_amatrix(bool atrans, pcamatrix a, pamatrix b)
{
	if (atrans) {
		assert(a->rows == b->cols);
		assert(b->cols == a->rows);

		copy_sub_amatrix(true, a, b);
	}
	else {
		assert(a->rows == b->rows);
		assert(a->cols == b->cols);

		copy_sub_amatrix(false, a, b);
	}
}

void
identity_amatrix(pamatrix a)
{
	int rows = a->rows;
	int cols = a->cols;
	int lda = a->ld;
	for (int j=0;j<cols;j++) {
		for (int i=0;i<rows;i++) {
			a->a[j*lda+i] = 0;
		}
	}
	for (int i=0;i<cols && i<rows;i++) {
		a->a[i*lda+i] = 1.0;
	}
}

pamatrix
init_sub_amatrix(pamatrix a, pamatrix src, int rows, int roff,
				int cols, int coff)
{
	assert(a!=NULL);
	assert(src != NULL);
	assert(roff + rows <= src->rows);
	assert(coff + cols <= src->cols);

	a->a = src->a + roff + src->ld * coff;
	a->ld = src->ld;
	a->rows = rows;
	a->cols = cols;
	a->owner = src;

	return a;
}

void
clear_amatrix(pamatrix a)
{
	int lda = a->ld;
	for (int i=0;i<a->rows;i++) {
		for (int j=0;j<a->cols;j++) {
			a->a[i+j*lda] = 0.0;
		}
	}
}

double
normfrob_amatrix(pcamatrix a)
{
	double sum;
	int i;
	sum = 0;
	for (i=0;i<a->rows;i++) {
		sum+= REAL_SQR(h2_nrm2(a->cols, (MKL_Complex16 *)a->a+i*a->ld, i_one));
	}

	return REAL_SQRT(sum);
}


void
addmul_amatrix(field alpha, bool atrans, pcamatrix a, bool btrans, pcamatrix b, pamatrix c)
{
	int m, n, k;
	CBLAS_TRANSPOSE transa, transb;

	if (atrans) {
		m = a->cols;
		k = a->rows;
		transa = h2_ctrans;
	}
	else {
		m = a->rows;
		k = a->cols;
		transa = h2_ntrans;
	}
	if (btrans) {
		n = b->rows;
		transb = h2_ctrans;
	}
	else  {
		n = b->cols;
		transb = h2_ntrans;
	}

	h2_gemm(h2_col_major, transa, transb, m, n, k, &alpha, a->a, a->ld, b->a, b->ld, &f_one, c->a, c->ld);

//	cout << b->ld << "," << c->ld << endl;
}
/*
void
addmul_amatrix(field alpha, bool atrans, pcamatrix a, bool btrans,
	       pcamatrix b, pamatrix c)
{
  int      rows, cols, mid;
  pcfield   aa = a->a;
  pcfield   ba = b->a;
  pfield    ca = c->a;
  int lda = a->ld;
  int ldb = b->ld;
  int ldc = c->ld;
  int i, j, k;

  if (atrans) {
    if (btrans) {
      assert(a->cols <= c->rows);
      assert(b->rows <= c->cols);
      assert(a->rows == b->cols);

      rows = a->cols;
      mid = a->rows;
      cols = b->rows;

      for (i = 0; i < rows; i++) {
	for (k = 0; k < cols; k++) {
	  for (j = 0; j < mid; j++) {
	    ca[i + k * ldc] += alpha
	      * CONJ(aa[j + i * lda]) * CONJ(ba[k + j * ldb]);
	  }
	}
      }
    }
    else {
      assert(a->cols <= c->rows);
      assert(b->cols <= c->cols);
      assert(a->rows == b->rows);

      rows = a->cols;
      mid = a->rows;
      cols = b->cols;

      for (k = 0; k < cols; k++) {
	for (i = 0; i < rows; i++) {
	  for (j = 0; j < mid; j++) {
	    ca[i + k * ldc] +=
	      alpha * CONJ(aa[j + i * lda]) * ba[j + k * ldb];
	  }
	}
      }
    }
  }
  else {
    if (btrans) {
      assert(a->rows <= c->rows);
      assert(b->rows <= c->cols);
      assert(a->cols == b->cols);

      rows = a->rows;
      mid = a->cols;
      cols = b->rows;

      for (j = 0; j < mid; j++) {
	for (k = 0; k < cols; k++) {
	  for (i = 0; i < rows; i++) {
	    ca[i + k * ldc] +=
	      alpha * aa[i + j * lda] * CONJ(ba[k + j * ldb]);
	  }
	}
      }
    }
    else {
      assert(a->rows <= c->rows);
      assert(b->cols <= c->cols);
      assert(a->cols == b->rows);

      rows = a->rows;
      mid = a->cols;
      cols = b->cols;

      for (k = 0; k < cols; k++) {
	for (j = 0; j < mid; j++) {
	  for (i = 0; i < rows; i++) {
	    ca[i + k * ldc] += alpha * aa[i + j * lda] * ba[j + k * ldb];
	  }
	}
      }
    }
  }
}
*/

void
resize_amatrix(pamatrix a, int rows, int cols)
{
	assert(a->owner == NULL);

	if (rows != a->rows || cols != a->cols) {
		delete [] a->a;
		a->a = new field[rows*cols]();
		a->rows = rows;
		a->cols = cols;
		a->ld = rows;
	}
}

void
scale_amatrix(field alpha, pamatrix a)
{
	int lda = a->ld;
	for (int j=0;j<a->cols;++j) {
		h2_scal(a->rows, &alpha, a->a + j*lda, i_one);
	}
}

void
add_amatrix(field alpha, bool atrans, pcamatrix a, pamatrix b)
{
	int lda = a->ld;
	int ldb = b->ld;
	int rows = a->rows;
	int cols = a->cols;

	if (atrans) {
		assert(rows <= b->cols);
		assert(cols <= b->rows);

		for (int j=0;j<cols;++j) {
			h2_gerc(h2_col_major, i_one, cols, &alpha, &f_one, i_one, a->a+j*lda, i_one, b->a + j, ldb);
		}
	}
	else {
		assert(rows <= b->rows);
		assert(cols <= b->cols);

		for (int j=0;j<cols;++j) {
			h2_axpy(rows, &alpha, a->a + j*lda, i_one, b->a+j*ldb, i_one);
		}
	}
}

pamatrix
new_sub_amatrix(pamatrix src, int rows, int roff, int cols, int coff)
{
	pamatrix a = new amatrix;
	init_sub_amatrix(a, src, rows, roff, cols, coff);
	return a;
}

void
conj_amatrix(pamatrix a)
{
	for (int j=0;j<a->cols;++j) {
		for (int i=0;i<a->rows;++i) {
			a->a[i+j*a->ld] = conj(a->a[i+j*a->ld]);
		}
	}
}

void
output_amatrix(const string &filename, pamatrix a)
{
	ofstream fp(filename.c_str(), ios::out);

	for (int j=0;j<a->cols;++j) {
		for (int i=0;i<a->rows;++i) {
			fp << a->a[i+j*a->ld] << " ";
		}
		fp << endl;
	}
}
