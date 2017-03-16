#ifndef AMATRIX_H
#define AMATRIX_H

#include "h2settings.h"
#include "realavector.h"
#include <assert.h>
#include <stdlib.h>
#include "blas.h"




typedef struct _amatrix amatrix;
typedef amatrix *pamatrix;
typedef const amatrix* pcamatrix;


struct _amatrix {
	field		*a;
    int         ld;
    int         rows;
    int         cols;
    void        *owner;
};


#include "avector.h"

void
print_amatrix(pamatrix a);

pamatrix
init_amatrix(pamatrix a, int rows, int cols);

void
uninit_amatrix(pamatrix a);

pamatrix
new_amatrix(int rows, int cols);

void
del_amatrix(pamatrix a);

pamatrix
new_zero_amatrix(int rows, int cols);

pamatrix
new_ones_amatrix(int rows, int cols);

void
addeval_amatrix_avector(field alpha, pamatrix am, char layout, char trans, pcavector x, field beta, pavector y);

void
init_eye_amatrix(pamatrix am, int n);

void
init_ones_amatrix(pamatrix am, int n);

void
init_zeros_amatrix(pamatrix am, int n);

void
init_copy_amatrix(pamatrix am, field *b);

pavector
init_avec_from_amat(pamatrix am, int basenum, int nu);

void
uninit_avec_from_amat(pavector pv);

void
copy_sub_amatrix(bool atrans, pcamatrix a, pamatrix b);

void
copy_amatrix(bool trans, pcamatrix a, pamatrix b);

void
identity_amatrix(pamatrix a);

pamatrix
init_sub_amatrix(pamatrix a, pamatrix src, int rows, int roff,
				int cols, int coff);

double
normfrob_amatrix(pcamatrix a);

void
addmul_amatrix(field alpha, bool atrans, pcamatrix a, bool btrans,
				pcamatrix b, pamatrix c);

void
clear_amatrix(pamatrix a);

void
resize_amatrix(pamatrix a, int rows, int cols);

void
scale_amatrix(field alpha, pamatrix a);

void
add_amatrix(field alpha, bool atrans, pcamatrix a, pamatrix b);

pamatrix
new_sub_amatrix(pamatrix src, int rows, int roff, int cols, int coff);

void
conj_amatrix(pamatrix a);

void
output_amatrix(const string &filename, pamatrix a);

#endif // _AMATRIX_H
