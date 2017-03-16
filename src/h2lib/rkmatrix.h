#ifndef RKMATRIX_H
#define RKMATRIX_H

#include "h2settings.h"
#include "amatrix.h"


typedef struct _rkmatrix rkmatrix;
typedef rkmatrix *prkmatrix;
typedef const rkmatrix *pcrkmatrix;

struct _rkmatrix {
    amatrix U;
    amatrix V;
    int k;
};


void
print_rkmatrix(prkmatrix r);

prkmatrix
new_rkmatrix(int rows, int cols, int k);

prkmatrix
init_rkmatrix(prkmatrix r, int rows, int cols, int k);

void
uninit_rkmatrix(prkmatrix r);

void
del_rkmatrix(prkmatrix r);

void
copy_rkmatrix(bool atrans, pcrkmatrix a, prkmatrix b);

void
setrank_rkmatrix(prkmatrix r, int k);

void
//rkmat_avec_prod(prkmatrix rkm, char trans, pcavector x, pavector y);
addeval_rkmatrix_avector(field alpha, prkmatrix rkm, char trans, pcavector x, field beta, pavector y);

pcrkmatrix
init_sub_rkmatrix(prkmatrix r, pcrkmatrix src, int rows, int roff,
				int cols, int coff);

pcrkmatrix
new_sub_rkmatrix(pcrkmatrix src, int rows, int roff, int cols, int coff);

void
resize_rkmatrix(prkmatrix r, int rows, int cols, int k);


#endif // _RKMATRIX_H

