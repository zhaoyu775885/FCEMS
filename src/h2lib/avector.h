#ifndef AVECTOR_H
#define AVECTOR_H


#include <mkl.h>
#include <assert.h>
#include "h2settings.h"

typedef struct _avector avector;
typedef avector* pavector;
typedef const avector* pcavector;

struct _avector {
    int dim;
    field *v;
//    double *v;
	const void *owner;
};

#include "amatrix.h"


pavector
init_avector(pavector v, int dim);

pavector
init_sub_avector(pavector sub, pcavector src, int size, int offset);

pavector
init_zero_avector(pavector v, int dim);

pavector
init_row_avector(pavector v, pamatrix src, int row);

pavector
init_pointer_avector(pavector v, pfield src, int dim);

void
uninit_avector(pavector v);

pavector
new_avector(int dim);

pavector
new_sub_avector(pcavector src, int dim, int off);

pavector
new_zeros_avector(int dim);

pavector
new_pointer_avector(pfield src, int dim);

void
del_avector(pavector v);

void
resize_avector(pavector v, int dim);

void
shrink_avector(pavector v, int dim);

pavector
new_ones_avector(int dim);

void
add_avec(field alpha, pcavector x, field beta, pavector y);

field
dot_prod_avec(char type, pavector x, pavector y);


void
clear_avector(pavector v);

void
fill_avector(pavector v, field x);

#include "amatrix.h"

pavector
init_column_avector(pavector v, pamatrix src, int row);
#endif
