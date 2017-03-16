#include "avector.h"

pavector
init_avector(pavector v, int dim)
{
	assert (v != NULL);

    v->dim = dim;
    v->v = (dim > 0) ? new field[dim]() : NULL;
#if 0
	for (int i=0;i<dim;i++) {
		v->v[i] = 0;
	}
#endif
	v->owner = NULL;

	return v;
}

pavector
init_sub_avector(pavector sub, pcavector src, int dim, int off)
{
	assert(sub != NULL);
	assert(src != NULL);
	assert(off+dim <= src->dim);

    sub->v = src->v + off;
    sub->dim = dim;
	sub->owner = src;

    return sub;
}

pavector
init_zero_avector(pavector v, int dim)
{
	assert(v != NULL);

	init_avector(v, dim);
	clear_avector(v);

    return v;
}

pavector
init_row_avector(pavector v, pamatrix src, int row)
{
	int lda;

	assert(v != NULL);
	assert(src != NULL);
	assert(row < src->rows);

	lda = src->ld;

	v->v = src->a + row * lda;
	v->dim = src->cols;
	v->owner = src;

	return v;
}

pavector
init_column_avector(pavector v, pamatrix src, int row)
{
	int lda;

	assert(v != NULL);
	assert(src != NULL);
	assert(row < src->rows);

	lda = src->ld;

	v->v = src->a + row * lda;
	v->dim = src->cols;
	v->owner = src;

	return v;
}

pavector
init_pointer_avector(pavector v, pfield src, int dim)
{
	assert (v != NULL);
	assert (dim == 0 || src != NULL);

	v->v = src;
	v->dim = dim;
	v->owner = src;

	return v;
}

void
uninit_avector(pavector v)
{
	if (! v->owner) {
		delete [] v->v;
	}
}

pavector
new_avector(int dim)
{
    pavector avec = new avector;

    init_avector(avec, dim);

    return avec;
}

pavector
new_sub_avector(pcavector src, int dim, int off)
{
	pavector v = new avector;

	v = init_sub_avector(v, src, dim, off);

	return v;
}

pavector
new_zeros_avector(int dim)
{
    pavector avec = new avector;

    init_zero_avector(avec, dim);

    return avec;
}

pavector
new_pointer_avector(pfield src, int dim)
{
	pavector v = new avector;

	init_pointer_avector(v, src, dim);

	return v;
}

void
del_avector(pavector v)
{
	uninit_avector(v);
	delete v;
}

void
resize_avector(pavector v, int dim)
{
	assert(v->owner == NULL);

	if (dim != v->dim) {
		delete [] v->v;
		v->v = new field[dim];
		v->dim = dim;
	}
}


pavector
new_ones_avector(int dim)
{
    pavector avec = new avector;
	//allocmem(sizeof(avector));

    init_avector(avec, dim);

    for (int i=0;i<dim;i++) {
        avec->v [i] = 1;
    }

    return avec;
}

void
shrink_avector(pavector v, int dim)
{
	assert(dim <= v->dim);

	v->dim = dim;
}



void
add_avec(field alpha, pcavector x, field beta, pavector y)
{
	assert(y->dim >= x->dim);

	h2_scal(y->dim, &beta, y->v, 1);
	h2_axpy(y->dim, &alpha, x->v, 1, y->v, 1);
}


field
dot_prod_avec(char type, pavector x, pavector y)
{
    field ans;
    if ('u' == type) {
        cblas_zdotu_sub(x->dim, x->v, 1, y->v, 1, &ans);
    }
    else if ('c' == type) {
        cblas_zdotc_sub(x->dim, x->v, 1, y->v, 1, &ans);
    }
    else {
        cout << "Error! Clarify the dot product type!\n" << endl;
    }

    return ans;
}

void
clear_avector(pavector v)
{
	for (int i=0;i<v->dim;i++) {
		v->v[i] = 0;
	}
}

void fill_avector(pavector v, field x)
{
	for (int i=0;i<v->dim;i++) {
		v->v[i] = x;
	}

}
