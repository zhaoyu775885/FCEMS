#include "rkmatrix.h"
#include <math.h>

void
print_rkmatrix(prkmatrix r)
{
	print_amatrix(&(r->U));
	print_amatrix(&(r->V));
}

prkmatrix
new_rkmatrix(int rows, int cols, int k)
{
    prkmatrix rkm = new rkmatrix;
    init_rkmatrix(rkm, rows, cols, k);
    return rkm;
}


prkmatrix
init_rkmatrix(prkmatrix r, int rows, int cols, int k)
{
    init_amatrix(&(r->U), rows, k);
    init_amatrix(&(r->V), cols, k);
    r->k = k;
	return r;
}

void
uninit_rkmatrix(prkmatrix r)
{
	uninit_amatrix(&(r->U));
	uninit_amatrix(&(r->V));
}

void
del_rkmatrix(prkmatrix r)
{
	uninit_rkmatrix(r);

	//delete r;
}

void
copy_rkmatrix(bool atrans, pcrkmatrix a, prkmatrix b)
{
	if (atrans) {
		assert(a->U.rows == b->V.rows);
		assert(a->V.rows == b->U.rows);

		if (b->k != a->k)
			setrank_rkmatrix(b, a->k);
		copy_amatrix(false, &a->U, &b->V);
		copy_amatrix(false, &a->V, &b->U);
	}
	else {
		assert(a->U.rows == b->U.rows);
		assert(a->V.rows == b->V.rows);

		if (b->k != a->k)
			setrank_rkmatrix(b, a->k);
		copy_amatrix(false, &a->U, &b->U);
		copy_amatrix(false, &a->V, &b->V);
	}
}

void
setrank_rkmatrix(prkmatrix r, int k)
{
	resize_amatrix(&r->U, r->U.rows, k);
	resize_amatrix(&r->V, r->V.rows, k);
	r->k = k;
}

void
addeval_rkmatrix_avector(field alpha, prkmatrix rkm, char trans, pcavector x, field beta, pavector y)
{
    /// Never forget the U and V 's layout type
	// U,V are all ColMajor order!!!
    pavector utemp, vtemp;

    if (trans == CblasNoTrans) {
        for (int i=0;i<rkm->k;i++) {

            vtemp = init_avec_from_amat(&(rkm->V), rkm->V.rows, i); // V is Col Major
            utemp = init_avec_from_amat(&(rkm->U), rkm->U.rows, i); // U is Col Major

            field tmp = dot_prod_avec('c', vtemp, (pavector)x);
            add_avec(tmp*alpha, utemp, beta, y);

			delete vtemp;
			delete utemp;
        }
    }
	else if (trans == CblasConjTrans) {
		for (int i=0;i<rkm->k;i++) {
            vtemp = init_avec_from_amat(&(rkm->V), rkm->V.rows, i);
            utemp = init_avec_from_amat(&(rkm->U), rkm->U.rows, i);
			field tmp = dot_prod_avec('c', utemp, (pavector)x);
			add_avec(tmp*alpha, vtemp, beta, y);

			delete vtemp;
			delete utemp;
		}

    }
	else if (trans == CblasTrans) {
		for (int i=0;i<rkm->k;i++) {
            vtemp = init_avec_from_amat(&(rkm->V), rkm->V.rows, i);
            utemp = init_avec_from_amat(&(rkm->U), rkm->U.rows, i);
			field tmp = dot_prod_avec('u', utemp, (pavector)x);
			pavector vttemp = new_avector(rkm->V.rows);
			for (int i=0;i<vttemp->dim;i++) {
				vttemp->v[i] = conj(vtemp->v[i]);
			}
			add_avec(tmp, vttemp, beta, y);

			delete vtemp;
			delete utemp;
			uninit_avector(vttemp);
		}
	}
    else {
        cout << "Error! Clarify the Trans Type" << endl;
    }
}

pcrkmatrix
init_sub_rkmatrix(prkmatrix r, pcrkmatrix src, int rows, int roff,
				int cols, int coff)
{
	prkmatrix wsrc = (prkmatrix) src;
	int k = src->k;

	init_sub_amatrix(&r->U, &wsrc->U, rows, roff, k, 0);
	init_sub_amatrix(&r->V, &wsrc->V, cols, coff, k, 0);
	r->k = k;

	return r;
}


pcrkmatrix
new_sub_rkmatrix(pcrkmatrix src, int rows, int roff, int cols, int coff)
{
	prkmatrix r = new rkmatrix;
	return init_sub_rkmatrix(r, src, rows, roff, cols, coff);
}

void
resize_rkmatrix(prkmatrix r, int rows, int cols, int k)
{
	resize_amatrix(&r->U, rows, k);
	resize_amatrix(&r->V, cols, k);
	r->k = k;
}
