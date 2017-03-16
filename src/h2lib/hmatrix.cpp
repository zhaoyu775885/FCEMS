#include "hmatrix.h"

void
print_hmatrix(phmatrix hm)
{
	if (hm->f != NULL) {
		print_amatrix(hm->f);
	}
	else if (hm->r != NULL) {
		print_rkmatrix(hm->r);
	}
}

phmatrix
init_hmatrix(phmatrix hm, pccluster rc, pccluster cc)
{
    hm->rc = rc;
    hm->cc = cc;

    hm->f = NULL;
    hm->r = NULL;

    hm->son = NULL;
    hm->rsons = 0;
    hm->csons = 0;

	hm->refs = 0;
	hm->desc = 0;

	return hm;
}

void
uninit_hmatrix(phmatrix hm)
{
	int rsons = hm->rsons;
	int csons = hm->csons;

	assert(hm->refs == 0);

	if (hm->son) {
		for (int j=0;j<csons;j++) {
			for (int i=0;i<rsons;i++) {
				unref_hmatrix(hm->son[i+j*rsons]);
			}
		}
		delete [] hm->son;
	}

	if (hm->f) {
		del_amatrix(hm->f);
	}

	if (hm->r) {
		del_rkmatrix(hm->r);
	}
}

phmatrix
new_hmatrix(pccluster rc, pccluster cc)
{
    phmatrix hm = new hmatrix;
    hm = init_hmatrix(hm, rc, cc);
    return hm;
}

phmatrix
new_rk_hmatrix(pccluster rc, pccluster cc, int k)
{
	phmatrix hm = new_hmatrix(rc, cc);

	hm->r = new_rkmatrix(rc->size, cc->size, k);

	hm->desc = 1;

	return hm;
}

phmatrix
new_full_hmatrix(pccluster rc, pccluster cc)
{
	phmatrix hm = new_hmatrix(rc, cc);

	hm->f = new_amatrix(rc->size, cc->size);

	hm->desc = 1;

	return hm;
}

phmatrix
new_super_hmatrix(pccluster rc, pccluster cc, int rsons, int csons)
{
	phmatrix hm = new_hmatrix(rc, cc);

	hm->rsons = rsons;
	hm->csons = csons;

	hm->son = new phmatrix[rsons*csons];

	for (int j=0;j<csons;j++) {
		for (int i=0;i<rsons;i++) {
			hm->son[i+j*rsons] = NULL;
		}
	}

	return hm;
}

void
del_hmatrix(phmatrix hm)
{
	uninit_hmatrix(hm);

	delete hm;
}

void
ref_hmatrix(phmatrix *ptr, phmatrix hm)
{
	if (*ptr) {
		unref_hmatrix(*ptr);
	}

	*ptr = hm;

	if (hm) {
		hm->refs++;
	}
}

void unref_hmatrix(phmatrix hm)
{
	assert(hm->refs > 0);

	hm->refs--;

	if (hm->refs == 0) {
		del_hmatrix(hm);
	}
}

phmatrix
build_from_block_hmatrix(pcblock b, int k)
{
	phmatrix h, h1;
	pblock b1;

	int rsons, csons;

	h = NULL;

	if (b->son) {
		rsons = b->rsons;
		csons = b->csons;

		h = new_super_hmatrix(b->rc, b->cc, rsons, csons);

		for (int j=0;j<csons;j++) {
			for (int i=0;i<rsons;i++) {
				b1 = b->son[i+j*rsons];
				h1 = build_from_block_hmatrix(b1, k);
				ref_hmatrix(h->son+j*rsons+i, h1);
			}
		}
	}
	else if (b->a) {
		//cout << ++rk_times << ": " << "build rk matrix" << endl;
		h = new_rk_hmatrix(b->rc, b->cc, k);
	}
	else {
		h = new_full_hmatrix(b->rc, b->cc);
	}

	return h;
}

bool
father_fmatrix(pchmatrix hm)
{
    int rsons = hm->rsons;
    int csons = hm->csons;
//    int count_rk = 0;
//    int count_full = 0;

    for (int j=0;j<csons;++j) {
        for (int i=0;i<rsons;++i) {
        #if 1
            if ( hm->son[i+j*rsons]->f == 0 )
                return false;
        #else
            if ( hm->son[i+j*rsons]->r != 0 )
                count_rk++;
            if ( hm->son[i+j*rsons]->f != 0 )
                count_full++;
        #endif // 0
        }
    }

//    if (count_rk + count_full == csons*rsons && count_rk<=1)
//        return true;

    return true;
}

bool
father_rkmatrix(pchmatrix hm)
{
    int rsons = hm->rsons;
    int csons = hm->csons;
//    int count_rk = 0;
//    int count_full = 0;

    for (int j=0;j<csons;++j) {
        for (int i=0;i<rsons;++i) {
#if 0
            if ( hm->son[i+j*rsons]->r != 0 )
                count_rk++;
            if ( hm->son[i+j*rsons]->f != 0 )
                count_full++;
#endif // 0
            if (hm->son[i+j*rsons]->r==0)
                return false;
        }
    }

    #if 0
    if (count_full<=1) {
        if (count_full+count_rk==csons*rsons)
            return true;
    }
    #endif // 0

    return true;
}

void
merge_fmatrix(phmatrix hm)
{
    int rsons = hm->rsons;
    int csons = hm->csons;

    if (hm->son != 0) {
        for (int j=0;j<csons;++j) {
            for (int i=0;i<rsons;++i) {
                merge_fmatrix(hm->son[i+j*rsons]);
            }
        }
        if ( father_fmatrix(hm) ) {
            for (int j=0;j<csons;++j) {
                for (int i=0;i<rsons;++i) {
                    unref_hmatrix(hm->son[i+j*rsons]);
                }
            }
            hm->f = new_amatrix(hm->rc->size, hm->cc->size);
            hm->r = 0;
            hm->son = 0;
            hm->rsons = 1;
            hm->csons = 1;
            hm->desc = 1;
        }
    }
}

void
merge_rkmatrix(phmatrix hm)
{
    int rsons = hm->rsons;
    int csons = hm->csons;

    if (hm->son != 0) {
        for (int j=0;j<csons;++j) {
            for (int i=0;i<rsons;++i) {
                merge_rkmatrix(hm->son[i+j*rsons]);
            }
        }
        if ( father_rkmatrix(hm) ) {
            for (int j=0;j<csons;++j) {
                for (int i=0;i<rsons;++i) {
                    unref_hmatrix(hm->son[i+j*rsons]);
                }
            }
            hm->r = new_rkmatrix(hm->rc->size, hm->cc->size, 0);
            hm->f = 0;
            hm->son = 0;
            hm->rsons = 1;
            hm->csons = 1;
            hm->desc = 1;
        }
    }
}

void
addeval_hmatrix_avector(field alpha, phmatrix hm, char trans, pcavector x, field beta, pavector y)
{
	pavector  xp, yp;
	avector   xtmp, ytmp;
	int      i, ip;

	assert(x->dim == hm->cc->size);
	assert(y->dim == hm->rc->size);

	/* Permutation of x */
	xp = init_avector(&xtmp, x->dim);
	for (i = 0; i < xp->dim; i++) {
		ip = hm->cc->idx[i];
		assert(ip < x->dim);
		xp->v[i] = x->v[ip];
	}

	/* Permutation of y */
	yp = init_avector(&ytmp, y->dim);
	for (i = 0; i < yp->dim; i++) {
		ip = hm->rc->idx[i];
		assert(ip < y->dim);
		yp->v[i] = y->v[ip];
	}

	/* Matrix-vector multiplication */
	//hmat_avec_prod(hm, trans, xp, yp);
	fastaddeval_hmatrix_avector(alpha, hm, trans, xp, beta, yp);

	/* Reverse permutation of y */
	for (i = 0; i < yp->dim; i++) {
		ip = hm->rc->idx[i];
		assert(ip < y->dim);
		y->v[ip] = yp->v[i];
	}

	uninit_avector(yp);
	uninit_avector(xp);
}

void
mvm_hmatrix_avector(phmatrix hm, char trans, pcavector x, pavector y)
{
	int rsons = hm->rsons;
	int csons = hm->csons;


	if(hm->f) {
		addeval_amatrix_avector(1.0, hm->f, h2_col_major, trans, x, 1.0, y);
	}
	else if (hm->r) {
		addeval_rkmatrix_avector(1.0, hm->r, trans, x, 1.0, y);
	}
	else {
		assert(hm->r==NULL && hm->f==NULL);
		if (h2_ntrans == trans) {
			int xoffset = 0;
			for (int j=0;j<csons;++j) {
				pavector xsub = new_sub_avector(x, hm->son[j*rsons]->cc->size, xoffset);
				int yoffset = 0;
				for (int i=0;i<hm->rsons;++i) {
					pavector ysub = new_sub_avector(y, hm->son[i]->rc->size, yoffset);
					mvm_hmatrix_avector(hm->son[i+j*rsons], trans, xsub, ysub);
					yoffset += hm->son[i]->rc->size;
					uninit_avector(ysub);
				}
				assert(yoffset == hm->rc->size);
				xoffset += hm->son[j*rsons]->cc->size;
				uninit_avector(xsub);
			}
			assert(xoffset == hm->cc->size);
		}
		else {
			int xoffset = 0;
			for (int i=0;i<rsons;++i) {
				pavector xsub = new_sub_avector(x, hm->son[i]->rc->size, xoffset);
				int yoffset = 0;
				for (int j=0;j<csons;++j) {
					pavector ysub = new_sub_avector(y, hm->son[j*rsons]->cc->size, yoffset);
					mvm_hmatrix_avector(hm->son[i+j*rsons], trans, xsub, ysub);
					yoffset += hm->son[j*rsons]->cc->size;
					uninit_avector(ysub);
				}
				assert(yoffset == hm->cc->size);
				xoffset += hm->son[i]->rc->size;
				uninit_avector(xsub);
			}
			assert(xoffset == hm->rc->size);
		}
	}
}

void
fastaddeval_hmatrix_avector(field alpha, phmatrix hm, char trans, pcavector x, field beta, pavector y)
{
	pavector newy = new_avector(y->dim);
	mvm_hmatrix_avector(hm, trans, x, newy);
	add_avec(alpha, newy, beta, y);

	/*
	for (int i=0;i<y->dim;i++) {
		y->v[i] = alpha * newy->v[i] + beta *y->v[i];
	}
	*/

#if 0
    if (hm->f) {
        addeval_amatrix_avector(alpha, hm->f, h2_col_major, trans, x, 0, y);
    }
    else if (hm->r) {
        addeval_rkmatrix_avector(alpha, hm->r, trans, x, 0, y);
    }
    else {
		assert(hm->r==NULL && hm->f==NULL);

        int yoffset = 0;
        if (CblasNoTrans == trans) {
            for (int i=0;i<hm->rsons;i++) {
                pavector ysub = new_sub_avector(y, hm->son[i]->rc->size, yoffset);
                pavector ytemp = new_avector(hm->son[i]->rc->size);
                int xoffset = 0;
                for (int j=0;j<hm->csons;j++) {
                    pavector xsub = new_sub_avector((pavector)x, hm->son[j*hm->rsons+i]->cc->size, xoffset);
                    fastaddeval_hmatrix_avector(alpha, hm->son[j*hm->rsons+i], trans, xsub, beta, ytemp);
                    add_avec(alpha, ytemp, 1, ysub);
                    xoffset += hm->son[j*hm->rsons+i]->cc->size;
					uninit_avector(xsub);
					// in urgent need to be optimized
                }
                yoffset += hm->son[i]->rc->size;
				uninit_avector(ysub);
				uninit_avector(ytemp);
                //free(ysub);
            }
        }
        else if (CblasTrans == trans || CblasConjTrans == trans) {
            for (int i=0;i<hm->csons;i++) {
				pavector ysub = new_sub_avector(y, hm->son[i*hm->rsons]->cc->size, yoffset);
				pavector ytemp = new_avector(hm->son[i*hm->rsons]->cc->size);
				int xoffset = 0;
				for (int j=0;j<hm->rsons;j++) {
					pavector xsub = new_sub_avector((pavector)x, hm->son[i*hm->rsons+j]->rc->size, xoffset);
					fastaddeval_hmatrix_avector(alpha, hm->son[i*hm->rsons+j], trans, xsub, beta, ytemp);
					add_avec(alpha, ytemp, 1, ysub);
					xoffset += hm->son[i*hm->rsons+j]->rc->size;
					uninit_avector(xsub);
				}
				yoffset += hm->son[i*hm->rsons]->cc->size;
				uninit_avector(ysub);
				uninit_avector(ytemp);
            }
        }
        else {
            cout<<"Error! Clarify the Trans Type\n"<<endl;
        }
    }
#endif
}

#if 1
int
hmat_cg(phmatrix hm, pavector b, pavector x)
{
    if (hm->rc->size != b->dim) {
        return -1;
    }
    else {
        double ep = 1e-2;
        pavector r = new_avector(b->dim);
        pavector p = new_avector(b->dim);
        pavector q = new_avector(b->dim);
        pavector tempq = new_avector(b->dim);
        pavector newb = new_avector(b->dim);

		pavector bp, xp;
		avector btmp, xtmp;
		int ip;
	/* Permutation of b */
		bp = init_avector(&btmp, b->dim);
		for (int i = 0; i < bp->dim; i++) {
			ip = hm->cc->idx[i];
			assert(ip < x->dim);
			bp->v[i] = b->v[ip];
		}
	/* Permutation of y */
		xp = init_avector(&xtmp, x->dim);

        fastaddeval_hmatrix_avector(1.0, hm, CblasConjTrans, bp, 0.0, newb);

        for (int i=0;i<bp->dim;i++) {
            r->v[i] = newb->v[i];
            xp->v[i] = 0;
        }

        int i;
        field alpha, beta;
        field rr1;
        field rr;
        for (i=0;i<40000;i++) {

            rr = dot_prod_avec('c', r, r);
            if (i == 0) {
                for (int j=0;j<bp->dim;j++) {
                    p->v[j] = r->v[j];
                }
            }
            else {
                beta = rr / rr1;
                add_avec(1, r, beta, p);
            }

            fastaddeval_hmatrix_avector(1.0, hm, CblasNoTrans, p, 0.0, tempq);
            fastaddeval_hmatrix_avector(1.0, hm, CblasConjTrans, tempq, 0.0, q);

            alpha = rr / dot_prod_avec('c', q, p);

            add_avec(alpha, p, 1, xp);
            add_avec(-alpha, q, 1, r);

            rr1 = rr;

            if (LAPACKE_zlange(LAPACK_ROW_MAJOR, 'f', 1, r->dim, (MKL_Complex16 *)r->v, r->dim) < ep) {
                break;
            }
        }
			/* Reverse permutation of y */
			for (int i = 0; i < xp->dim; i++) {
				ip = hm->rc->idx[i];
				assert(ip < xp->dim);
				x->v[ip] = xp->v[i];
			}
			uninit_avector(bp);
			uninit_avector(xp);
        return i;
    }
}
#else
int
hmat_cg(phmatrix hm, pavector b, pavector x)
{
    if (hm->rc->size != b->dim) {
        return -1;
    }
    else {
        double ep = 1e-2;
        pavector r = new_avector(b->dim);
        pavector p = new_avector(b->dim);
        pavector q = new_avector(b->dim);
        pavector tempq = new_avector(b->dim);
        pavector newb = new_avector(b->dim);

        addeval_hmatrix_avector(f_one, hm, CblasConjTrans, b, newb);

        for (int i=0;i<b->dim;i++) {
            r->v[i] = newb->v[i];
            x->v[i] = 0;
        }

        int i;
        field alpha, beta;
        field rr1;
        field rr;
        for (i=0;i<10001;i++) {

            rr = dot_prod_avec('c', r, r);

            if (i == 0) {
                for (int j=0;j<b->dim;j++) {
                    p->v[j] = r->v[j];
                }
            }
            else {
                beta = rr / rr1;
                add_avec(1, r, beta, p);
            }

            addeval_hmatrix_avector(f_one, hm, CblasNoTrans, p, tempq);

            addeval_hmatrix_avector(f_one, hm, CblasConjTrans, tempq, q);


            alpha = rr / dot_prod_avec('c', q, p);

            add_avec(alpha, p, 1, x);
            add_avec(-alpha, q, 1, r);

            rr1 = rr;

            if (LAPACKE_zlange(LAPACK_ROW_MAJOR, 'f', 1, r->dim, (MKL_Complex16 *)r->v, r->dim) < ep) {
                break;
            }
        }
        return i;
    }
}

#endif

void
update_hmatrix(phmatrix hm)
{
	int desc;
	int rsons, csons;

	desc = 1;

	if (hm->son) {
		rsons = hm->rsons;
		csons = hm->csons;

		for (int j=0;j<csons;++j) {
			for (int i=0;i<rsons;++i) {
				desc += hm->son[i+j*rsons]->desc;
			}
		}
	}

	hm->desc = desc;
}


void
clear_hmatrix(phmatrix hm)
{
	if (hm->son) {
		int rsons = hm->rsons;
		int csons = hm->csons;

		for (int j=0;j<csons;j++) {
			for (int i=0;i<rsons;i++) {
				clear_hmatrix(hm->son[i+j*rsons]);
			}
		}
	}
	else if(hm->r) {
		setrank_rkmatrix(hm->r, 0);
	}
	else {
		assert(hm->f);
		clear_amatrix(hm->f);
	}
}

phmatrix
clone_hmatrix(pchmatrix src)
{
	const int rsons = src->rsons;
	const int csons = src->csons;

	phmatrix hm1, hm;

	if(src->son != NULL) {
		hm = new_super_hmatrix(src->rc, src->cc, rsons, csons);
		for (int j=0;j<csons;++j) {
			for (int i=0;i<rsons;++i) {
				hm1 = clone_hmatrix(src->son[i+j*rsons]);
				ref_hmatrix(hm->son +i +j*rsons, hm1);
			}
		}
		update_hmatrix(hm);
	}
	else if (src->r != NULL) {
		hm = new_rk_hmatrix(src->rc, src->cc, src->r->k);
		copy_amatrix(false, &src->r->U, &hm->r->U);
		copy_amatrix(false, &src->r->V, &hm->r->V);
	}
	else {
		assert(src->f != NULL);
		hm = new_full_hmatrix(src->rc, src->cc);
		copy_amatrix(false, src->f, hm->f);
	}

	update_hmatrix(hm);

	return hm;
}

phmatrix
clonestructure_hmatrix(pchmatrix src)
{
	const int rsons = src->rsons;
	const int csons = src->csons;

	phmatrix hm1, hm;

	if(src->son != NULL) {
		hm = new_super_hmatrix(src->rc, src->cc, rsons, csons);
		for (int j=0;j<csons;++j) {
			for (int i=0;i<rsons;++i) {
				hm1 = clonestructure_hmatrix(src->son[i+j*rsons]);
				ref_hmatrix(hm->son +i +j*rsons, hm1);
			}
		}
		update_hmatrix(hm);
	}
	else if (src->r != NULL) {
		hm = new_rk_hmatrix(src->rc, src->cc, src->r->k);
	}
	else {
		assert(src->f != NULL);
		hm = new_full_hmatrix(src->rc, src->cc);
	}

	update_hmatrix(hm);

	return hm;
}

unsigned int
count_hmatrix_mem(pchmatrix hm)
{
    unsigned int mem_size(0);

    if (hm->f) {
        mem_size = hm->f->cols * hm->f->rows;
    }
    else if (hm->r) {
        mem_size = hm->r->U.cols*hm->r->U.rows + hm->r->V.cols*hm->r->V.rows;
    }
    else if (hm->son) {
        for (int j=0;j<hm->csons;++j) {
            for (int i=0;i<hm->rsons;++i) {
                mem_size += count_hmatrix_mem(hm->son[i+j*hm->rsons]);
            }
        }
    }

    return mem_size;
}

