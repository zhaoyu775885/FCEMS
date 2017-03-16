#include "block.h"

int depth = 0;

bool
admissible_2_cluster(pcluster rc, pcluster cc, double* data)
{
	bool is_admissible;

	double eta = *(double *)data;

	double diamt, diams, dist;

	diamt = diams = dist = 0;

	for (int i=0;i<rc->dim;i++) {
		double temp;

		temp = rc->bmax[i] - rc->bmin[i];
		diamt += temp * temp;

		temp = cc->bmax[i] - cc->bmin[i];
		diams += temp * temp;

		temp = max<double>( max<double>(rc->bmin[i]-cc->bmax[i], cc->bmin[i]-rc->bmax[i]), 0);

		dist += temp * temp;
	}

	if (max<double>(diamt, diams) < eta * eta * dist) {
		is_admissible = true;
	}
	else {
		is_admissible = false;
	}

	return is_admissible;
}


bool
admissible_sphere_cluster(pccluster rc, pccluster cc, double *data)
{
    double eta;
    double diamt, diams, dist, a;
    int dim = rc->dim;

    assert(cc->dim == dim);

    diamt = 0.0;
    diams = 0.0;
    dist = 0.0;
    eta = *(double *) data;

    diamt = getdiam_2_cluster(rc);
    diams = getdiam_2_cluster(cc);

    for (int i = 0; i < dim; ++i) {
        double temp = ((rc->bmin[i] + rc->bmax[i]) - (cc->bmin[i] + cc->bmax[i])) * 0.5;
        dist += temp*temp;
    }
    dist = sqrt(dist) - 0.5 * (diamt + diams);

    a = max(diamt, diams);

    if (a <= eta * dist) {
        return true;
    }
    else {
        return false;
    }
}


pblock
new_block(pcluster rowcluster, pcluster clocluster, bool a, int rsons, int csons)
{
	pblock b = new block;
	b->rc = rowcluster;
	b->cc = clocluster;
	b->a = a;
	b->rsons = rsons;
	b->csons = csons;
	b->son = NULL;

	if (rsons > 0 && csons > 0) {
		b->son = new pblock [rsons*csons];
	}

	for (int i=0;i<rsons*csons;i++) {
		b->son[i] = NULL;
	}

	return b;
}

void
del_block(pblock b)
{
	for (int j=0;j<b->csons;j++) {
		for (int i=0;i<b->rsons;i++) {
			del_block(b->son[i+j*(b->rsons)]);
		}
	}

	delete [] b->son;
	delete b;
}

void
update_block(pblock b)
{
	int rsons = b->rsons;
	int csons = b->csons;
	int desc;

	desc = 1;
	for (int j=0;j<csons;j++) {
		for (int i=0;i<rsons;i++) {
			desc += b->son[i+j*rsons]->desc;
		}
	}

	b->desc = desc;
}

pblock
build_block(pcluster rc, pcluster cc, double *eta, admissible admis)
{
	pblock b;

	int rsons, csons;

	bool a = admis(rc, cc, eta);

	if (a == false) {
		if (rc->nsons == 0) {
			if (cc->nsons == 0) {
				rsons = 0;
				csons = 0;
				b = new_block(rc, cc, a, rsons, csons);
			}
			else {
				rsons = 1;
				csons = cc->nsons;
				b = new_block(rc, cc, a, rsons, csons);
				depth++;
				//cout << depth << endl;
				for (int i=0;i<csons;i++) {
					b->son[i] = build_block(rc, cc->son[i], eta, admis);
				}
				depth--;
			}
		}
		else {
			if (cc->nsons == 0) {
				rsons = rc->nsons;
				csons = 1;
				b = new_block(rc, cc, a, rsons, csons);
				depth++;
				//cout << depth << endl;
				for (int i=0;i<rsons;i++) {
					b->son[i] = build_block(rc->son[i], cc, eta, admis);
				}
				depth--;
			}
			else {
				rsons = rc->nsons;
				csons = cc->nsons;
				b = new_block(rc, cc, a, rsons, csons);
				depth++;
				//cout << depth << endl;
				for (int j=0;j<csons;j++) {
					for (int i=0;i<rsons;i++) {
						b->son[i+j*rsons] = build_block(rc->son[i], cc->son[j], eta, admis);
					}
				}
				depth--;
			}
		}
	}
	else {
		assert(a == true);
		rsons = 0;
		csons = 0;
		b = new_block(rc, cc, a, rsons, csons);
		//cout << "ad: depth=" << depth << endl;
	}

	update_block(b);

	return b;
}

#if 1

void
init_block_cluster_csp(pblock b)
{
    b->rc->csp_count = 0;
    if (b->son) {
        for (int i=0;i<b->rsons;++i) {
            for (int j=0;j<b->csons;++j) {
                init_block_cluster_csp(b->son[i+j*b->rsons]);
            }
        }
    }
}

void
count_csp(pblock b)
{
    b->rc->csp_count++;
    if (b->son) {
        for (int i=0;i<b->rsons;++i) {
            for (int j=0;j<b->csons;++j) {
                count_csp(b->son[i+j*b->rsons]);
            }
        }
    }
}

void
traverse_cluster(pccluster clu, int &csp)
{
    if (clu->csp_count > csp) csp = clu->csp_count;
    if (clu->son) {
        for (int i=0;i<clu->nsons;++i) {
            traverse_cluster(clu->son[i], csp);
        }
    }
}

int &get_csp(pblock b, int &csp)
{
    init_block_cluster_csp(b);
    count_csp(b);
    traverse_cluster(b->rc, csp);
    return csp;
}

#endif
