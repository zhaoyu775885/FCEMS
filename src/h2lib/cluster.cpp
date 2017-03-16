#include "cluster.h"
#include <math.h>
#include <assert.h>
#include <queue>
#include "amatrix.h"
#include "eigensolvers.h"

using namespace std;

int cluster_depth = 0;

pcluster
new_cluster(int size, int* idx, int nsons, int dim)
{
	pcluster p = new cluster;

	p->size = size;
	p->idx = idx;
	p->nsons = nsons;
	p->dim = dim;
	p->bmin = new double [dim];
	p->bmax = new double [dim];

	if (nsons > 0) {
		p->son = new pcluster [nsons];

		for (int i=0;i<nsons;i++) {
			p->son[i] = NULL;
		}
	}
	else {
		p->son = NULL;
	}
	p->desc = 1;

	return p;
}


void
update_cluster(pcluster pc)
{
    for (int i=0;i<pc->nsons;++i) {
        pc->desc += pc->son[i]->desc;
    }
}

void
del_cluster(pcluster p)
{
	if (p->nsons > 0) {
		for (int i=0;i<p->nsons;i++) {
			del_cluster(p->son[i]);
		}
		delete [] p->son;
	}
	delete [] p->bmin;
	delete [] p->bmax;
	delete p;
}

pcluster
build_regular_cluster(pclustergeometry cf, int size, int* label_tree, int leaf_threshold, int direction)
{
	pcluster p;

	int newd;
	int size0, size1;

	double lower, middle, upper;

	update_bbox_for_clustergeometry(cf, size, label_tree);

	if (size > leaf_threshold) {
		size0 = size1 = 0;

		/* 用以紧缩计算空间，保证后续的空间二分操作 */
		/*  避免在z方向进行二分，因为z向尺度很小   */
		newd = (direction+1) % cf->dim;

		middle = cf->hmax[direction] - cf->hmin[direction];

		if ( middle > 0.0 ) {
			middle = (cf->hmax[direction] + cf->hmin[direction]) / 2.0;
			for (int i=0;i<size;i++) {
				if (cf->x[label_tree[i]][direction] < middle) {
					swap<int>(label_tree[i], label_tree[size0]);
					size0++;
				}
				else {
					size1++;
				}
				// 分均分配点数是不合理的，因为判断admissible是根据空间信息来弄的
			}
			if (size0 > 0) {
				if (size1 > 0) {
					p = new_cluster(size, label_tree, 2, cf->dim);

					lower = cf->hmin[direction];
					upper = cf->hmax[direction];
					cf->hmax[direction] = middle;
					p->son[0] = build_regular_cluster(cf, size0, label_tree, leaf_threshold, newd);

					cf->hmin[direction] = middle;
					cf->hmax[direction] = upper;
					p->son[1] = build_regular_cluster(cf, size1, label_tree+size0, leaf_threshold, newd);

					cf->hmin[direction] = lower;
					update_bbox_for_cluster(p);
					/* 递归更新cluster tree的bounding box信息 */
				}
				else {
					p = new_cluster(size, label_tree, 1, cf->dim);
					upper = cf->hmax[direction];
					cf->hmax[direction] = middle;
					p->son[0] = build_regular_cluster(cf, size, label_tree, leaf_threshold, newd);
					cf->hmax[direction] = upper;
					update_bbox_for_cluster(p);
				}
			}
			else {
				p = new_cluster(size, label_tree, 1, cf->dim);
				lower = cf->hmin[direction];
				cf->hmin[direction] = middle;
				p->son[0] = build_regular_cluster(cf, size, label_tree, leaf_threshold, newd);
				cf->hmin[direction] = lower;
				update_bbox_for_cluster(p);
			}
		}
		else {
			assert(fabs(middle) <= 1e-10);
			p = new_cluster(size, label_tree, 1, cf->dim);
			p->son[0] = build_regular_cluster(cf, size, label_tree, leaf_threshold, newd);
			update_bbox_for_cluster(p);
		}
	}
	else {
		p = new_cluster(size, label_tree, 0, cf->dim);
		update_sbox_for_cluster(cf, p);
		/* 用特征点信息来实现二分，最终还是用全体节点更新support的bounding box信息*/
	}
	update_cluster(p);

	return p;
}

pcluster
build_adaptive_cluster(pclustergeometry cf, int size, int *idx, int clf)
{
    pcluster  t;

    int direction;
    int size0, size1;
    int i, j;
    double a, m;

    if (size > clf) {
        update_bbox_for_clustergeometry(cf, size, idx);

        /* compute the direction of partition */
        direction = 0;
        a = cf->hmax[0] - cf->hmin[0];

        for (j = 1; j < cf->dim; j++) {
            m = cf->hmax[j] - cf->hmin[j];
            if (a < m) {
                a = m;
                direction = j;
            }
        }

        /* build sons */
        if (a > 0.0) {
//            m = (cf->hmax[direction] + cf->hmin[direction]) / 2.0;

            double pos(0);
            for (int i=0;i<size;++i) {
                pos += cf->x[idx[i]][direction];
            }
            m = pos / size;

            size0 = 0;
            size1 = 0;

            for (i = 0; i < size; i++) {
                if (cf->x[idx[i]][direction] < m) {
                    j = idx[i];
                    idx[i] = idx[size0];
                    idx[size0] = j;
                    size0++;
                }
                else {
                    size1++;
                }
            }

            t = new_cluster(size, idx, 2, cf->dim);

            t->son[0] = build_adaptive_cluster(cf, size0, idx, clf);
            t->son[1] = build_adaptive_cluster(cf, size1, idx + size0, clf);

            update_bbox_for_cluster(t);
        }
        else {
            assert(a == 0.0);
            t = new_cluster(size, idx, 0, cf->dim);
            update_sbox_for_cluster(cf, t);
        }
    }
    else {
        update_bbox_for_clustergeometry(cf, size, idx);
        t = new_cluster(size, idx, 0, cf->dim);
        update_sbox_for_cluster(cf, t);
    }

    update_cluster(t);

    return t;
}

pcluster
build_pca_cluster(pclustergeometry cf, int size, int * idx, int clf)
{
    const int dim = cf->dim;

    pamatrix  C, Q;
    pavector  v;
    prealavector lambda;

    field *x, *y;
    field w;
    int i, j, k, size0, size1;

    pcluster  t;

    size0 = 0;
    size1 = 0;

    if (size > clf) {
        x = new field [dim];
        y = new field [dim];

        /* determine weight of current cluster */
        w = 0.0;
        for (i = 0; i < size; ++i) {
            w += cf->w[idx[i]];
        }
        w = 1.0 / w;

        for (j = 0; j < dim; ++j) {
            x[j] = 0.0;
        }

        /* determine center of mass */
        for (i = 0; i < size; ++i) {
            for (j = 0; j < dim; ++j) {
                x[j] += cf->w[idx[i]] * cf->x[idx[i]][j];
            }
        }

        for (j = 0; j < dim; ++j) {
            x[j] *= w;
        }

        C = new_zero_amatrix(dim, dim);
        Q = new_zero_amatrix(dim, dim);
        lambda = new_realavector(dim);

        /* setup covariance matrix */
        for (i = 0; i < size; ++i) {

            for (j = 0; j < dim; ++j) {
                y[j] = cf->x[idx[i]][j] - x[j];
            }

            for (j = 0; j < dim; ++j) {

                for (k = 0; k < dim; ++k) {
                    C->a[j + k * C->ld] += cf->w[idx[i]] * y[j] * y[k];
                }
            }
        }

        /* get eigenvalues and eigenvectors of covariance matrix */
        eig_amatrix(C, lambda, Q);

        /* get eigenvector from largest eigenvalue */
        v = new_avector(0);
        init_column_avector(v, Q, dim - 1);

        /* separate cluster with v as separation-plane */
        for (i = 0; i < size; ++i) {
            /* x_i - X */
            for (j = 0; j < dim; ++j) {
                y[j] = cf->x[idx[i]][j] - x[j];
            }

            /* <y,v> */
            w = 0.0;
            for (j = 0; j < dim; ++j) {
                w += y[j] * v->v[j];
            }

            if (abs(w) >= 0.0) {
                j = idx[i];
                idx[i] = idx[size0];
                idx[size0] = j;
                size0++;
            }
            else {
                size1++;
            }
        }

        assert(size0 + size1 == size);

        del_amatrix(Q);
        del_amatrix(C);
        del_realavector(lambda);
        del_avector(v);
        delete [] x;
        delete [] y;

        /* recursion */
        if (size0 > 0) {
            if (size1 > 0) {
                t = new_cluster(size, idx, 2, cf->dim);

                t->son[0] = build_pca_cluster(cf, size0, idx, clf);
                t->son[1] = build_pca_cluster(cf, size1, idx + size0, clf);

                update_bbox_for_cluster(t);
            }
            else {
                t = new_cluster(size, idx, 1, cf->dim);
                t->son[0] = build_pca_cluster(cf, size0, idx, clf);

                update_bbox_for_cluster(t);
            }
        }
        else {
            assert(size1 > 0);
            t = new_cluster(size, idx, 1, cf->dim);
            t->son[0] = build_pca_cluster(cf, size1, idx, clf);

            update_bbox_for_cluster(t);
        }
    }
    else {
        t = new_cluster(size, idx, 0, cf->dim);
        update_sbox_for_cluster(cf, t);
    }

    update_cluster(t);

    return t;
}

pcluster
build_cluster(pclustergeometry cf, int size, int * idx, int clf, clustermode mode)
{
    pcluster  t;

    if (mode == H2_ADAPTIVE) {
        t = build_adaptive_cluster(cf, size, idx, clf);
    }
    else if (mode == H2_REGULAR) {
        update_bbox_for_clustergeometry(cf, size, idx);
        t = build_regular_cluster(cf, size, idx, clf, 0);
    }
    else if (mode == H2_PCA) {
        t = build_pca_cluster(cf, size, idx, clf);
    }
//    else {
//        assert(mode == H2_SIMSUB);
//        update_point_bbox_clustergeometry(cf, size, idx);
//        t = build_simsub_cluster(cf, size, idx, clf);
//    }

    return t;
}

void
print_cluster_dim(pccluster t)
{
    queue<pccluster> preque;
    preque.push(t);
    cout << time << " layer: " << t->size << endl;
//    pccluster tmp = preque.front();
    while (!preque.empty()) {
        int n = preque.size();
        for (int i=0;i<n;++i) {
            pccluster tmp = preque.front();
            cout << tmp->size << " ";
            for (int j=0;j<tmp->nsons;++j) {
                preque.push(tmp->son[j]);
            }
            preque.pop();
        }
        cout << endl;
    }
}

int
getdepth_cluster(pccluster t)
{
    int depth(0);
    int son_depth(0);

    if (t->son) {
        for (int i=0;i<t->nsons;++i) {
            son_depth = getdepth_cluster(t->son[i]);
            depth = max(depth, son_depth);
        }
        depth++;
    }

    return depth;
}

double
getdiam_2_cluster(pccluster t)
{
    double diam2(0);

    for (int i=0;i<t->dim;++i)
        diam2 += (t->bmax[i]-t->bmin[i])*(t->bmax[i]-t->bmin[i]);

    return sqrt(diam2);
}
