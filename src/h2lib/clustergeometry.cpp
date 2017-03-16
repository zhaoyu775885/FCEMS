/* ------------------------------------------------------------
 This is the file "clustergeometry.c" of the H2Lib package.
 All rights reserved, Knut Reimer 2009
 ------------------------------------------------------------ */

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "clustergeometry.h"

/* ------------------------------------------------------------
 Constructors and destructors
 ------------------------------------------------------------ */

pclustergeometry
new_clustergeometry(int dim, int nidx)
{
    pclustergeometry cf = new clustergeometry;

    cf->dim = dim;
    cf->nidx = nidx;
    cf->x = new double* [nidx];
    cf->hmin = new double [3];
    cf->hmax = new double [3];
    cf->smin = new double* [nidx];
    cf->smax = new double* [nidx];
    cf->w = new double [nidx];

    return cf;
}

void
del_clustergeometry(pclustergeometry cf)
{
  free(cf->buf);
  free(cf->smax);
  free(cf->smin);
  free(cf->w);
  free(cf->hmax);
  free(cf->hmin);
  free(cf->x);
  free(cf);
}



/* ------------------------------------------------------------
 Auxiliary routines
 ------------------------------------------------------------ */

void
update_bbox_for_clustergeometry(pclustergeometry cf, int size, int * idx)
{
    for (int j = 0; j < cf->dim; j++) {
        cf->hmin[j] = cf->x[idx[0]][j];
        cf->hmax[j] = cf->x[idx[0]][j];
    }

    for (int i = 1; i < size; i++) {
        for (int j = 0; j < cf->dim; j++) {
            if (cf->x[idx[i]][j] < cf->hmin[j]) {
                cf->hmin[j] = cf->x[idx[i]][j];
            }
            if (cf->x[idx[i]][j] > cf->hmax[j]) {
                cf->hmax[j] = cf->x[idx[i]][j];
            }
        }
    }
}

void
update_sbox_for_cluster(pclustergeometry cf, pcluster t)
{
    for (int j = 0; j < cf->dim; j++) {
        t->bmin[j] = cf->hmin[j];
        t->bmax[j] = cf->hmax[j];
    }

  /**
  * Take care! Still need to be revised, ensure to use the
  * correct parameter to update bounding box.
    for (i = 1; i < t->size; i++) {
        for (j = 0; j < cf->dim; j++) {
            if (cf->smin[t->idx[i]][j] < t->bmin[j]) {
                t->bmin[j] = cf->smin[t->idx[i]][j];
            }
            if (cf->smax[t->idx[i]][j] > t->bmax[j]) {
                t->bmax[j] = cf->smax[t->idx[i]][j];
            }
        }
    }
  **/
}

void
update_bbox_for_cluster(pcluster t)
{
    for (int j = 0; j < t->dim; j++) {
        t->bmin[j] = t->son[0]->bmin[j];
        t->bmax[j] = t->son[0]->bmax[j];
    }
    for (int i = 1; i < t->nsons; i++) {
        for (int j = 0; j < t->dim; j++) {
            t->bmin[j] = min<double>(t->bmin[j], t->son[i]->bmin[j]);
            t->bmax[j] = max<double>(t->bmax[j], t->son[i]->bmax[j]);
        }
    }
}
