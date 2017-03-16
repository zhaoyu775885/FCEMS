#ifndef CLUSTER_H
#define CLUSTER_H

typedef struct _cluster cluster;

typedef cluster* pcluster;

typedef const cluster* pccluster;

#include "clustergeometry.h"
#include "h2settings.h"


struct _cluster {
	int dim;

	int size;
	int *idx;

	double* bmin;
	double* bmax;

	int nsons;
	pcluster* son;

	int desc;
	int csp_count;
};

pcluster
new_cluster(int size, int *idx, int nsons, int dim);

void
del_cluster(pcluster p);

pcluster
copy_cluster(pcluster p);

typedef enum {
    H2_ADAPTIVE,
    H2_REGULAR,
    H2_SIMSUB,
    H2_PCA
} clustermode;

pcluster
build_adaptive_cluster(pclustergeometry cf, int size, int* idx, int clf);

pcluster
build_regular_cluster(pclustergeometry cf, int size, int* idx, int clf, int direction);

pcluster
build_simsub_cluster(pclustergeometry cf, int size, int *idx, int clf);

pcluster
build_pca_cluster(pclustergeometry cf, int size, int *idx, int clf);

pcluster
build_cluster(pclustergeometry cf, int size, int *idx, int clf, clustermode mode);


/* ------------------------------------------------------------
 * Statistics
 * ------------------------------------------------------------ */
void
print_cluster_dim(pccluster t);

int
getdepth_cluster(pccluster t);

int
getmindepth_cluster(pccluster t);

/* ------------------------------------------------------------
 * Structure Adaptation
 * ------------------------------------------------------------ */
void
extend_cluster(pcluster t, uint depth);

void
cut_cluster(pcluster t, uint depth);

void
balance_cluster(pcluster t, uint depth);

void
coarsen_cluster(pcluster t, uint minsize);

double
getdiam_2_cluster(pccluster t);

#endif
