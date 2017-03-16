#ifndef BLOCK_H
#define BLOCK_H

typedef struct _block block;

typedef block* pblock;

typedef const block* pcblock;

#include <algorithm>
#include <assert.h>
#include "cluster.h"
#include "h2settings.h"

using namespace std;

struct _block {
	pcluster rc;
	pccluster cc;

	bool a;

	pblock* son;

	int rsons;
	int csons;

	int desc;
};

typedef bool (*admissible)(pcluster rc, pcluster cc, double *data);

extern int depth;

pblock
new_block(pcluster rc, pcluster cc, bool a, int rsons, int csons);

bool
admissible_2_cluster(pcluster rc, pcluster cc, double* data);

bool
admissible_sphere_cluster(pcluster rc, pcluster cc, double *data);

pblock
build_block(pcluster rc, pcluster cc, double *eta, admissible admis);

void
init_block_cluster_csp(pblock b);

void
count_csp(pblock b);

void
traverse_cluster(pccluster clu, int &csp);

int &get_csp(pblock b, int &csp);
#endif
