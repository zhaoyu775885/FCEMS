#ifndef HMATRIX_H
#define HMATRIX_H

extern int hmatrix_full_mem_size;
extern int hmatrix_rk_mem_size;

typedef struct _hmatrix hmatrix;

typedef hmatrix* phmatrix;

typedef const hmatrix* pchmatrix;

#include "h2settings.h"
#include "avector.h"
#include "amatrix.h"
#include "rkmatrix.h"
#include "cluster.h"
#include "block.h"

struct _hmatrix {
    pccluster rc;          //
    pccluster cc;          //

    prkmatrix r;            //
    pamatrix f;             //

    phmatrix *son;          //

    int rsons;               //
    int csons;               //

	int refs;
	int desc;
};

void
print_hmatrix(phmatrix hm);

phmatrix
init_hmatrix(phmatrix hm, pccluster rc, pccluster cc);

void
uninit_hmatrix(phmatrix hm);

void
del_hmatrix(phmatrix hm);

void
ref_hmatrix(phmatrix *ptr, phmatrix hm);

void
unref_hmatrix(phmatrix hm);

phmatrix
new_hmatrix(pccluster rc, pccluster cc);

phmatrix
new_rk_hmatrix(pccluster rc, pccluster cc, int k);

phmatrix
new_full_hmatrix(pccluster rc, pccluster cc);

phmatrix
new_super_hmatrix(pccluster rc, pccluster cc, int rsons, int csons);

phmatrix
new_main_hmatrix(pcluster rclu, pcluster cclu, int rblk, int cblk);

int
nearRegion(int i, int j);

phmatrix
build_from_block_hmatrix(pcblock b, int k);

bool
father_fmatrix(pchmatrix hm);

void
merge_fmatrix(phmatrix hm);

bool
father_rkmatrix(pchmatrix hm);

void
merge_rkmatrix(phmatrix hm);

void
//hmat_avec_prod(phmatrix hm, char trans, pcavector x, pavector y);
fastaddeval_hmatrix_avector(field alpha, phmatrix hm, char trans, pcavector x, field beta, pavector y);

void
mvm_hmatrix_avector(phmatrix hm, char trans, pcavector x, pavector y);

void
//mvm_hmatrix(phmatrix hm, char trans, pcavector x, pavector y);
addeval_hmatrix_avector(field alpha, phmatrix hm, char trans, pcavector x, pavector y);

int
hmat_cg(phmatrix hm, pavector x, pavector y);

void
update_hmatrix(phmatrix hm);

void
clear_hmatrix(phmatrix hm);

phmatrix
clone_hmatrix(pchmatrix src);

phmatrix
clonestructure_hmatrix(pchmatrix src);

unsigned int
count_hmatrix_mem(pchmatrix hm);
#endif // HMATRIX_H
