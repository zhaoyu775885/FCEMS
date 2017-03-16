#ifndef TRUNCATION_H
#define TRUNCATION_H

typedef struct _truncmode truncmode;

typedef truncmode *ptruncmode;

typedef const truncmode *pctruncmode;

#include "h2settings.h"
#include "realavector.h"

struct _truncmode {
	bool frobenius;
	bool absolute;
	bool blocks;

	double zeta_level;
	double zeta_age;
};

ptruncmode
new_truncmode();

void
del_truncmode(ptruncmode tm);

ptruncmode
new_releucl_truncmode();

ptruncmode
new_relfrob_truncmode();

ptruncmode
new_blockreleucl_truncmode();

ptruncmode
new_blockrelfrob_truncmode();

ptruncmode
new_abseucl_truncmode();

int 
findrank_truncmode(pctruncmode tm, double eps, pcrealavector sigma);

#endif
