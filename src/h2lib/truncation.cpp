#include "truncation.h"

ptruncmode
new_truncmode()
{
	ptruncmode tm = new truncmode;

	tm->frobenius = true;
	tm->absolute = false;
	tm->blocks = false;
	tm->zeta_level = 1.0;
	tm->zeta_age = 1.0;

	return tm;
}

void
del_truncmode(ptruncmode tm)
{
	delete tm;
}

ptruncmode
new_releucl_truncmode()
{
	ptruncmode tm = new_truncmode();

	return tm;
}

ptruncmode
new_relfrob_truncmode()
{
	ptruncmode tm = new_truncmode();

	tm->frobenius = true;

	return tm;
}

ptruncmode
new_blockreleucl_truncmode()
{
	ptruncmode tm = new_truncmode();
	
	tm->blocks = true;
	tm->absolute = true;
	tm->zeta_age = 1.5;

	return tm;
}

ptruncmode
new_abseucl_truncmode()
{
	ptruncmode tm = new_truncmode();

	tm->absolute = true;
	tm->zeta_level = 1.2;
	tm->zeta_age = 1.5;

	return tm;
}

int
findrank_truncmode(pctruncmode tm, double eps, pcrealavector sigma)
{
	double norm, sum;
	int k;

	if (tm && tm->frobenius) {
		if (tm->absolute) {
			sum = 0.0;
			k = sigma->dim;
			while (k>0 && (sum += fabs(sigma->v[k-1])) <= eps*eps )
				--k;
		}
		else {
			norm = 0.0;
			for (k=0;k<sigma->dim;k++)
				norm += fabs(sigma->v[k]);
			sum = 0.0;
			k = sigma->dim;
			while ( k>0 && (sum += fabs(sigma->v[k-1])) <=eps*eps*norm )
				--k;
		}
	}
	else {
		if (tm && tm->absolute) {
			k = 0;
			while (k<sigma->dim && sigma->v[k]>eps)
				++k;
		}
		else {
			k = 0;
			while (k<sigma->dim && sigma->v[k]>eps*sigma->v[0])
				++k;
		}
	}

	return k;
}
