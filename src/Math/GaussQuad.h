#ifndef QUADRATURE_H
#define QUADRATURE_H

#include "../settings.h"

#define MAX_GAUSSNUM 25

#if TST
#define DEFAULT_TRI_GAUSS_RULE 0
#else
#define DEFAULT_TRI_GAUSS_RULE 1
#endif // TST

int const triRules = 10;
extern int const tri_gauss_num[triRules];
extern double const tri_gauss_weight[triRules][25];
extern double const tri_gauss_coef[triRules][25][3];

int const rectRules = 10;
extern int const rect_gauss_num[rectRules];
extern double const rect_gauss_weight[rectRules][25];
extern double const rect_gauss_coef[rectRules][25][4];

#endif // _QUADRATURE_H
