#ifndef LINALG_H
#define LINALG_H

#include "settings.h"

void leqSolve_z(DComplex *a, DComplex *b, int n, int nb);
void matInv_z(DComplex *a, int n);
void mmp_z(DComplex *a, CBLAS_TRANSPOSE atrans, DComplex *b, CBLAS_TRANSPOSE btrans, DComplex *c, int m, int k, int n);

#endif // _LINALG_H
