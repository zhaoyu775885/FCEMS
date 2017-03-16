#ifndef ACA_H
#define ACA_H

#include <math.h>
#include <mkl.h>
#include "hmatrix.h"
#include "../Intg.h"

prkmatrix
aca(DComplex *Z, int m, int n);

void
mom_aca(const BasisFunc &basis_space, phmatrix hm, double wavenum, Intg intg);

void
mom_aca_plus(const BasisFunc &basis_space, phmatrix hm, double wavenum, Intg intg);

int getRowACA(DComplex *R, int m, int n, int curColNum, int *flagRow);
int getColACA(DComplex *R, int m, int n, int curRowNum, int *flagCol);
int getRefCol(int flagCol[], const int n, const int k);
int getRefRow(DComplex val[], int flagRow[], const int m, const int k);
int getRowACA_plus(DComplex *R, int m, int n, int curColNum, int *flagRow);
int
get_rowidx_acaplus(DComplex *col_val, int m, int *flagRow);
int
get_colidx_acaplus(DComplex *row_val, int n, int *flagCol);
#endif // _ACA_H
