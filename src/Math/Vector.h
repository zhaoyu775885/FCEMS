#ifndef VECTOR_H
#define VECTOR_H

#include "../settings.h"

void numProduct(double *r, double num, int n);

void vectorSub(double *r1, double *r2, double *r, int n);
void vectorSub_cc(DComplex *r1, DComplex *r2, DComplex *r, int n);

double dotProduct(double *r1, double *r2, int n);
DComplex dotProduct_cd(double *r1, DComplex *r2, int n);
DComplex dotProduct_cc(DComplex *z1, DComplex *z2, int n);

void crossProduct(double *r1, double *r2, double *r);
void crossProduct_cd(double *r1, DComplex *r2, DComplex *r);

double norm(double *r, int n);

void normalizeVector(double *r, int n);

#endif
