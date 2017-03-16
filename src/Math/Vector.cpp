#include "Vector.h"

void numProduct(double *r, double num, int n)
{
    for (int i=0;i<n;i++) {
        r[i] *= num;
    }
}

void vectorSub(double *r1, double *r2, double *r, int n)
{
    int i;
    for (i=0;i<n;i++) {
        r[i] = r1[i] - r2[i];
    }
}

void vectorSub_cc(DComplex *r1, DComplex *r2, DComplex *r, int n)
{
    int i;
    for (i=0;i<n;i++) {
        r[i] = r1[i] - r2[i];
    }
}


double dotProduct(double *r1, double *r2, int n)
{
    int i = 0;
    double sum = 0;
    while (i<n) {
        sum += r1[i]*r2[i];
        i++;
    }
    return sum;
}

DComplex dotProduct_cd(double *r1, DComplex *r2, int n)
{
    int i = 0;
    DComplex sum = 0;
    while(i<n) {
        sum += r1[i]*r2[i];
        i++;
    }
    return sum;
}

DComplex dotProduct_cc(DComplex *z1, DComplex *z2, int n)
{
	int i=0;
	DComplex z=0;
	while (i<n) {
		z += z1[i] * z2[i];
		i++;
	}
	return z;
}

void crossProduct(double *r1, double *r2, double *r)
{
    r[0] = r1[1]*r2[2] - r1[2]*r2[1];
    r[1] = r1[2]*r2[0] - r1[0]*r2[2];
    r[2] = r1[0]*r2[1] - r1[1]*r2[0];
}

void crossProduct_cd(double *r1, DComplex *r2, DComplex *r)
{
    r[0] = r1[1]*r2[2] - r1[2]*r2[1];
    r[1] = r1[2]*r2[0] - r1[0]*r2[2];
    r[2] = r1[0]*r2[1] - r1[1]*r2[0];
}

double norm(double *r, int n)
{
    return(sqrt(dotProduct(r, r, n)));
}

void normalizeVector(double *r, int n)
{
    int i;
    double len = norm(r, n);
    for (i=0;i<n;i++) {
        r[i] /= len;
    }
}





