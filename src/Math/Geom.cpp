#include "Geom.h"

void centerTrigon(double *a, double *b, double *c, double *center)
{
    center[0] = (a[0] + b[0] + c[0])/3;
    center[1] = (a[1] + b[1] + c[1])/3;
    center[2] = (a[2] + b[2] + c[2])/3;
}


double lenEdge(double *p1, double *p2, int n)
{
	int i=0;
	double temp, sum = 0;
	while(i<n) {
		temp = p1[i] - p2[i];
		sum += temp*temp;
		i++;
	}
	return sqrt(sum);
}


double areaFace(double *a, double *b, double *c)
{
    double ab, bc, ca, temp;
    ab = lenEdge(a, b, 3);
    bc = lenEdge(b, c, 3);
    ca = lenEdge(c, a, 3);
    temp = (ab + bc + ca)/2.0;
    return sqrt((temp-ab)*(temp-bc)*(temp-ca)*temp);
}
