#ifndef _GEOM_H
#define _GEOM_H

#include <math.h>

void centerTrigon(double *a, double *b, double *c, double *center);
double lenEdge(double *p1, double *p2, int n);
double areaFace(double *a, double *b, double *c);
#endif //_GEOM_H
