#ifndef SINGULARINTEGRAL_H
#define SINGULARINTEGRAL_H
#include "Point3d.h"
#include "Vector.h"
#include "Geom.h"

void
Isp3_and_Isp(double *Isp3, double *Isp,
            const Point3D &pr1, const Point3D &pr2, const Point3D &pr3, const Point3D &pr0);
void
Isp_and_Is(double *Isp, double *Is,
            const Point3D &pr1, const Point3D &pr2, const Point3D &pr3, const Point3D &pr0);

void numProduct(double *r, double num, int n);

#endif // _SINGULARINTEGRAL_H
