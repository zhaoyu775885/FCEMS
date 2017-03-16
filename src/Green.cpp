#include "Green.h"

DComplex
Green0(const Point3D &r1, const Point3D &r2, const DComplex k)
{
    // Green Function in Free Space
    double distance = (r1-r2).length();
    return ( (exp(-CJ*k*distance))/(4*PI*distance) );
}

Complex3D
nablaGreen0(const Point3D &r1, const Point3D &r2, const DComplex k)
{
// Green Function in Free Space
    Point3D r0 = r1 - r2;
    double R = r0.length();
    DComplex coef = -(1.0+CJ*k*R)/(R*R*R)*exp(-CJ*k*R)/(4*PI);
    return coef * r0;
}
