#include "Point3d.h"

Complex3D
operator*(const Point3D &d, const DComplex &coef)
{
    return Complex3D(d[0]*coef, d[1]*coef, d[2]*coef);
}

Complex3D
operator*(const DComplex &coef, const Point3D &d)
{
    return Complex3D(d[0]*coef, d[1]*coef, d[2]*coef);
}

Complex3D
operator*(const Complex3D &c, const double &coef)
{
    return Complex3D(c[0]*coef, c[1]*coef, c[2]*coef);
}

Complex3D
operator*(const double &coef, const Complex3D &c)
{
    return Complex3D(c[0]*coef, c[1]*coef, c[2]*coef);
}

DComplex
dotprod(const Point3D &d, const Complex3D &c)
{
    return d[0]*c[0]+d[1]*c[1]+d[2]*c[2];
}

DComplex
dotprod(const Complex3D &c, const Point3D &d)
{
    return d[0]*c[0]+d[1]*c[1]+d[2]*c[2];
}

//DComplex
//dotprod(const Complex3D &c, const Complex3D &d)
//{
//    return d[0]*c[0]+d[1]*c[1]+d[2]*c[2];
//}

Complex3D
crossprod(const Point3D &p1, const Complex3D &p2)
{
    Complex3D newp(p1[1]*p2[2]-p1[2]*p2[1], p1[2]*p2[0]-p1[0]*p2[2], p1[0]*p2[1]-p1[1]*p2[0]);
    return newp;
}

Complex3D
crossprod(const Complex3D &p1, const Point3D &p2)
{
    Complex3D newp(p1[1]*p2[2]-p1[2]*p2[1], p1[2]*p2[0]-p1[0]*p2[2], p1[0]*p2[1]-p1[1]*p2[0]);
    return newp;
}

ostream& operator<< (ostream &os, const Point3D &point)
{
    for (int i=0;i<3;i++) os << scientific << point[i] << " ";
    os << endl;
    return os;
}

Point3D
get_trigon_ctr(const Point3D &p0, const Point3D &p1, const Point3D &p2)
{
    Point3D ctr;
    for (int i=0;i<3;++i) ctr[i] = (p0[i]+p1[i]+p2[i]) / 3.0;
    return ctr;
}

Point3D
get_trigon_nvec(const Point3D &p0, const Point3D &p1, const Point3D &p2)
{
    Point3D v_01(p1-p0);
    Point3D v_12(p2-p1);
    Point3D nvec(crossprod(v_01, v_12));
    return nvec.normalize();
}

double
get_trigon_area(const Point3D &p0, const Point3D &p1, const Point3D &p2)
{
    Point3D v_01(p1-p0);
    Point3D v_12(p2-p1);
    Point3D nvec(crossprod(v_01, v_12));
    return nvec.length()*0.5;
}

