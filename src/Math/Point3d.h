#ifndef POINT3D_H
#define POINT3D_H

#include <vector>
#include <string>
#include <assert.h>
#include "../settings.h"
#include "../Basics/routine.h"
#include "Vector3D.h"

using namespace std;



#if 0
class Point3D
{
public:
    Point3D() : coord(3, 0) {}
    Point3D(const vector<double> &v) : coord(v) { assert(v.size()==3); }
    Point3D(const double &x, const double &y, const double &z) : coord(3, 0) { coord[0] = x; coord[1] = y; coord[2] = z; }
    Point3D(const string &s);

    double &operator[] (int n) {return coord[n];}
    const double &operator[] (int n) const {return coord[n];}

    Point3D& operator+=(const Point3D &p)
    { coord[0] += p[0]; coord[1] += p[1]; coord[2] += p[2]; return *this; }
    Point3D& operator-=(const Point3D &p)
    { coord[0] -= p[0]; coord[1] -= p[1]; coord[2] -= p[2]; return *this; }
    Point3D& operator*=(const double &coef)
    { coord[0] *= coef; coord[1] *= coef; coord[2] *= coef; return *this; }
    Point3D& operator/=(const double &coef)
    { coord[0] /= coef; coord[1] /= coef; coord[2] /= coef; return *this; }

    Point3D &normalize() { *this /= this->length(); return *this; }
    double length() const
    {
        double len_2 = 0;
        for (int i=0;i<3;++i) len_2 += coord[i]*coord[i];
        return sqrt(len_2);
    }

    friend ostream& operator<< (ostream &os, const Point3D &point);

private:
    vector<double> coord;
};

inline Point3D
operator+(const Point3D &p1, const Point3D &p2)
{
    Point3D point(p1[0]+p2[0], p1[1]+p2[1], p1[2]+p2[2]);
    return point;
}

inline Point3D
operator-(const Point3D &p1, const Point3D &p2)
{
    Point3D point(p1[0]-p2[0], p1[1]-p2[1], p1[2]-p2[2]);
    return point;
}

inline Point3D
operator*(const Point3D &p1, const double &multiplier)
{
    Point3D newp(p1[0]*multiplier, p1[1]*multiplier, p1[2]*multiplier);
    return newp;
}

inline Point3D
operator*(const double &multiplier, const Point3D &p1)
{
    Point3D newp(p1[0]*multiplier, p1[1]*multiplier, p1[2]*multiplier);
    return newp;
}

inline Point3D
operator/(const Point3D &p1, const double &divisor)
{
    Point3D newp(p1[0]/divisor, p1[1]/divisor, p1[2]/divisor);
    return newp;
}

inline Point3D
operator*(const Point3D &p1, const Point3D &p2)
{
    Point3D newp(p1[1]*p2[2]-p1[2]*p2[1], p1[2]*p2[0]-p1[0]*p2[2], p1[0]*p2[1]-p1[1]*p2[0]);
    return newp;
}

inline double
dotprod(const Point3D &p1, const Point3D &p2)
{
    double sum = 0;
    for (int i=0;i<3;i++) sum += p1[i] * p2[i];
    return sum;
}

inline Point3D
crossprod(const Point3D &p1, const Point3D &p2)
{
    return p1*p2;
}
#else

typedef Vector3D<double> Point3D;
typedef Vector3D<DComplex> Complex3D;

Complex3D
operator*(const Point3D &d, const DComplex &coef);

Complex3D
operator*(const DComplex &coef, const Point3D &d);

Complex3D
operator*(const Complex3D &d, const double &coef);

Complex3D
operator*(const double &coef, const Complex3D &d);

DComplex
dotprod(const Point3D &d, const Complex3D &c);

DComplex
dotprod(const Complex3D &c, const Point3D &d);

//DComplex
//dotprod(const Complex3D &c, const Complex3D &d);

Complex3D
crossprod(const Point3D &d, const Complex3D &c);

Complex3D
crossprod(const Complex3D &c, const Point3D &d);

Point3D
get_trigon_ctr(const Point3D &p0, const Point3D &p1, const Point3D &p2);

Point3D
get_trigon_nvec(const Point3D &p0, const Point3D &p1, const Point3D &p2);

double
get_trigon_area(const Point3D &p0, const Point3D &p1, const Point3D &p2);

#endif // 0
#endif // _POINT3D_H
