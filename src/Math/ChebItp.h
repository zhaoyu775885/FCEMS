#ifndef CHEBITP_H
#define CHEBITP_H

#include "../Basics/Constants.h"

class ChebItp1D {
public:
    ChebItp1D(const int number) : n(number) {};
    double
    get_cheb_point(const int i, const double a, const double b)
    {
        double ang = (2.0*i+1.0)/(2.0*n+2.0)*PI;
        return cos(ang)*(b-a)*0.5 + (a+b)*0.5;
    }

private:
    int n;
};

class ChebItp3D {
public:
    ChebItp3D(const int n1, const int n2, const int n3) :
    cheb_x(n1), cheb_y(n2), cheb_z(n3) {}
    double
    get_cheb_point_x(const int ix, const double a, const double b)
    { return cheb_x.get_cheb_point(ix, a, b);}
    double
    get_cheb_point_y(const int iy, const double a, const double b)
    { return cheb_y.get_cheb_point(iy, a, b);}
    double
    get_cheb_point_z(const int iz, const double a, const double b)
    { return cheb_z.get_cheb_point(iz, a, b);}
private:
    ChebItp1D cheb_x;
    ChebItp1D cheb_y;
    ChebItp1D cheb_z;
};

#endif // CHEBITP_H
