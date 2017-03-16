#ifndef VECTOR3D_H
#define VECTOR3D_H

#include <vector>
#include <string>
#include <assert.h>
#include "../settings.h"
#include "../Basics/routine.h"

using namespace std;

template <typename T>
class Vector3D
{
public:
    Vector3D() : coord(3, 0) {}
    Vector3D(const vector<T> &vec) : coord(vec) { assert(vec.size()==3); }
    Vector3D(const T &x, const T &y, const T &z) : coord(3, 0) { coord[0] = x; coord[1] = y; coord[2] = z; }
    Vector3D(const Vector3D &old) : coord(3, 0) {coord[0]=old[0]; coord[1]=old[1]; coord[2]=old[2];}
    Vector3D(const string &s, bool iris_mesh_line);

    T &operator[] (int i) {return coord[i];}
    const T &operator[] (int i) const {return coord[i];}

    Vector3D& operator+=(const Vector3D &p)
    {for (typename vector<T>::size_type i=0;i!=coord.size();++i) coord[i] += p[i];return *this;}
    Vector3D& operator-=(const Vector3D &p)
    {for (typename vector<T>::size_type i=0;i!=coord.size();++i) coord[i] -= p[i];return *this;}
    Vector3D& operator*=(const T &coef)
    {for (typename vector<T>::size_type i=0;i!=coord.size();++i) coord[i] *= coef;return *this;}
    Vector3D& operator/=(const T &coef)
    {assert(abs(coef)!=0);for (typename vector<T>::size_type i=0;i!=coord.size();++i) coord[i] /= coef;return *this;}


    double length() const { return sqrt(abs(coord[0]*coord[0]+coord[1]*coord[1]+coord[2]*coord[2])); }
    Vector3D &normalize() { *this /= this->length(); return *this; }

    template <typename U> friend ostream& operator<<(ostream &os, const Vector3D<U> &point);
//    friend ostream& operator<< (ostream &os, const Vector3D &point);

private:
    vector<T> coord;
};

template <typename T>
Vector3D<T>::Vector3D(const string &node_line, bool iris_mesh_line)
{
    if (iris_mesh_line) {
        vector<string> str_list = StringSplit(node_line, ' ');
        for (size_t k=1;k<str_list.size();k++) {
            coord.push_back(1e-6 * Str2Num<double>(str_list[k]));
        }
        assert(coord.size() == 3);
    }
    else {
        vector<string> str_list = StringSplit(node_line, ' ');
        for (size_t k=0;k<3;k++) {
            coord.push_back(Str2Num<double>(str_list[k]));
        }
        assert(coord.size() == 3);
    }
}

template <typename T>
inline Vector3D<T> operator+(const Vector3D<T> &p1, const Vector3D<T> &p2)
{
    Vector3D<T> p0(p1);
    p0+=p2;
    return p0;
}

template <typename T>
inline Vector3D<T> operator-(const Vector3D<T> &p1, const Vector3D<T> &p2)
{
    Vector3D<T> p0(p1);
    p0-=p2;
    return p0;
}

template <typename T>
inline Vector3D<T> operator*(const Vector3D<T> &p1, const T &multiplier)
{
    Vector3D<T> newp(p1);
    newp *= multiplier;
    return newp;
}

template <typename T>
inline Vector3D<T> operator*(const T &multiplier, const Vector3D<T> &p1)
{
    Vector3D<T> newp(p1);
    newp *= multiplier;
    return newp;
}

template <typename T>
inline Vector3D<T> operator/(const Vector3D<T> &p1, const T &divisor)
{
    Vector3D<T> newp(p1);
    newp /= divisor;
    return newp;
}

template <typename T>
inline Vector3D<T> operator*(const Vector3D<T> &p1, const Vector3D<T> &p2)
{
    Vector3D<T> newp(p1[1]*p2[2]-p1[2]*p2[1], p1[2]*p2[0]-p1[0]*p2[2], p1[0]*p2[1]-p1[1]*p2[0]);
    return newp;
}

template <typename T>
inline T dotprod(const Vector3D<T> &p1, const Vector3D<T> &p2)
{
    T sum = 0;
    for (int i=0;i<3;i++) sum += p1[i] * p2[i];
    return sum;
}

template <typename T>
inline Vector3D<T> crossprod(const Vector3D<T> &p1, const Vector3D<T> &p2)
{
    return p1*p2;
}

template <typename U>
ostream& operator<<(ostream &os, const Vector3D<U> &point)
{
    for (int i=0;i<3;i++) os << scientific << point[i] << " ";
    os << endl;
    return os;
}

#endif // _POINT3D_H

