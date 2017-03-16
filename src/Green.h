#ifndef GREEN_H
#define GREEN_H

#include "Basics/Constants.h"
#include "settings.h"
#include "Math/Point3d.h"

DComplex
Green0(const Point3D &r1, const Point3D &r2, const DComplex k);

Complex3D
nablaGreen0(const Point3D &r1, const Point3D &r2, const DComplex k);

#endif // GREEN_H
