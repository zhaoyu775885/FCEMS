#ifndef PATCH_H
#define PATCH_H

#include <vector>
#include <string>
#include "../Basics/routine.h"
#include "../Math/Point3d.h"

typedef enum { Tri = 3, Rect = 4 } Shape;

class Patch
{
private:
    vector<int> info;
    Shape shape;
    int label;

public:
    Patch(string s, bool iris_label);
    Shape get_shape() const { return shape; }
    int operator[](int n) const {return info[n];}
    int get_label() const {return label;}

    friend ostream& operator<<(ostream& os, Patch &patch);
};

#endif // _PATCH_H
