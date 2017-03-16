#ifndef BASISFUNC_H
#define BASISFUNC_H

#include "Mesh/Mesh.h"

typedef enum {Half=0, Rwg = 1, Rtp = 2, Mix = 3} BFType;

typedef struct {
    Point3D v0, v1, v2;
    Point3D ctr;
    Point3D nvec;
    double len;
    double area;
    int face_label;
} halfrwgBF;

typedef struct{
    Point3D v1, v2, v3, v4;
	Point3D vct;
	Point3D nvec;
	double lx;  //lx 是基函数长度
	double ly;  //ly 是公共边长度
	int face_label;
} halfrtpBF;

struct HalfBF
{
    Shape type;
    int face_label;
    halfrwgBF *hRwg;
    halfrtpBF *hRtp;
};

struct FullBF
{
    BFType type;
    HalfBF *pHalf;
    HalfBF *nHalf;
};

class BasisFunc
{
public:
    BasisFunc(const Mesh &mesh);

    const FullBF &operator[](const int i) const {return bf_database[i];}

    const Mesh & get_mesh() const {return g_mesh;}

    int get_common_edge_num() const {return common_edge_num;}
    int get_dof() const {return bf_database.size();}

    const vector<FullBF> &get_basis_space() const {return bf_database;}
    const FullBF& get_basis(const int i) const {return bf_database[i];}

    FullBF buildBF(const Patch &face_a, int index_a, const Patch &face_b, int index_b);
    FullBF buildBF(const Patch &face_a, int index_a);
    HalfBF buildPortBF(const Patch &face, int index);
    HalfBF * buildHalfBF(const Patch &face, int index);
    halfrwgBF * buildHalfRwg(const Patch &face, int index) const;
    halfrtpBF * buildHalfRtp(const Patch &face, int index) const;
//    int judge_port(const Patch &face, int &index);

    Point3D calPos_hrwg(halfrwgBF *hrbf, double const *coef);

private:
    const Mesh &g_mesh;
    vector<FullBF> bf_database;
    int common_edge_num;
};
#endif // _BASEFUNCTION_H
