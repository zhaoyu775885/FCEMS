#ifndef INTG_H
#define INTG_H

#include "BasisFunction.h"
#include "Math/SingularIntegral.h"
#include "Green.h"

#include "Basics/Constants.h"
#include "Math/GaussQuad.h"
#include "Math/Point3d.h"

typedef DComplex (*Intg)(const FullBF&, const FullBF&, const double , const double);

DComplex Intg_L (const HalfBF * hbf1, const HalfBF * hbf2, const double freq, const DComplex er=1.0);
DComplex Intg_KN(const HalfBF * hbf1, const HalfBF * hbf2, const double freq, const DComplex er=1.0);
DComplex Intg_KT(const HalfBF * hbf1, const HalfBF * hbf2, const double freq, const DComplex er=1.0);
double Intg_hbf (const HalfBF * hbf1, const HalfBF * hbf2);
void
LO_normal(const HalfBF *hbf1, const HalfBF *hbf2, const double freq,
        const DComplex er, const Point3D &, Complex3D &l1, DComplex &l2);
void
LO_normal(const HalfBF *hbf, const double freq, const DComplex er,
        const Point3D &obs_gauss_point, Complex3D &l1, DComplex &l2);
Complex3D
LO_vec_ntest(const HalfBF *hbf, const double freq,
             const DComplex er, const Point3D &obs_gauss_point);

Complex3D
LO_vec_ttest(const HalfBF *hbf, const double freq,
             const DComplex er, const Point3D &obs_gauss_point);

Complex3D
KO_normal(const HalfBF *hbf1, const HalfBF *hbf2, const double freq,
        const DComplex er, const Point3D &obs_gauss_point);
Complex3D
KO_normal(const HalfBF *hbf, const double freq,
          const DComplex er, const Point3D &obs_gauss_point);

class Integration
{
public:
    Integration(const BasisFunc &basisfunc, const double frequency) :
        func_space(basisfunc), basis_space(func_space.get_basis_space()),
        freq(frequency), er_in(1.0), er_ex(1.0) {}
    Integration(const BasisFunc &basisfunc, const double frequency,
                const DComplex in, const DComplex ex) :
        func_space(basisfunc), basis_space(func_space.get_basis_space()),
        freq(frequency), er_in(in), er_ex(ex) {}
    DComplex efie(const int obs, const int src, const DComplex er=1.0) const ;
//    DComplex efie(const int obs, const int src, const bool out=true) const ;
    DComplex mfie(const int obs, const int src, const DComplex er=1.0) const ;
//    DComplex mfie(const int obs, const int src, const bool out=true) const ;
    DComplex cfie(const int obs, const int src, const DComplex er=1.0) const ;
    DComplex pmchwt(const int obs, const int src, const DComplex er=1.0) const;
    DComplex xfie(const int obs, const int src, const char c, const DComplex er=1.0) const;
    double get_freq() const {return freq;}

    void single_layer_intg(const int , const Point3D &, Complex3D &, DComplex &, const bool out=true) const;
    Complex3D
    single_layer_intg_m1(const int obs, const Point3D &src_point, const bool=true) const;
    Complex3D
    single_layer_intg_m2(const int obs, const Point3D &src_point, const bool=true) const;
    Complex3D
    single_layer_intg_m3(const int obs, const Point3D &src_point, const bool=true) const;
    Complex3D
    single_layer_intg_m4(const int obs, const Point3D &src_point, const bool=true) const;

private:
    const BasisFunc &func_space;
    const vector<FullBF> &basis_space;
    const double freq;
    const DComplex er_in;
    const DComplex er_ex;
};

DComplex innerProduct_UE(const FullBF &rbf, double freq);
DComplex innerProduct_UTH(const FullBF &rbf, double freq);
DComplex innerProduct_UH(const FullBF &rbf, double freq);

inline Point3D &
get_gauss_point(const HalfBF *hbf, double const *coef, Point3D &point)
{
    if (hbf->type == Tri) {
        point = coef[0]*hbf->hRwg->v0 +
                coef[1]*hbf->hRwg->v1 +
                coef[2]*hbf->hRwg->v2;
    }
    else if (hbf->type == Rect) {
        point = coef[0]*hbf->hRtp->v1 +
                coef[1]*hbf->hRtp->v2 +
                coef[2]*hbf->hRtp->v3 +
                coef[3]*hbf->hRtp->v4;
    }
    return point;
}


#endif // INTG_H
