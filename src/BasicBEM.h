#ifndef BASICBEM_H
#define BASICBEM_H

#include "BasisFunction.h"
#include "Intg.h"
#include "LinAlg.h"

class BasicBem
{
public:
    BasicBem(const BasisFunc &bf) : func_space(bf),Z(0), U(0),
                    nrow(func_space.get_dof()), nrhs(1) {}

    void build_euqation(const double freq, const char c='e');
    void fill_zmat(const double freq, const char c);
    void fill_umat(const double freq, const char c);
    void solve();
    const DComplex *export_data() const {return U;}

    void write_zmat(const string &filename);
    void write_umat(const string &filename);
    const int get_nrow() const {return nrow;}
    const int get_nrhs() const {return nrhs;}
    const BasisFunc &get_func_space() const {return func_space;}

private:
    const BasisFunc &func_space;
    DComplex *Z;
    DComplex *U;
    int nrow;
    int nrhs;
};

#endif // FILLMAT_H
