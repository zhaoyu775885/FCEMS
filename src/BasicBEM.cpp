#include "BasicBEM.h"

void
BasicBem::build_euqation(const double freq, const char c)
{
    if (c=='p' || c=='P') {
        nrow *= 2;
        Z = new DComplex [nrow*nrow];
        U = new DComplex [nrow*nrhs];
    }
    else {
        Z = new DComplex [nrow*nrow];
        U = new DComplex [nrow*nrhs];
    }

    fill_zmat(freq, c);
    cout << "Build Z over!" << endl;
    fill_umat(freq, c);
    cout << "Build U over!" << endl;
}

void
BasicBem::fill_zmat(const double freq, const char c)
{
    int nEdge = func_space.get_dof();

    Integration intg(func_space, freq, 4.0, 1.0);

    if (c=='e' || c=='E') {
        #pragma omp parallel for
        for (int j=0;j<nEdge;++j) {
            for (int i=0;i<nEdge;++i) {
                Z[i+nEdge*j] = intg.efie(i, j);
            }
        }
    }
    else if (c=='m' || c=='M') {
        #pragma omp parallel for
        for (int j=0;j<nEdge;++j) {
            for (int i=0;i<nEdge;++i) {
                Z[i+nEdge*j] = intg.mfie(i, j);
            }
        }
    }
    else if (c=='c' || c=='C') {
        #pragma omp parallel for
        for (int j=0;j<nEdge;++j) {
            for (int i=0;i<nEdge;++i) {
                Z[i+nEdge*j] = intg.cfie(i, j);
            }
        }
    }
    else if (c=='p' || c=='P') {
        DComplex eta1_2 = ETA0*ETA0;
        DComplex eta2_2 = eta1_2/4.0;
        #pragma omp parallel for
        for (int j=0;j<nEdge;++j) {
            for (int i=0;i<nEdge;++i) {
                DComplex l1 = intg.efie(i, j, 1);
                DComplex l2 = intg.efie(i, j, 4);
                DComplex k1 = intg.pmchwt(i, j, 1);
                DComplex k2 = intg.pmchwt(i, j, 4);

                Z[i+nrow*j] = l1 + l2;
                Z[i+nrow*(j+nEdge)] = -(k1+k2);
                Z[i+nEdge+nrow*j] = -Z[i+nrow*(j+nEdge)];
                Z[i+nEdge+nrow*(j+nEdge)] = l1/eta1_2 + l2/eta2_2;
            }
        }
    }
}

void
BasicBem::fill_umat(const double freq, const char c)
{
    int nEdge = func_space.get_dof();
    for (int i=0;i<nEdge;i++) {
        if (c=='e' || c=='E') {
            U[i] = innerProduct_UE(func_space[i], freq);
        }
        else if (c=='m' || c=='M') {
            U[i] = innerProduct_UH(func_space[i], freq);
        }
        else if (c=='c' || c=='C') {
            U[i] = innerProduct_UE(func_space[i], freq) +
                ETA0 * innerProduct_UH(func_space[i], freq);
        }
        else if (c=='p' || c=='P') {
            U[i] = innerProduct_UE(func_space[i], freq);
            U[i+nEdge] = innerProduct_UTH(func_space[i], freq);
        }
    }
}

void
BasicBem::solve()
{
    leqSolve_z(Z, U, nrow, nrhs);
}

void
BasicBem::write_umat(const string &filename)
{
    ofstream fp(filename.c_str());
    for (int i=0;i<nrow;i++) {
        for (int j=0;j<nrhs;j++) {
            fp << U[i+j*nrow] << " ";
        }
        fp << endl;
    }
    fp.close();
}

void
BasicBem::write_zmat(const string &filename)
{
    ofstream fp(filename.c_str());
    for (int i=0;i<nrow;i++) {
        for (int j=0;j<1;j++) {
            fp << Z[i+nrow*j] << endl;
        }
    }
    fp.close();
}
