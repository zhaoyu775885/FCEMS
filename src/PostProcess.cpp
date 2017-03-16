#include "PostProcess.h"

void
PostProcess::gen_rcs(const double beg, const double end, const int count)
{
    assert(beg>=0);
    assert(end<=360);
    assert(beg<end);
    assert(count>0);

    const DComplex *ecurrent, *mcurrent;
    int e_size, m_size;

    if (type=='p' || type=='P') {
        assert (rows%2==0);
        e_size = rows / 2;
        m_size = rows / 2;
        ecurrent = data;
        mcurrent = data + e_size;
    }
    else {
        ecurrent = data;
        mcurrent = 0;
        e_size = rows;
        m_size = 0;
    }

    const int rcs_size = count+1;
    const double step = (end-beg) / count;
    deg_vec.resize(rcs_size);
    rcs_vec.resize(rcs_size);

    for (int i=0;i<rcs_size;i++) deg_vec[i] = (beg + i*step) * (PI/180.0);

    double omega = 2*PI*freq;
    double wavenum = omega/SPEED_OF_LIGHT;

    for(int i=0;i<rcs_size;i++) {
        DComplex eyes[3*3] = {1, 0, 0, 0, 1, 0, 0, 0, 1};
        DComplex ar[3];
        DComplex arar[3*3];
        DComplex templ[3] = {0};
        DComplex tempk[3] = {0};
        Complex3D tmpl(0, 0, 0);
        Complex3D tmpk(0, 0, 0);
        DComplex es_by_r[3];
        ar[0] = sin(deg_vec[i]);
        ar[1] = 0;
        ar[2] = cos(deg_vec[i]);

        mmp_z(ar, CblasNoTrans, ar, CblasNoTrans, arar, 3, 1, 3);

        Complex3D ar3D(ar[0], ar[1], ar[2]);

        for(int j=0;j<e_size;j++) {
            tmpl += (-CJ*omega*MU0)/(8*PI) * ecurrent[j]*basis_space[j].pHalf->hRwg->len* (
                    (basis_space[j].pHalf->hRwg->ctr-basis_space[j].pHalf->hRwg->v0)*exp(CJ*wavenum*dotprod(basis_space[j].pHalf->hRwg->ctr, ar3D)) +
                    (basis_space[j].nHalf->hRwg->v0-basis_space[j].nHalf->hRwg->ctr)*exp(CJ*wavenum*dotprod(basis_space[j].nHalf->hRwg->ctr, ar3D)) );
        }

        for (int j=0;j<m_size;++j) {
            Complex3D tempP = crossprod(basis_space[j].pHalf->hRwg->ctr-basis_space[j].pHalf->hRwg->v0, ar3D);
            Complex3D tempN = crossprod(basis_space[j].nHalf->hRwg->ctr-basis_space[j].nHalf->hRwg->v0, ar3D);

            tmpk += (-CJ*wavenum)/(8*PI)*basis_space[j].pHalf->hRwg->len*mcurrent[j] * (
                    tempP*exp(CJ*wavenum*dotprod(basis_space[j].pHalf->hRwg->ctr, ar3D)) -
                    tempN*exp(CJ*wavenum*dotprod(basis_space[j].nHalf->hRwg->ctr, ar3D))
                    );
        }

        templ[0] = tmpl[0];templ[1] = tmpl[1];templ[2] = tmpl[2];
        tempk[0] = tmpk[0];tempk[1] = tmpk[1];tempk[2] = tmpk[2];
        DComplex alpha1(-1), alpha2(1);

        cblas_zaxpy(9, &alpha1, arar, 1, eyes, 1);

        mmp_z(eyes, CblasNoTrans, templ, CblasNoTrans, es_by_r, 3, 3, 1);

        cblas_zaxpy(3, &alpha2, tempk, 1, es_by_r, 1);

        double tmp = rcs_vec[i];

        cblas_zdotc_sub(3, es_by_r, 1, es_by_r, 1, &tmp);

        rcs_vec[i] = 10*log10(tmp*4*PI/(Em*Em));
    }
}

void
PostProcess::write_rcs(const string &filename)
{
    ofstream fp(filename.c_str());
    for (size_t i=0;i<deg_vec.size();i++) {
        fp << deg_vec[i]*180/PI << " " << rcs_vec[i] << endl;
    }
    cout << deg_vec[0]*180/PI << " " << rcs_vec[0] << endl;
    fp.close();
}

