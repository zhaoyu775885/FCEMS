#include "Intg.h"

DComplex
Integration::xfie(const int obs, const int src, const char c, const DComplex er) const
{
    if (c=='c' || c=='C') {
        return this->cfie(obs, src, er);
    }
    else if (c=='e' || c=='E') {
        return this->efie(obs, src, er);
    }
    else if (c=='m' || c=='M') {
        return this->mfie(obs, src, er);
    }
    else if (c=='p' || c=='P') {
        return this->pmchwt(obs, src);
    }
    else if (c=='w' || c=='W') {
        return this->efie(obs, src, er_in)+this->efie(obs, src, er_ex);
    }
    else if (c=='x' || c=='X') {
        return -this->pmchwt(obs, src, er_in)-this->pmchwt(obs, src, er_ex);
    }
    else if (c=='y' || c=='Y') {
        return this->pmchwt(obs, src, er_in)+this->pmchwt(obs, src, er_ex);
    }
    else if (c=='z' || c=='Z') {
        DComplex eta0_2 = ETA0*ETA0;
        return (this->efie(obs, src, er_in)*er_in+this->efie(obs, src, er_ex)*er_ex)/eta0_2;
    }
    return 0;
}

//DComplex
//Integration::efie(const int obs, const int src, const bool out) const
//{
//    if (out) {
//        return this->efie(obs, src, er_ex);
//    }
//    else {
//        return this->efie(obs, src, er_in);
//    }
//}

DComplex
Integration::efie(const int obs, const int src, const DComplex er) const
{
	DComplex iP_pp(0), iP_pn(0), iP_np(0), iP_nn(0);
	DComplex result;

    const FullBF &basis_obs = basis_space[obs];
    const FullBF &basis_src = basis_space[src];

	if (basis_obs.pHalf) {
	    if (basis_src.pHalf) iP_pp = Intg_L(basis_obs.pHalf, basis_src.pHalf, freq, er);
	    if (basis_src.nHalf) iP_pn = Intg_L(basis_obs.pHalf, basis_src.nHalf, freq, er);
	}
	if (basis_obs.nHalf) {
	    if (basis_src.pHalf) iP_np = Intg_L(basis_obs.nHalf, basis_src.pHalf, freq, er);
	    if (basis_src.nHalf) iP_nn = Intg_L(basis_obs.nHalf, basis_src.nHalf, freq, er);
	}
	result = iP_pp - iP_pn - iP_np + iP_nn;

	return result;
}

//DComplex
//Integration::mfie(const int obs, const int src, const bool out) const
//{
//    if (out) {
//        return this->mfie(obs, src, er_ex);
//    }
//    else {
//        return this->mfie(obs, src, er_in);
//    }
//}

DComplex
Integration::mfie(const int obs, const int src, const DComplex er) const
{
	DComplex iP_pp(0), iP_pn(0), iP_np(0), iP_nn(0);
	DComplex result;

    const FullBF &basis_obs = basis_space[obs];
    const FullBF &basis_src = basis_space[src];

	if (basis_obs.pHalf) {
	    if (basis_src.pHalf) iP_pp = Intg_KN(basis_obs.pHalf, basis_src.pHalf, freq, er) +
                                0.5*Intg_hbf(basis_obs.pHalf, basis_src.pHalf);
	    if (basis_src.nHalf) iP_pn = Intg_KN(basis_obs.pHalf, basis_src.nHalf, freq, er) +
                                0.5*Intg_hbf(basis_obs.pHalf, basis_src.nHalf);
	}
	if (basis_obs.nHalf) {
	    if (basis_src.pHalf) iP_np = Intg_KN(basis_obs.nHalf, basis_src.pHalf, freq, er) +
                                0.5*Intg_hbf(basis_obs.nHalf, basis_src.pHalf);
	    if (basis_src.nHalf) iP_nn = Intg_KN(basis_obs.nHalf, basis_src.nHalf, freq, er) +
                                0.5*Intg_hbf(basis_obs.nHalf, basis_src.nHalf);
	}
	result = iP_pp - iP_pn - iP_np + iP_nn;

	return result;
}

DComplex
Integration::cfie(const int obs, const int src, const DComplex er) const
{
    DComplex efie_value = this->efie(obs, src, er);
    DComplex mfie_value = this->mfie(obs, src, er);
    DComplex eta = ETA0;

    DComplex result = efie_value + eta*mfie_value;
	return result;
}

//DComplex
//Integration::pmchwt(const int obs, const int src, const bool out) const
//{
//    if (out) {
//        pmchwt(obs, src, er_ex);
//    }
//    else {
//        pmchwt(obs, src, er_in);
//    }
//}

DComplex
Integration::pmchwt(const int obs, const int src, const DComplex er) const
{
	DComplex iP_pp(0), iP_pn(0), iP_np(0), iP_nn(0);
	DComplex result;

    const FullBF &basis_obs = basis_space[obs];
    const FullBF &basis_src = basis_space[src];

	if (basis_obs.pHalf) {
	    if (basis_src.pHalf) iP_pp = Intg_KT(basis_obs.pHalf, basis_src.pHalf, freq, er);
	    if (basis_src.nHalf) iP_pn = Intg_KT(basis_obs.pHalf, basis_src.nHalf, freq, er);
	}
	if (basis_obs.nHalf) {
	    if (basis_src.pHalf) iP_np = Intg_KT(basis_obs.nHalf, basis_src.pHalf, freq, er);
	    if (basis_src.nHalf) iP_nn = Intg_KT(basis_obs.nHalf, basis_src.nHalf, freq, er);
	}
	result = iP_pp - iP_pn - iP_np + iP_nn;

	return result;
}

void
Integration::single_layer_intg(const int obs, const Point3D &src_point,
                    Complex3D &vecT, DComplex &scaT, const bool out) const
{
    const FullBF &basis_obs = basis_space[obs];
    Complex3D vecTermPos(0, 0, 0), vecTermNeg(0, 0, 0);
    DComplex scaTermPos(0), scaTermNeg(0);

    if (basis_obs.pHalf) {
        if (out)    LO_normal(basis_obs.pHalf, freq, er_ex, src_point, vecTermPos, scaTermPos);
        else        LO_normal(basis_obs.pHalf, freq, er_in, src_point, vecTermPos, scaTermPos);
        const double coef = basis_obs.pHalf->hRwg->len;
        vecTermPos *= 0.5*coef;
        scaTermPos *= coef;
    }
    if (basis_obs.nHalf) {
        if (out)    LO_normal(basis_obs.nHalf, freq, er_ex, src_point, vecTermNeg, scaTermNeg);
        else        LO_normal(basis_obs.nHalf, freq, er_in, src_point, vecTermNeg, scaTermNeg);
        const double coef = basis_obs.nHalf->hRwg->len;
        vecTermNeg *= 0.5*coef;
        scaTermNeg *= coef;
    }

    vecT = vecTermPos - vecTermNeg;
    scaT = scaTermPos - scaTermNeg;
}

Complex3D
Integration::single_layer_intg_m1(const int obs, const Point3D &src_point, const bool out) const
{
    const FullBF &basis_obs = basis_space[obs];
    Complex3D vecTermPos(0, 0, 0), vecTermNeg(0, 0, 0);

    if (basis_obs.pHalf) {
        if (out)    vecTermPos = LO_vec_ntest(basis_obs.pHalf, freq, er_ex, src_point);
        else        vecTermPos = LO_vec_ntest(basis_obs.pHalf, freq, er_in, src_point);
        vecTermPos *= 0.5*basis_obs.pHalf->hRwg->len;
    }
    if (basis_obs.nHalf) {
        if (out)    vecTermNeg = LO_vec_ntest(basis_obs.nHalf, freq, er_ex, src_point);
        else        vecTermNeg = LO_vec_ntest(basis_obs.nHalf, freq, er_in, src_point);
        vecTermNeg *= 0.5*basis_obs.nHalf->hRwg->len;
    }

    return vecTermPos - vecTermNeg;
}

Complex3D
Integration::single_layer_intg_m2(const int obs, const Point3D &obs_point, const bool out) const
{
    const FullBF &basis_obs = basis_space[obs];
    Complex3D vecTermPos(0, 0, 0), vecTermNeg(0, 0, 0);

    if (basis_obs.pHalf) {
        if (out)    vecTermPos = KO_normal(basis_obs.pHalf, freq, er_ex, obs_point);
        else        vecTermPos = KO_normal(basis_obs.pHalf, freq, er_in, obs_point);
        vecTermPos *= 0.5*basis_obs.pHalf->hRwg->len;
    }
    if (basis_obs.nHalf) {
        if (out)    vecTermNeg = KO_normal(basis_obs.nHalf, freq, er_ex, obs_point);
        else        vecTermNeg = KO_normal(basis_obs.nHalf, freq, er_in, obs_point);
        vecTermNeg *= 0.5*basis_obs.nHalf->hRwg->len;
    }

    return vecTermPos - vecTermNeg;
}


Complex3D
Integration::single_layer_intg_m3(const int obs, const Point3D &src_point, const bool out) const
{
    const FullBF &basis_obs = basis_space[obs];
    Complex3D vecTermPos(0, 0, 0), vecTermNeg(0, 0, 0);

    if (basis_obs.pHalf) {
        if (out)    vecTermPos = LO_vec_ttest(basis_obs.pHalf, freq, er_ex, src_point);
        else        vecTermPos = LO_vec_ttest(basis_obs.pHalf, freq, er_in, src_point);
        vecTermPos *= 0.5*basis_obs.pHalf->hRwg->len;
    }
    if (basis_obs.nHalf) {
        if (out)   vecTermNeg = LO_vec_ttest(basis_obs.nHalf, freq, er_ex, src_point);
        else       vecTermNeg = LO_vec_ttest(basis_obs.nHalf, freq, er_in, src_point);
        vecTermNeg *= 0.5*basis_obs.nHalf->hRwg->len;
    }

    return vecTermPos - vecTermNeg;
}

Complex3D
Integration::single_layer_intg_m4(const int obs, const Point3D &obs_point, const bool out) const
{
    const FullBF &basis_obs = basis_space[obs];
    Complex3D vecTermPos(0, 0, 0), vecTermNeg(0, 0, 0);

    if (basis_obs.pHalf) {
        if (out)    vecTermPos = KO_normal(basis_obs.pHalf, freq, er_ex, obs_point);
        else        vecTermPos = KO_normal(basis_obs.pHalf, freq, er_in, obs_point);
        vecTermPos *= 0.5*basis_obs.pHalf->hRwg->len;
    }
    if (basis_obs.nHalf) {
        if (out)    vecTermNeg = KO_normal(basis_obs.nHalf, freq, er_ex, obs_point);
        else        vecTermNeg = KO_normal(basis_obs.nHalf, freq, er_in, obs_point);
        vecTermNeg *= 0.5*basis_obs.nHalf->hRwg->len;
    }

    return vecTermPos - vecTermNeg;
}

#if 0
DComplex efie(const FullBF &obs, const FullBF &src, const double freq, const double er)
{
	DComplex iP_pp(0), iP_pn(0), iP_np(0), iP_nn(0);
	DComplex result;

	if (obs.pHalf) {
	    if (src.pHalf) iP_pp = Intg_L(obs.pHalf, src.pHalf, freq);
	    if (src.nHalf) iP_pn = Intg_L(obs.pHalf, src.nHalf, freq);
	}
	if (obs.nHalf) {
	    if (src.pHalf) iP_np = Intg_L(obs.nHalf, src.pHalf, freq);
	    if (src.nHalf) iP_nn = Intg_L(obs.nHalf, src.nHalf, freq);
	}
	result = iP_pp - iP_pn - iP_np + iP_nn;

	return result;
}

DComplex mfie(const FullBF &obs, const FullBF &src, const double freq, const double er)
{
	DComplex iP_pp(0), iP_pn(0), iP_np(0), iP_nn(0);
	DComplex result;

	if (obs.pHalf) {
	    if (src.pHalf) iP_pp = Intg_KN(obs.pHalf, src.pHalf, freq) +
                                0.5*Intg_hbf(obs.pHalf, src.pHalf);
	    if (src.nHalf) iP_pn = Intg_KN(obs.pHalf, src.nHalf, freq) +
                                0.5*Intg_hbf(obs.pHalf, src.nHalf);
	}
	if (obs.nHalf) {
	    if (src.pHalf) iP_np = Intg_KN(obs.nHalf, src.pHalf, freq) +
                                0.5*Intg_hbf(obs.nHalf, src.pHalf);
	    if (src.nHalf) iP_nn = Intg_KN(obs.nHalf, src.nHalf, freq) +
                                0.5*Intg_hbf(obs.nHalf, src.nHalf);
	}
	result = iP_pp - iP_pn - iP_np + iP_nn;

	return result;
}

DComplex cfie(const FullBF &obs, const FullBF &src, const double freq, const double er)
{
    DComplex efie_value = efie(obs, src, freq);
    DComplex mfie_value = mfie(obs, src, freq);
    DComplex eta = ETA0;

    DComplex result = efie_value + eta*mfie_value;
	return result;
}

DComplex
pmchwt(const FullBF &obs, const FullBF &src, const double freq,
        const double er_in, const double er_ex)
{
	DComplex iP_pp(0), iP_pn(0), iP_np(0), iP_nn(0);
	DComplex result;

	if (obs.pHalf) {
	    if (src.pHalf) iP_pp = Intg_KN(obs.pHalf, src.pHalf, freq);
	    if (src.nHalf) iP_pn = Intg_KN(obs.pHalf, src.nHalf, freq);
	}
	if (obs.nHalf) {
	    if (src.pHalf) iP_np = Intg_KN(obs.nHalf, src.pHalf, freq);
	    if (src.nHalf) iP_nn = Intg_KN(obs.nHalf, src.nHalf, freq);
	}

	result = iP_pp - iP_pn - iP_np + iP_nn;

	return result;
}
#endif // 0



DComplex
Intg_L(const HalfBF * hbf1, const HalfBF * hbf2, const double freq, const DComplex er)
{
    /// 需知， 函数指针可以解决函数重载的实现问题
    DComplex result(0), vecResult(0), scaResult(0);
	double omega(2*PI*freq);
	DComplex epsilon_r(er);
    DComplex coefA(CJ*omega*MU0), coefV(-CJ/(epsilon_r*EPSILON0*omega));

    int gauss_rule;
    int gauss_num;
    double const *gauss_weight=0;
    double const *gauss_coef[MAX_GAUSSNUM];
    double len1 = 0, len2 = 0;

    len1 = (hbf1->type == Tri) ? 0.5*hbf1->hRwg->len : hbf1->hRtp->ly;
    len2 = (hbf2->type == Tri) ? 0.5*hbf2->hRwg->len : hbf2->hRtp->ly;

    if (hbf1->type == Tri) {
        gauss_rule = DEFAULT_TRI_GAUSS_RULE;
        gauss_num = tri_gauss_num[gauss_rule];
        gauss_weight = tri_gauss_weight[gauss_rule];
        for (int i=0;i<gauss_num;i++) gauss_coef[i] = tri_gauss_coef[gauss_rule][i];

        for (int i=0;i<gauss_num;i++) {
            DComplex scaT(0);
            Complex3D vecT(0, 0, 0);
            Point3D pg;
            pg = get_gauss_point(hbf1, gauss_coef[i], pg);

            LO_normal(hbf1, hbf2, freq, epsilon_r, pg, vecT, scaT);

            vecResult += dotprod(pg - hbf1->hRwg->v0, vecT) * gauss_weight[i];
            scaResult += scaT * gauss_weight[i];
        }
    }

    result = len1 * len2 * (coefA*vecResult + 4.0*coefV*scaResult);
	return result;
}

DComplex
Intg_KN(const HalfBF * hbf1, const HalfBF * hbf2, const double freq, const DComplex er)
{
    DComplex result(0);

    int gauss_rule;
    int gauss_num;
    double const *gauss_weight=0;
    double const *gauss_coef[MAX_GAUSSNUM];
    double len1 = hbf1->hRwg->len, len2 = hbf2->hRwg->len;

    if (hbf1->type==Tri) {
        const halfrwgBF *hRwg = hbf1->hRwg;

        gauss_rule = DEFAULT_TRI_GAUSS_RULE;
        gauss_num = tri_gauss_num[gauss_rule];
        gauss_weight = tri_gauss_weight[gauss_rule];
        for (int i=0;i<gauss_num;i++) gauss_coef[i] = tri_gauss_coef[gauss_rule][i];

        for (int i=0;i<gauss_num;++i) {
            Complex3D vecKO;
            Point3D gp;
            gp = get_gauss_point(hbf1, gauss_coef[i], gp);

            vecKO = KO_normal(hbf1, hbf2, freq, er, gp);

            Point3D tmp1 = gp - hRwg->v0;
            Point3D nvec = hRwg->nvec;
            Point3D tmp2 = crossprod(tmp1, nvec);

            result += dotprod(tmp2, vecKO) * gauss_weight[i];
        }

        result *= len1 * len2 * 0.25;
    }

    return result;
}

DComplex
Intg_KT(const HalfBF * hbf1, const HalfBF * hbf2, const double freq, const DComplex er)
{
    DComplex result(0);

    int gauss_rule;
    int gauss_num;
    double const *gauss_weight=0;
    double const *gauss_coef[MAX_GAUSSNUM];
    double len1 = hbf1->hRwg->len, len2 = hbf2->hRwg->len;

    if (hbf1->type==Tri) {
        const halfrwgBF *hRwg = hbf1->hRwg;

        gauss_rule = DEFAULT_TRI_GAUSS_RULE;
        gauss_num = tri_gauss_num[gauss_rule];
        gauss_weight = tri_gauss_weight[gauss_rule];
        for (int i=0;i<gauss_num;i++) gauss_coef[i] = tri_gauss_coef[gauss_rule][i];

        for (int i=0;i<gauss_num;++i) {
            Complex3D vecKO;
            Point3D gp;
            gp = get_gauss_point(hbf1, gauss_coef[i], gp);

            vecKO = KO_normal(hbf1, hbf2, freq, er, gp);

            Point3D tmp2 = gp - hRwg->v0;

            result += dotprod(tmp2, vecKO) * gauss_weight[i];
        }

        result *= len1 * len2 * 0.25;
    }

    return result;
}

double
Intg_hbf(const HalfBF * hbf1, const HalfBF * hbf2)
{
	double ans(0);

	if (hbf1->face_label == hbf2->face_label) {
        if (hbf1->type == Tri) {
            const halfrwgBF *hrwg1 = hbf1->hRwg;
            const halfrwgBF *hrwg2 = hbf2->hRwg;
            ans = dotprod(hrwg1->ctr-hrwg1->v0, hrwg2->ctr-hrwg2->v0);
            ans *= hrwg1->len * hrwg2->len / (4.0*hrwg1->area);
        }
	}

	return ans;
}

DComplex
innerProduct_UE(const FullBF &rbf, double freq)
{
    double wavenum = 2*PI*freq /SPEED_OF_LIGHT;
    Point3D incDirt(0, 0, -1);
    Point3D polarDirt(1, 0, 0);

    DComplex result;

    result = (rbf.pHalf->hRwg->len/2.0) * Em *
    (dotprod(rbf.pHalf->hRwg->ctr-rbf.pHalf->hRwg->v0, polarDirt) * exp(-CJ*wavenum*dotprod(incDirt, rbf.pHalf->hRwg->ctr))-
     dotprod(rbf.nHalf->hRwg->ctr-rbf.nHalf->hRwg->v0, polarDirt) * exp(-CJ*wavenum*dotprod(incDirt, rbf.nHalf->hRwg->ctr)));

    return result;
}
DComplex
innerProduct_UTH(const FullBF &rbf, double freq)
{
    double wavenum = 2*PI*freq /SPEED_OF_LIGHT;
    Point3D incDirt(0, 0, -1);
    Point3D polarDirt(0, -1, 0);

    DComplex result;

    result = (rbf.pHalf->hRwg->len/2.0) * Hm *
    (dotprod(rbf.pHalf->hRwg->ctr-rbf.pHalf->hRwg->v0, polarDirt) * exp(-CJ*wavenum*dotprod(incDirt, rbf.pHalf->hRwg->ctr))-
     dotprod(rbf.nHalf->hRwg->ctr-rbf.nHalf->hRwg->v0, polarDirt) * exp(-CJ*wavenum*dotprod(incDirt, rbf.nHalf->hRwg->ctr)));

    return result;
}

DComplex
innerProduct_UH(const FullBF &rbf, double freq)
{
    double wavenum = 2*PI*freq /SPEED_OF_LIGHT;
    Point3D incDirt(0, 0, -1);
    Point3D polarDirt(0, -1, 0);
    const halfrwgBF *prwg = rbf.pHalf->hRwg;
    const halfrwgBF *nrwg = rbf.nHalf->hRwg;

    Point3D temp1 = prwg->ctr - prwg->v0;
    Point3D temp2 = nrwg->ctr - nrwg->v0;

    Point3D trouble1 = crossprod(prwg->nvec, polarDirt);
    Point3D trouble2 = crossprod(nrwg->nvec, polarDirt);

    DComplex value1 = dotprod(temp1, trouble1)*exp(-CJ*wavenum*dotprod(incDirt, prwg->ctr));
    DComplex value2 = dotprod(temp2, trouble2)*exp(-CJ*wavenum*dotprod(incDirt, nrwg->ctr));

    DComplex result = prwg->len*0.5*Hm* (value1 - value2);

    return result;
}

void
LO_normal(const HalfBF *hbf1, const HalfBF *hbf2, const double freq,
        const DComplex er, const Point3D &obs_gauss_point, Complex3D &l1, DComplex &l2)
{
    DComplex wavenum = 2*PI*freq*sqrt(er)/SPEED_OF_LIGHT;

    int gauss_rule = DEFAULT_TRI_GAUSS_RULE;
    int gauss_num = tri_gauss_num[gauss_rule];
    double const *gauss_weight=tri_gauss_weight[gauss_rule];
    double const *gauss_coef[MAX_GAUSSNUM] = {0};
    for (int i=0;i<gauss_num;i++) gauss_coef[i] = tri_gauss_coef[gauss_rule][i];

    const Point3D ac(hbf1->hRwg->ctr), bc(hbf2->hRwg->ctr);

    if (hbf1->hRwg->face_label == hbf2->hRwg->face_label) {
        double Isp[3], Is;
        Complex3D temp(0, 0, 0);
        Isp_and_Is(Isp, &Is, hbf2->hRwg->v0, hbf2->hRwg->v1, hbf2->hRwg->v2, obs_gauss_point);

        for (int j=0;j<gauss_num;++j) {
            Point3D src_gauss_point;
            src_gauss_point = get_gauss_point(hbf2, gauss_coef[j], src_gauss_point);
            temp -= CJ * wavenum * (src_gauss_point-hbf2->hRwg->v0) * gauss_weight[j];
        }

        for (int i=0;i<3;i++) {
            l1[i] = 1/(4*PI) * (1/(hbf2->hRwg->area)*(Isp[i] + Is*(obs_gauss_point[i]-(hbf2->hRwg->v0[i]))) + temp[i]);
        }
        l2 = (1.0/(4*PI)) * (-CJ*wavenum + 1.0*Is/hbf2->hRwg->area);
    }
    else {
        for (int i=0;i<gauss_num;++i) {
            Point3D src_gauss_point;
            src_gauss_point = get_gauss_point(hbf2, gauss_coef[i], src_gauss_point);
            DComplex green_value = Green0(obs_gauss_point, src_gauss_point, wavenum);
            l1 += (gauss_weight[i]*green_value)*(src_gauss_point-hbf2->hRwg->v0);
            l2 += gauss_weight[i]*green_value;
        }
    }
}

void
LO_normal(const HalfBF *hbf, const double freq,
        const DComplex er, const Point3D &obs_gauss_point, Complex3D &l1, DComplex &l2)
{
    DComplex wavenum = 2*PI*freq*sqrt(er)/SPEED_OF_LIGHT;

//    cout << "in LO ntest " << endl;

    int gauss_rule = DEFAULT_TRI_GAUSS_RULE;
    int gauss_num = tri_gauss_num[gauss_rule];
    double const *gauss_weight=tri_gauss_weight[gauss_rule];
    double const *gauss_coef[MAX_GAUSSNUM] = {0};
    for (int i=0;i<gauss_num;i++) gauss_coef[i] = tri_gauss_coef[gauss_rule][i];

    for (int i=0;i<gauss_num;++i) {
        Point3D src_gauss_point;
        src_gauss_point = get_gauss_point(hbf, gauss_coef[i], src_gauss_point);
        DComplex green_value = Green0(obs_gauss_point, src_gauss_point, wavenum);
        l1 += (gauss_weight[i]*green_value)*(src_gauss_point-hbf->hRwg->v0);
        l2 += gauss_weight[i]*green_value;
    }
}

Complex3D
LO_vec_ntest(const HalfBF *hbf, const double freq,
             const DComplex er, const Point3D &obs_gauss_point)
{
    DComplex wavenum = 2*PI*freq*sqrt(er)/SPEED_OF_LIGHT;
    Complex3D ans(0, 0, 0);

    int gauss_rule = DEFAULT_TRI_GAUSS_RULE;
    int gauss_num = tri_gauss_num[gauss_rule];
    double const *gauss_weight=tri_gauss_weight[gauss_rule];
    double const *gauss_coef[MAX_GAUSSNUM] = {0};
    for (int i=0;i<gauss_num;i++) gauss_coef[i] = tri_gauss_coef[gauss_rule][i];

    for (int i=0;i<gauss_num;++i) {
        Point3D src_gauss_point;
        src_gauss_point = get_gauss_point(hbf, gauss_coef[i], src_gauss_point);
        DComplex green_value = Green0(obs_gauss_point, src_gauss_point, wavenum);
        ans += (gauss_weight[i]*green_value)*crossprod(src_gauss_point-hbf->hRwg->v0, hbf->hRwg->nvec);
    }

    return ans;
}

Complex3D
LO_vec_ttest(const HalfBF *hbf, const double freq,
             const DComplex er, const Point3D &obs_gauss_point)
{
    DComplex wavenum = 2*PI*freq*sqrt(er)/SPEED_OF_LIGHT;
    Complex3D ans(0, 0, 0);

    int gauss_rule = DEFAULT_TRI_GAUSS_RULE;
    int gauss_num = tri_gauss_num[gauss_rule];
    double const *gauss_weight=tri_gauss_weight[gauss_rule];
    double const *gauss_coef[MAX_GAUSSNUM] = {0};
    for (int i=0;i<gauss_num;i++) gauss_coef[i] = tri_gauss_coef[gauss_rule][i];

    for (int i=0;i<gauss_num;++i) {
        Point3D src_gauss_point;
        src_gauss_point = get_gauss_point(hbf, gauss_coef[i], src_gauss_point);
        DComplex green_value = Green0(obs_gauss_point, src_gauss_point, wavenum);
        ans += (gauss_weight[i]*green_value)*(src_gauss_point-hbf->hRwg->v0);
    }

    return ans;
}

Complex3D
KO_normal(const HalfBF *hbf1, const HalfBF *hbf2, const double freq,
        const DComplex er, const Point3D &obs_gauss_point)
{
    const DComplex wavenum(2*PI*freq*sqrt(er)/SPEED_OF_LIGHT);
    Complex3D ans(0, 0, 0);

    int gauss_rule;
    int gauss_num;
    double const *gauss_weight;
    double const *gauss_coef[MAX_GAUSSNUM] = {0};

    if (hbf1->face_label == hbf2->face_label) {
        const halfrwgBF * hRwg = hbf2->hRwg;
        double Isp3[3], Isp[3];

        Isp3_and_Isp(Isp3, Isp, hRwg->v0, hRwg->v1, hRwg->v2, obs_gauss_point);
        Complex3D Isp3_oop(Isp3[0], Isp3[1], Isp3[2]);
        Complex3D Isp_oop(Isp[0], Isp[1], Isp[2]);
        Isp3_oop = Isp3_oop + 0.5*wavenum*wavenum*Isp_oop;

        Point3D temp = hRwg->v0 - hRwg->ctr;
        Complex3D majorPart = crossprod(temp, Isp3_oop);

        temp = obs_gauss_point - hRwg->v0;
        Point3D reduntPart = crossprod(temp, hRwg->nvec);

        Complex3D majorPart_c(majorPart[0], majorPart[1], majorPart[2]);
        Complex3D reduntPart_c(reduntPart[0], reduntPart[1], reduntPart[2]);

//            ans = (1.0/(4*PI) * majorPart_c - 0.5*flag*reduntPart_c)/(DComplex)hRwg->area;
//            外部用 1/2 的 奇异值消除
        ans = (1.0/(4*PI) * majorPart_c)/(DComplex)hRwg->area;
    }
    else {
            const halfrwgBF * hRwg = hbf2->hRwg;

            gauss_rule = DEFAULT_TRI_GAUSS_RULE;
            gauss_num = tri_gauss_num[gauss_rule];
            gauss_weight=tri_gauss_weight[gauss_rule];
            for (int i=0;i<gauss_num;i++) gauss_coef[i] = tri_gauss_coef[gauss_rule][i];

            for (int j=0;j<gauss_num;j++) {
                Point3D src_gauss_point;
                src_gauss_point = get_gauss_point(hbf2, gauss_coef[j], src_gauss_point);
                Point3D temp2 = src_gauss_point-hRwg->v0;
                Complex3D temp1 = nablaGreen0(obs_gauss_point, src_gauss_point, wavenum);
                Complex3D temp3 = crossprod(temp2, temp1);
                ans += temp3 * gauss_weight[j];
            }
    }

    return ans;
}

Complex3D
KO_normal(const HalfBF *hbf, const double freq,
        const DComplex er, const Point3D &obs_gauss_point)
{
    const DComplex wavenum(2*PI*freq*sqrt(er)/SPEED_OF_LIGHT);
    Complex3D ans(0, 0, 0);

    int gauss_rule = DEFAULT_TRI_GAUSS_RULE;
    int gauss_num = tri_gauss_num[gauss_rule];
    double const *gauss_weight=tri_gauss_weight[gauss_rule];
    double const *gauss_coef[MAX_GAUSSNUM] = {0};
    for (int i=0;i<gauss_num;i++) gauss_coef[i] = tri_gauss_coef[gauss_rule][i];

    const halfrwgBF * hRwg = hbf->hRwg;

    for (int j=0;j<gauss_num;j++) {
        Point3D src_gauss_point;
        src_gauss_point = get_gauss_point(hbf, gauss_coef[j], src_gauss_point);
        Point3D temp2 = src_gauss_point-hRwg->v0;
        Complex3D temp1 = nablaGreen0(obs_gauss_point, src_gauss_point, wavenum);
        Complex3D temp3 = crossprod(temp2, temp1);
        ans += temp3 * gauss_weight[j];
    }

    return ans;
}

Complex3D
KO_normal_with_n(const HalfBF *hbf, const HalfBF *hbf_obs, const double freq,
        const DComplex er, const Point3D &obs_gauss_point)
{
    const DComplex wavenum(2*PI*freq*sqrt(er)/SPEED_OF_LIGHT);
    Complex3D ans(0, 0, 0);

    int gauss_rule = DEFAULT_TRI_GAUSS_RULE;
    int gauss_num = tri_gauss_num[gauss_rule];
    double const *gauss_weight=tri_gauss_weight[gauss_rule];
    double const *gauss_coef[MAX_GAUSSNUM] = {0};
    for (int i=0;i<gauss_num;i++) gauss_coef[i] = tri_gauss_coef[gauss_rule][i];

    const halfrwgBF * hRwg = hbf->hRwg;

    for (int j=0;j<gauss_num;j++) {
        Point3D src_gauss_point;
        src_gauss_point = get_gauss_point(hbf, gauss_coef[j], src_gauss_point);
        Point3D temp2 = src_gauss_point-hRwg->v0;
        Complex3D temp1 = nablaGreen0(obs_gauss_point, src_gauss_point, wavenum);
        Complex3D temp3 = crossprod(temp2, temp1);
        ans += temp3 * gauss_weight[j];
    }

    return crossprod(hbf_obs->hRwg->nvec, ans);
}

