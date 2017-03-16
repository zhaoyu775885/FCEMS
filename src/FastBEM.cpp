#include "FastBEM.h"
#include <algorithm>
#include "Math/ChebItp.h"
#include "h2lib/eigensolvers.h"
#include "h2lib/factorizations.h"
#include <chrono>
#include <queue>

typedef std::chrono::high_resolution_clock Clock;

#define TST_HMATRIX_ERROR 0
#define USE_ACA 0
#define USE_HCA 1
#define USE_FULL_ACA 0
#define USE_MERGE 0
#define INTP_ORDER 4
#define SHEET_SIZE 1E-6
#define USE_KERNEL_SVD 0
#define SEPERATELY 0
#define COST_REDUCTION 0

int adepth = 0;
int print = 0;
double total_diff(6);
const DComplex epsilon_ex(1.0);
const DComplex epsilon_in(4.0);

FastBem::FastBem(const BasisFunc &func_space, const FastBemConf &config) :
            func_space(func_space), g_cluster(0), g_hmatrix(0),
            uvec(0), nrow(func_space.get_dof()), nrhs(1),
            leaf_size_threshold(config.max_leaf_size), kesi(config.eta),
            aca_threshold(config.threshold)
{
    pclustergeometry pcg = new_clustergeometry(3, nrow);
	for (int i=0;i<nrow;i++) {
		pcg->x[i] = new double [3];
		Point3D point_mid;
		if (func_space[i].pHalf->type == Tri) {
            const halfrwgBF *hRwg = func_space[i].pHalf->hRwg;
            point_mid = (hRwg->v1 + hRwg->v2) / 2.0;
		}
		else if (func_space[i].pHalf->type == Rect) {
            const halfrtpBF *hRtp = func_space[i].pHalf->hRtp;
            point_mid = (hRtp->v1 + hRtp->v2) / 2.0;
		}
        pcg->x[i][0] = point_mid[0];
        pcg->x[i][1] = point_mid[1];
        pcg->x[i][2] = point_mid[2];
	}
	index_tree = new int [nrow];
	aux_index_tree = new int [2*nrow];
	for (int i=0;i<nrow;i++) index_tree[i] = i;

	clustermode CMode = H2_ADAPTIVE;
	g_cluster = build_cluster(pcg, nrow, index_tree, config.max_leaf_size, CMode);

	/// The clustering should began at the longest coordinate.
	/// Like, [-1, 1]x[-2, 2]x[-4, 4]
	/// In this case, we should start at Z axis.
	for (int i=0;i<nrow;i++) aux_index_tree[i] = aux_index_tree[i+nrow] = index_tree[i];
	aux_cluster = build_regular_cluster(pcg, 2*nrow, aux_index_tree, config.max_leaf_size, 0);

	max_rank_in_rk = 0;
    frob_norm_2_total = 0;
    frob_norm_error_2_abs = 0;
    frob_norm_error_2_rel = 0;

    svd_threshold = sqrt(aca_threshold);
    #if !USE_ACA
    aca_threshold *= 0.1;
    #endif // USE_ACA

    #if USE_ACA
    cout << "use pure ACA" << endl;
    #elif USE_HCA
    cout << "use HCA" << endl;
    #endif // USE_ACA
}

void
FastBem::build_equation(const double freq, const char c)
{
    build_system_matrix(freq, c);
    build_excitation(freq, c);
}

void
FastBem::build_system_matrix(const double freq, const char c)
{
    pblock g_block = build_block(g_cluster, g_cluster, &kesi, admissible_2_cluster);

    Integration intg(func_space, freq, epsilon_in, epsilon_ex);

    if (c=='p' || c=='P') {
        for (int i=0;i<nrow;++i) {
            aux_index_tree[i] = index_tree[i];
            aux_index_tree[i+nrow] = index_tree[i] + nrow;
        }
        nrow *= 2;
        g_hmatrix = new_super_hmatrix(aux_cluster, aux_cluster, 2, 2);
        g_hmatrix->son[0] = build_from_block_hmatrix(g_block, 0);
        g_hmatrix->son[1] = build_from_block_hmatrix(g_block, 0);
        g_hmatrix->son[2] = build_from_block_hmatrix(g_block, 0);
        g_hmatrix->son[3] = build_from_block_hmatrix(g_block, 0);
        g_hmatrix->son[0]->refs = 1;
        g_hmatrix->son[1]->refs = 1;
        g_hmatrix->son[2]->refs = 1;
        g_hmatrix->son[3]->refs = 1;
        g_hmatrix->refs = 1;

        #if !USE_ACA

        assemble_hmatrix(g_hmatrix->son[0], g_hmatrix->son[1], g_hmatrix->son[2], g_hmatrix->son[3], intg);

        #else
        assemble_hmatrix(g_hmatrix->son[0], intg, 'w');
        assemble_hmatrix(g_hmatrix->son[1], intg, 'y');
        assemble_hmatrix(g_hmatrix->son[2], intg, 'x');
        assemble_hmatrix(g_hmatrix->son[3], intg, 'z');
        #endif // HCA_PMCHWT
    }
    else {
        g_hmatrix = build_from_block_hmatrix(g_block, 0);
        g_hmatrix->refs = 1;

        #if 0
        cout << "begin merge matrix" << endl;
        merge_rkmatrix(g_hmatrix);
        merge_fmatrix(g_hmatrix);
        cout << "merge over" << endl;
        #endif // 0

        assemble(intg, c);
    }

    #if COST_REDUCTION
    cost_reduction(1e-1);
    #endif // COST_REDUCTION

    #if TST_HMATRIX_ERROR
    if (frob_norm_2_total != 0.0) frob_norm_error_2_rel = sqrt(frob_norm_error_2_abs/frob_norm_2_total);
    #endif // TST_HMATRIX_ERROR
}

void
FastBem::build_excitation(const double freq, const char c)
{
    uvec = new_avector(nrow);

    int nEdge = func_space.get_dof();

    for (int i=0;i<nEdge;i++) {
        if (c=='e' || c=='E') {
            uvec->v[i] = innerProduct_UE(func_space[i], freq);
        }
        else if (c=='m' || c=='M') {
            uvec->v[i] = innerProduct_UH(func_space[i], freq);
        }
        else if (c=='c' || c=='C') {
            uvec->v[i] = innerProduct_UE(func_space[i], freq) +
                         innerProduct_UH(func_space[i], freq) * ETA0;
        }
        else if (c=='p' || c=='P') {
            uvec->v[i] = innerProduct_UE(func_space[i], freq);
            uvec->v[i+nEdge] = innerProduct_UTH(func_space[i], freq);
        }
    }
}

void
FastBem::assemble(const Integration &intg, const char c)
{
    #if SEPERATELY
	auto t1 = Clock::now();
    assemble_hmatrix_fmatrix(g_hmatrix, intg, c);
    auto t2 = Clock::now();
    cout << "assemble full matrices : " <<std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count() / 1000.0
         << " s" << endl;
    assemble_hmatrix_rkmatrix(g_hmatrix, intg, c);
    auto t3 = Clock::now();
    cout << "assemble rk matrices : " <<std::chrono::duration_cast<std::chrono::milliseconds>(t3 - t2).count() / 1000.0
         << " s" << endl;
    #else
    assemble_hmatrix(g_hmatrix, intg, c);
    #endif // SEPERATELY
}

phmatrix
FastBem::assemble_hmatrix_fmatrix(phmatrix hm, const Integration &intg, const char c)
{
	int rsons = hm->rsons;
	int csons = hm->csons;

	if (hm->son) {
		adepth++;
		for (int j=0;j<csons;j++) {
			for (int i=0;i<rsons;i++) {
				assemble_hmatrix(hm->son[i+j*rsons], intg, c);
			}
		}
		adepth--;
	}
	else if (hm->f){
		assert(hm->r == NULL);
		hm = assemble_fmatrix(hm, intg, c);
	}

	return hm;
}

phmatrix
FastBem::assemble_hmatrix_rkmatrix(phmatrix hm, const Integration &intg, const char c)
{
	int rsons = hm->rsons;
	int csons = hm->csons;

	if (hm->son) {
		adepth++;
		for (int j=0;j<csons;j++) {
			for (int i=0;i<rsons;i++) {
				assemble_hmatrix(hm->son[i+j*rsons], intg, c);
			}
		}
		adepth--;
	}
	else if (hm->r) {
        assert(hm->f == NULL);
		hm = assemble_rkmatrix(hm, intg, c);
	}

	return hm;
}

phmatrix
FastBem::assemble_hmatrix(phmatrix hm, const Integration &intg, const char c)
{
	int rsons = hm->rsons;
	int csons = hm->csons;

	if (hm->son) {
		adepth++;
		for (int j=0;j<csons;j++) {
			for (int i=0;i<rsons;i++) {
				assemble_hmatrix(hm->son[i+j*rsons], intg, c);
			}
		}
		adepth--;
	}
	else if (hm->r) {
		hm = assemble_rkmatrix(hm, intg, c);
	}
	else {
		assert(hm->f != NULL);
		hm = assemble_fmatrix(hm, intg, c);
	}

	return hm;
}

void
FastBem::assemble_hmatrix(phmatrix hm0, phmatrix hm1, phmatrix hm2, phmatrix hm3, const Integration &intg)
{
	int rsons = hm0->rsons;
	int csons = hm0->csons;

	if (hm0->son) {
		adepth++;
		for (int j=0;j<csons;j++) {
			for (int i=0;i<rsons;i++) {
				assemble_hmatrix(hm0->son[i+j*rsons], hm1->son[i+j*rsons], hm2->son[i+j*rsons], hm3->son[i+j*rsons], intg);
			}
		}
		adepth--;
	}
	else if (hm0->r) {
		assemble_rkmatrix(hm0, hm1, hm2, hm3, intg);
	}
	else {
		assert(hm0->f != NULL);
		assemble_fmatrix(hm0, hm1, hm2, hm3, intg);
	}
}

phmatrix
FastBem::assemble_rkmatrix(phmatrix hm, const Integration &intg, const char c)
{
    #if USE_ACA
	aca(hm, intg, c);
	ptruncmode tm = new_truncmode();
	trunc_rkmatrix(tm, svd_threshold, hm->r);
	#elif USE_HCA
	hca(hm, intg, c);
	#endif // USE_ACA

	if (max_rank_in_rk < hm->r->k) max_rank_in_rk = hm->r->k;
//	cout << "After trunction, the rank= " << max_rank_in_rk << endl<<endl;
	return hm;
}

void
FastBem::assemble_rkmatrix(phmatrix hm0, phmatrix hm1, phmatrix hm2, phmatrix hm3, const Integration &intg)
{
	hca(hm0, hm1, hm2, hm3, intg);

//	if (max_rank_in_rk < hm->r->k) max_rank_in_rk = hm->r->k;
//	cout << "After trunction, the rank= " << max_rank_in_rk << endl<<endl;
}

phmatrix
FastBem::assemble_fmatrix(phmatrix hm, const Integration &intg, const char c)
{
    pamatrix f = hm->f;
    pccluster rc = hm->rc;
    pccluster cc = hm->cc;

    #pragma omp parallel for
	for (int j=0;j<f->cols;j++) {
		for (int i=0;i<f->rows;i++) {
            f->a[i+j*f->ld] = intg.xfie(rc->idx[i], cc->idx[j], c);
		}
	}

	for (int j=0;j<f->cols;j++) {
		for (int i=0;i<f->rows;i++) {
            #if TST_HMATRIX_ERROR
            DComplex tst = f->a[i+j*f->ld];
            frob_norm_2_total += abs(tst*tst);
            #endif // TST_HMATRIX_ERROR
		}
	}

	return hm;
}

void
FastBem::assemble_fmatrix(phmatrix hm0, phmatrix hm1, phmatrix hm2, phmatrix hm3, const Integration &intg)
{
    pamatrix f = hm0->f;
    pccluster rc = hm0->rc;
    pccluster cc = hm0->cc;

    #pragma omp parallel for
	for (int j=0;j<f->cols;j++) {
		for (int i=0;i<f->rows;i++) {
            #if 1
            DComplex a1 = intg.efie(rc->idx[i], cc->idx[j], epsilon_in);
            DComplex a2 = intg.efie(rc->idx[i], cc->idx[j], epsilon_ex);
            DComplex b1 = intg.pmchwt(rc->idx[i], cc->idx[j], epsilon_in);
            DComplex b2 = intg.pmchwt(rc->idx[i], cc->idx[j], epsilon_ex);
            hm0->f->a[i+j*f->ld] = a1+a2;
            hm1->f->a[i+j*f->ld] = b1+b2;
            hm2->f->a[i+j*f->ld] = -b1-b2;
            DComplex eta0_2 = ETA0*ETA0;
            hm3->f->a[i+j*f->ld] = (a1*epsilon_in+a2*epsilon_ex) / eta0_2;
            #else
            hm0->f->a[i+j*f->ld] = intg.xfie(rc->idx[i], cc->idx[j], 'w');
            hm1->f->a[i+j*f->ld] = intg.xfie(rc->idx[i], cc->idx[j], 'y');
            hm2->f->a[i+j*f->ld] = intg.xfie(rc->idx[i], cc->idx[j], 'x');
            hm3->f->a[i+j*f->ld] = intg.xfie(rc->idx[i], cc->idx[j], 'z');
            #endif // 1
		}
	}

	for (int j=0;j<f->cols;j++) {
		for (int i=0;i<f->rows;i++) {
            #if TST_HMATRIX_ERROR
            DComplex tst0 = f->a[i+j*f->ld];
            DComplex tst1 = hm1->f->a[i+j*f->ld];
            DComplex tst2 = hm2->f->a[i+j*f->ld];
            DComplex tst3 = hm3->f->a[i+j*f->ld];
            frob_norm_2_total += abs(tst0*tst0) + abs(tst1*tst1) + abs(tst2*tst2) + abs(tst3*tst3);
            #endif // TST_HMATRIX_ERROR
		}
	}
}

void
FastBem::precondition()
{
    phmatrix g_hmatrix_inv = clone_hmatrix(g_hmatrix);
    phmatrix fuck_hmat = clone_hmatrix(g_hmatrix);
    ptruncmode tm = new_truncmode();

    phmatrix aux_hmatrix = clonestructure_hmatrix(g_hmatrix);
    clear_hmatrix(aux_hmatrix);
    invert_hmatrix(g_hmatrix_inv, aux_hmatrix, tm, 1e-3);
    clear_hmatrix(g_hmatrix);
    addmul_hmatrix(1.0, false, g_hmatrix_inv, false, fuck_hmat, tm, 1e-3, g_hmatrix);
//
//    g_hmatrix = aux_hmatrix;

    pavector newuvec = new_avector(uvec->dim);
    for (int i=0;i<nrow;++i) {
        newuvec->v[i] = uvec->v[g_cluster->idx[i]];
    }
    for (int i=0;i<nrow;++i) {
        uvec->v[i] = newuvec->v[i];
    }

    fastaddeval_hmatrix_avector(1.0, g_hmatrix_inv, CblasNoTrans, uvec, 0.0, newuvec);

    for (int i=0;i<nrow;++i) {
        uvec->v[g_cluster->idx[i]] = newuvec->v[i];
    }

    del_hmatrix(g_hmatrix_inv);
//    del_hmatrix(aux_hmatrix);
}

void
FastBem::direct_solve_vec()
{
//    cout << "begin direct H-LU decomposition..." << endl;
    ptruncmode tm = new_truncmode();
    lrdecomp_hmatrix(g_hmatrix, tm, 1e-2);
//    cout << "H-LU decomposition finished." << endl;
    lrsolve_hmatrix_avector(false, g_hmatrix, uvec);
//    cout << "direct method finished." << endl;
}

void
FastBem::iterative_solve_vec()
{
    pavector newuvec = new_avector(uvec->dim);
    for (int i=0;i<nrow;++i) {
        newuvec->v[i] = uvec->v[i];
    }

    cout << "Begin CG iteration..." << endl;
    int steps = hmat_cg(g_hmatrix, newuvec, uvec);
    cout << "CG iteration steps = " << steps << endl;
    cout << "iterative method finished." << endl;
}

void
FastBem::print_zmat(const string filename) const
{
    int n = g_cluster->size;
    ofstream fp(filename.c_str());
    for (int i=0;i<n;++i)
        fp << g_hmatrix->f->a[i] << endl;
    fp.close();
}

void
FastBem::print_umat(const string filename) const
{
    int n = g_cluster->size;
    ofstream fp(filename.c_str());
    for (int i=0;i<n;++i)
        fp << uvec->v[i] << endl;
    fp.close();
}

void
FastBem::hca(phmatrix hm0, phmatrix hm1, phmatrix hm2, phmatrix hm3, const Integration &intg)
{
//    cout << hca_times++ << ": use HCA" << endl;

    /** parameter for integral **/
    const double freq(intg.get_freq());
    const double omega(2*PI*freq);
    const DComplex wavenum_ex = omega*(1.0/Constants::SPEED_OF_LIGHT)*sqrt(epsilon_ex);
    const DComplex wavenum_in = omega*(1.0/Constants::SPEED_OF_LIGHT)*sqrt(epsilon_in);
    const double lambda_ex(2*Constants::PI/wavenum_ex.real());
    const double lambda_in(2*Constants::PI/wavenum_in.real());
    const DComplex coefA_ex(CJ*omega*MU0), coefV_ex(-CJ/(epsilon_ex*EPSILON0*omega));
    const DComplex coefA_in(CJ*omega*MU0), coefV_in(-CJ/(epsilon_in*EPSILON0*omega));

    /** the bounding box is crucial for hca **/
    const double *rbox_min = hm0->rc->bmin;
    const double *rbox_max = hm0->rc->bmax;
    const double *cbox_min = hm0->cc->bmin;
    const double *cbox_max = hm0->cc->bmax;

    /** fixed points number may be too small or too large
     ** use the space resolution method to select the size adaptively **/
    double rxa = rbox_min[0], rxb = rbox_max[0];
    double rya = rbox_min[1], ryb = rbox_max[1];
    double rza = rbox_min[2], rzb = rbox_max[2];
    double cxa = cbox_min[0], cxb = cbox_max[0];
    double cya = cbox_min[1], cyb = cbox_max[1];
    double cza = cbox_min[2], czb = cbox_max[2];
    const double rbox_xlen = rxb - rxa;
    const double rbox_ylen = ryb - rya;
    const double rbox_zlen = rzb - rza;
    const double cbox_xlen = cxb - cxa;
    const double cbox_ylen = cyb - cya;
    const double cbox_zlen = czb - cza;

    const int r_nx = get_inter_num(rbox_xlen, lambda_ex);
    const int r_ny = get_inter_num(rbox_ylen, lambda_ex);
    const int r_nz = get_inter_num(rbox_zlen, lambda_ex);
//    cout << r_nx << ", " << r_ny << ", " << r_nz << endl;
    /** build the point set for row cluster **/
    vector<Point3D> rbox_nodes_set;
    const int rbox_nodes_num(r_nx*r_ny*r_nz);
    ChebItp3D row_cheb(r_nx, r_ny, r_nz);
    for (int i=0;i<r_nx;++i) {
        double px = row_cheb.get_cheb_point_x(i, rxa, rxb);
        for (int j=0;j<r_ny;++j) {
            double py = row_cheb.get_cheb_point_y(j, rya, ryb);
            for (int k=0;k<r_nz;++k) {
                double pz = row_cheb.get_cheb_point_z(k, rza, rzb);
                rbox_nodes_set.push_back(Point3D(px, py, pz));
            }
        }
    }

    /** build the point set for col cluster **/

    const int c_nx = get_inter_num(cbox_xlen, lambda_ex);
    const int c_ny = get_inter_num(cbox_ylen, lambda_ex);
    const int c_nz = get_inter_num(cbox_zlen, lambda_ex);
//    cout << c_nx << ", " << c_ny << ", " << c_nz << endl;
    vector<Point3D> cbox_nodes_set;
    const int cbox_nodes_num(c_nx*c_ny*c_nz);
    ChebItp3D col_cheb(c_nx, c_ny, c_nz);
    for (int i=0;i<c_nx;++i) {
        double px = col_cheb.get_cheb_point_x(i, cxa, cxb);
        for (int j=0;j<c_ny;++j) {
            double py = col_cheb.get_cheb_point_y(j, cya, cyb);
            for (int k=0;k<c_nz;++k) {
                double pz = col_cheb.get_cheb_point_z(k, cza, czb);
                cbox_nodes_set.push_back(Point3D(px, py, pz));
            }
        }
    }

    /** **************************************** **/
    /** interpolation points generation finished **/
    /** **************************************** **/

    vector<int> row_idx_ex; /// 我们有row & col Flags
    vector<int> col_idx_ex; /// So, 不需要额外加变量了? No.

    /** parameters for ACA **/
    const char norm = 'F';
    int nrowclu(rbox_nodes_num);
    int ncolclu(cbox_nodes_num);
    int tmpRank_ex(0);
    double frbNormZ, frbNormZ_2, frbNormU, frbNormV;

    /** record the current row and col index **/
    int curRowNum, curColNum;

    /**  allocate temp space  **/
    vector<bool> flagRow_ex(nrowclu, true);
    vector<bool> flagCol_ex(ncolclu, true);
    vector<DComplex> resRowVal_ex(ncolclu, 0);
    vector<DComplex> resColVal_ex(nrowclu, 0);

    DComplex *table_ex = new DComplex [nrowclu*ncolclu];

    /** the auxiliary rank-k matrix**/
    /** Attention: it is for the Green's Function matrix, not the final rank-k hmatrix **/
    const prkmatrix rkmat_ex = new_rkmatrix(rbox_nodes_num, cbox_nodes_num, 0);

    /**Start from the 0th row**/
    curRowNum = 0;
    flagRow_ex[curRowNum] = false;
    #pragma omp parallel for
    for (int j=0;j<ncolclu;j++) {
        table_ex[curRowNum+j*nrowclu] = resRowVal_ex[j] = Green0(rbox_nodes_set[curRowNum], cbox_nodes_set[j], wavenum_ex);
    }

    curColNum = get_next_idx_aca(resRowVal_ex, flagCol_ex);
    if (curColNum == -1) return ;
    #pragma omp parallel for
    for (int i=0;i<nrowclu;i++) {
        table_ex[i+curColNum*nrowclu] = resColVal_ex[i] = Green0(rbox_nodes_set[i], cbox_nodes_set[curColNum], wavenum_ex);
    }
    if ( abs( resRowVal_ex[curColNum] ) < 1e-20 )  return ; // too close to zero is regarded as zero matrix
    DComplex temp_max_ex(resRowVal_ex[curColNum]);

    row_idx_ex.push_back(curRowNum);   // record which row is selected in ACA
    col_idx_ex.push_back(curColNum);   // record which col is selected in ACA

    rkmat_ex->V.a = new DComplex [ncolclu];
    rkmat_ex->U.a = new DComplex [nrowclu];

    for (int j=0;j<ncolclu;j++) rkmat_ex->V.a[j] = resRowVal_ex[j]/temp_max_ex;
    for (int i=0;i<nrowclu;i++) rkmat_ex->U.a[i] = resColVal_ex[i];

    frbNormV = LAPACKE_zlange(LAPACK_COL_MAJOR, norm, ncolclu, 1, (MKL_Complex16 *)rkmat_ex->V.a, ncolclu);
    frbNormU = LAPACKE_zlange(LAPACK_COL_MAJOR, norm, nrowclu, 1, (MKL_Complex16 *)rkmat_ex->U.a, nrowclu);
    frbNormZ = frbNormU * frbNormV;
    frbNormZ_2 = frbNormZ * frbNormZ;

    while ( frbNormU*frbNormV > aca_threshold*frbNormZ ) {

        tmpRank_ex++;

        if((curRowNum = get_next_idx_aca(resColVal_ex, flagRow_ex)) == -1) {
            tmpRank_ex--;
            break;
        }

        #pragma omp parallel for
        for (int j=0;j<ncolclu;j++) {
            resRowVal_ex[j] = Green0(rbox_nodes_set[curRowNum], cbox_nodes_set[j], wavenum_ex);
            table_ex[curRowNum+j*nrowclu] = resRowVal_ex[j];
            for (int i=0;i<tmpRank_ex;i++) {
                resRowVal_ex[j] -= rkmat_ex->U.a[curRowNum+i*nrowclu]*rkmat_ex->V.a[j+i*ncolclu];
            }
        }

        if ((curColNum = get_next_idx_aca(resRowVal_ex, flagCol_ex)) == -1) {
            tmpRank_ex--;
            break;
        }

        #pragma omp parallel for
        for (int i=0;i<nrowclu;i++) {
            resColVal_ex[i] = Green0(rbox_nodes_set[i], cbox_nodes_set[curColNum], wavenum_ex);
            table_ex[i+curColNum*nrowclu] = resColVal_ex[i];
            for (int j=0;j<tmpRank_ex;j++) {
                resColVal_ex[i] -= rkmat_ex->U.a[i+j*nrowclu] * rkmat_ex->V.a[curColNum+j*ncolclu];
            }
        }

        if ( abs( resRowVal_ex[curColNum] ) < 1e-20 ) {
            tmpRank_ex--;
            break;
        }
        temp_max_ex = resRowVal_ex[curColNum];

        rkmat_ex->V.a = (DComplex *) realloc(rkmat_ex->V.a, (tmpRank_ex+1)*ncolclu*sizeof(DComplex));
        rkmat_ex->U.a = (DComplex *) realloc(rkmat_ex->U.a, (tmpRank_ex+1)*nrowclu*sizeof(DComplex));

        row_idx_ex.push_back(curRowNum);
        col_idx_ex.push_back(curColNum);

        #pragma omp parallel for
        for (int j=0;j<ncolclu;j++) rkmat_ex->V.a[j+tmpRank_ex*ncolclu] = resRowVal_ex[j]/temp_max_ex;
        #pragma omp parallel for
        for (int i=0;i<nrowclu;i++) rkmat_ex->U.a[i+tmpRank_ex*nrowclu] = resColVal_ex[i];

        frbNormV = LAPACKE_zlange(LAPACK_COL_MAJOR, norm, ncolclu, 1, (MKL_Complex16 *)rkmat_ex->V.a+tmpRank_ex*ncolclu, ncolclu);
        frbNormU = LAPACKE_zlange(LAPACK_COL_MAJOR, norm, nrowclu, 1, (MKL_Complex16 *)rkmat_ex->U.a+tmpRank_ex*nrowclu, nrowclu);
        for (int i=0;i<tmpRank_ex;i++) {
            DComplex c1, c2;
            cblas_zdotc_sub(nrowclu, rkmat_ex->U.a+i*nrowclu, 1, rkmat_ex->U.a+tmpRank_ex*nrowclu, 1, &c1);
            cblas_zdotc_sub(ncolclu, rkmat_ex->V.a+i*ncolclu, 1, rkmat_ex->V.a+tmpRank_ex*ncolclu, 1, &c2);
            frbNormZ_2 += 2*abs( c1 )*abs( c2 );
        }
        frbNormZ_2 += frbNormU*frbNormU*frbNormV*frbNormV;
        frbNormZ = sqrt(frbNormZ_2);
    }

    tmpRank_ex ++;

    cout << "m=" << nrowclu << ", " << "n=" << ncolclu << ": " << tmpRank_ex << endl;

    /** **************************************** **/
    /** interpolation points generation finished **/
    /** **************************************** **/

    vector<int> row_idx_in; /// 我们有row & col Flags
    vector<int> col_idx_in; /// So, 不需要额外加变量了? No.

    /** parameters for ACA **/


    int tmpRank_in(0);
    frbNormZ = frbNormZ_2 = frbNormU = frbNormV = 0;

    /**  allocate temp space  **/
    vector<bool> flagRow_in(nrowclu, true);
    vector<bool> flagCol_in(ncolclu, true);
    vector<DComplex> resRowVal_in(ncolclu, 0);
    vector<DComplex> resColVal_in(nrowclu, 0);

    DComplex *table_in = new DComplex [nrowclu*ncolclu];

    /** the auxiliary rank-k matrix**/
    /** Attention: it is for the Green's Function matrix, not the final rank-k hmatrix **/
    const prkmatrix rkmat_in = new_rkmatrix(rbox_nodes_num, cbox_nodes_num, 0);

    /**Start from the 0th row**/
    curRowNum = 0;
    flagRow_in[curRowNum] = false;
    #pragma omp parallel for
    for (int j=0;j<ncolclu;j++) {
        table_in[curRowNum+j*nrowclu] = resRowVal_in[j] = Green0(rbox_nodes_set[curRowNum], cbox_nodes_set[j], wavenum_in);
    }

    curColNum = get_next_idx_aca(resRowVal_in, flagCol_in);
    if (curColNum == -1) return ;
    #pragma omp parallel for
    for (int i=0;i<nrowclu;i++) {
        table_in[i+curColNum*nrowclu] = resColVal_in[i] = Green0(rbox_nodes_set[i], cbox_nodes_set[curColNum], wavenum_in);
    }
    if ( abs( resRowVal_in[curColNum] ) < 1e-20 )  return ; // too close to zero is regarded as zero matrix
    DComplex temp_max_in(resRowVal_in[curColNum]);

    row_idx_in.push_back(curRowNum);   // record which row is selected in ACA
    col_idx_in.push_back(curColNum);   // record which col is selected in ACA

    rkmat_in->V.a = new DComplex [ncolclu];
    rkmat_in->U.a = new DComplex [nrowclu];

    for (int j=0;j<ncolclu;j++) rkmat_in->V.a[j] = resRowVal_in[j]/temp_max_in;
    for (int i=0;i<nrowclu;i++) rkmat_in->U.a[i] = resColVal_in[i];

    frbNormV = LAPACKE_zlange(LAPACK_COL_MAJOR, norm, ncolclu, 1, (MKL_Complex16 *)rkmat_in->V.a, ncolclu);
    frbNormU = LAPACKE_zlange(LAPACK_COL_MAJOR, norm, nrowclu, 1, (MKL_Complex16 *)rkmat_in->U.a, nrowclu);
    frbNormZ = frbNormU * frbNormV;
    frbNormZ_2 = frbNormZ * frbNormZ;

    while ( frbNormU*frbNormV > aca_threshold*frbNormZ ) {

        tmpRank_in++;

        if((curRowNum = get_next_idx_aca(resColVal_in, flagRow_in)) == -1) {
            tmpRank_in--;
            break;
        }

        #pragma omp parallel for
        for (int j=0;j<ncolclu;j++) {
            resRowVal_in[j] = Green0(rbox_nodes_set[curRowNum], cbox_nodes_set[j], wavenum_in);
            table_in[curRowNum+j*nrowclu] = resRowVal_in[j];
            for (int i=0;i<tmpRank_in;i++) {
                resRowVal_in[j] -= rkmat_in->U.a[curRowNum+i*nrowclu]*rkmat_in->V.a[j+i*ncolclu];
            }
        }

        if ((curColNum = get_next_idx_aca(resRowVal_in, flagCol_in)) == -1) {
            tmpRank_in--;
            break;
        }

        #pragma omp parallel for
        for (int i=0;i<nrowclu;i++) {
            resColVal_in[i] = Green0(rbox_nodes_set[i], cbox_nodes_set[curColNum], wavenum_in);
            table_in[i+curColNum*nrowclu] = resColVal_in[i];
            for (int j=0;j<tmpRank_in;j++) {
                resColVal_in[i] -= rkmat_in->U.a[i+j*nrowclu] * rkmat_in->V.a[curColNum+j*ncolclu];
            }
        }

        if ( abs( resRowVal_in[curColNum] ) < 1e-20 ) {
            tmpRank_in--;
            break;
        }
        temp_max_in = resRowVal_in[curColNum];

        rkmat_in->V.a = (DComplex *) realloc(rkmat_in->V.a, (tmpRank_in+1)*ncolclu*sizeof(DComplex));
        rkmat_in->U.a = (DComplex *) realloc(rkmat_in->U.a, (tmpRank_in+1)*nrowclu*sizeof(DComplex));

        row_idx_in.push_back(curRowNum);
        col_idx_in.push_back(curColNum);

        #pragma omp parallel for
        for (int j=0;j<ncolclu;j++) rkmat_in->V.a[j+tmpRank_in*ncolclu] = resRowVal_in[j]/temp_max_in;
        #pragma omp parallel for
        for (int i=0;i<nrowclu;i++) rkmat_in->U.a[i+tmpRank_in*nrowclu] = resColVal_in[i];

        frbNormV = LAPACKE_zlange(LAPACK_COL_MAJOR, norm, ncolclu, 1, (MKL_Complex16 *)rkmat_in->V.a+tmpRank_in*ncolclu, ncolclu);
        frbNormU = LAPACKE_zlange(LAPACK_COL_MAJOR, norm, nrowclu, 1, (MKL_Complex16 *)rkmat_in->U.a+tmpRank_in*nrowclu, nrowclu);
        for (int i=0;i<tmpRank_in;i++) {
            DComplex c1, c2;
            cblas_zdotc_sub(nrowclu, rkmat_in->U.a+i*nrowclu, 1, rkmat_in->U.a+tmpRank_in*nrowclu, 1, &c1);
            cblas_zdotc_sub(ncolclu, rkmat_in->V.a+i*ncolclu, 1, rkmat_in->V.a+tmpRank_in*ncolclu, 1, &c2);
            frbNormZ_2 += 2*abs( c1 )*abs( c2 );
        }
        frbNormZ_2 += frbNormU*frbNormU*frbNormV*frbNormV;
        frbNormZ = sqrt(frbNormZ_2);
    }

    tmpRank_in ++;

    cout << "m=" << nrowclu << ", " << "n=" << ncolclu << ": " << tmpRank_in << endl;

    nrowclu = hm0->rc->size;
    ncolclu = hm0->cc->size;

    /** actually, the table can be avoided **/
    DComplex *gs_mat_ex = new DComplex [tmpRank_ex*tmpRank_ex];
    #pragma omp parallel for
    for (int i=0;i<tmpRank_ex;++i) {
        for (int j=0;j<tmpRank_ex;++j) {
            gs_mat_ex[i+j*tmpRank_ex] = table_ex[row_idx_ex[i]+col_idx_ex[j]*rbox_nodes_num];
        }
    }/// This S matrix is still low rank, so direct inverse is not feasible.
    delete [] table_ex;

    DComplex *gs_mat_in = new DComplex [tmpRank_in*tmpRank_in];
    #pragma omp parallel for
    for (int i=0;i<tmpRank_in;++i) {
        for (int j=0;j<tmpRank_in;++j) {
            gs_mat_in[i+j*tmpRank_in] = table_in[row_idx_in[i]+col_idx_in[j]*rbox_nodes_num];
        }
    }/// This S matrix is still low rank, so direct inverse is not feasible.
    delete [] table_in;

    matInv_z(gs_mat_ex, tmpRank_ex);
    matInv_z(gs_mat_in, tmpRank_in);

    Complex3D *newa_vec_e_ex = new Complex3D [tmpRank_ex*nrowclu];
    Complex3D *newb_vec_e_ex = new Complex3D [tmpRank_ex*ncolclu];
    prkmatrix rk_vec_x_e_ex = new_rkmatrix(nrowclu, ncolclu, tmpRank_ex);
    prkmatrix rk_vec_y_e_ex = new_rkmatrix(nrowclu, ncolclu, tmpRank_ex);
    prkmatrix rk_vec_z_e_ex = new_rkmatrix(nrowclu, ncolclu, tmpRank_ex);
    prkmatrix rk_sca_e_ex = new_rkmatrix(nrowclu, ncolclu, tmpRank_ex);
    prkmatrix tmp_rk_e_ex = new_rkmatrix(nrowclu, ncolclu, 0);

    Complex3D *newa_vec_e_in = new Complex3D [tmpRank_in*nrowclu];
    Complex3D *newb_vec_e_in = new Complex3D [tmpRank_in*ncolclu];
    prkmatrix rk_vec_x_e_in = new_rkmatrix(nrowclu, ncolclu, tmpRank_in);
    prkmatrix rk_vec_y_e_in = new_rkmatrix(nrowclu, ncolclu, tmpRank_in);
    prkmatrix rk_vec_z_e_in = new_rkmatrix(nrowclu, ncolclu, tmpRank_in);
    prkmatrix rk_sca_e_in = new_rkmatrix(nrowclu, ncolclu, tmpRank_in);
    prkmatrix tmp_rk_e_in = new_rkmatrix(nrowclu, ncolclu, 0);

    Complex3D *newa_vec_m_ex = new Complex3D [tmpRank_ex*nrowclu];
    Complex3D *newb_vec_m_ex = new Complex3D [tmpRank_ex*ncolclu];
    prkmatrix rk_vec_x_m_ex = new_rkmatrix(nrowclu, ncolclu, tmpRank_ex);
    prkmatrix rk_vec_y_m_ex = new_rkmatrix(nrowclu, ncolclu, tmpRank_ex);
    prkmatrix rk_vec_z_m_ex = new_rkmatrix(nrowclu, ncolclu, tmpRank_ex);
    prkmatrix tmp_rk_m_ex = new_rkmatrix(nrowclu, ncolclu, 0);

    Complex3D *newa_vec_m_in = new Complex3D [tmpRank_in*nrowclu];
    Complex3D *newb_vec_m_in = new Complex3D [tmpRank_in*ncolclu];
    prkmatrix rk_vec_x_m_in = new_rkmatrix(nrowclu, ncolclu, tmpRank_in);
    prkmatrix rk_vec_y_m_in = new_rkmatrix(nrowclu, ncolclu, tmpRank_in);
    prkmatrix rk_vec_z_m_in = new_rkmatrix(nrowclu, ncolclu, tmpRank_in);
    prkmatrix tmp_rk_m_in = new_rkmatrix(nrowclu, ncolclu, 0);

    /// For the external, build the elementary info
    #pragma omp parallel for
    for (int j=0;j<tmpRank_ex;++j) {
        for (int i=0;i<nrowclu;++i) {
            intg.single_layer_intg(hm0->rc->idx[i], cbox_nodes_set[col_idx_ex[j]],
                                newa_vec_e_ex[i+nrowclu*j], rk_sca_e_ex->U.a[i+nrowclu*j], true);
            newa_vec_m_ex[i+nrowclu*j] = intg.single_layer_intg_m3(hm0->rc->idx[i], cbox_nodes_set[col_idx_ex[j]], true);
        }
        for (int i=0;i<ncolclu;++i) {
            intg.single_layer_intg(hm0->cc->idx[i], rbox_nodes_set[row_idx_ex[j]],
                                newb_vec_e_ex[i+ncolclu*j], rk_sca_e_ex->V.a[i+ncolclu*j], true);
            newb_vec_m_ex[i+ncolclu*j] = intg.single_layer_intg_m4(hm0->cc->idx[i], rbox_nodes_set[row_idx_ex[j]], true);
        }
    }

    /// For the internal, build the elementary info
    #pragma omp parallel for
    for (int j=0;j<tmpRank_in;++j) {
        for (int i=0;i<nrowclu;++i) {
            intg.single_layer_intg(hm0->rc->idx[i], cbox_nodes_set[col_idx_in[j]],
                                newa_vec_e_in[i+nrowclu*j], rk_sca_e_in->U.a[i+nrowclu*j], false);
            newa_vec_m_in[i+nrowclu*j] = intg.single_layer_intg_m3(hm0->rc->idx[i], cbox_nodes_set[col_idx_in[j]], false);
        }
        for (int i=0;i<ncolclu;++i) {
            intg.single_layer_intg(hm0->cc->idx[i], rbox_nodes_set[row_idx_in[j]],
                                newb_vec_e_in[i+ncolclu*j], rk_sca_e_in->V.a[i+ncolclu*j], false);
            newb_vec_m_in[i+ncolclu*j] = intg.single_layer_intg_m4(hm0->cc->idx[i], rbox_nodes_set[row_idx_in[j]], false);
        }
    }

    /// move the values into rk matrices, for external
    #pragma omp parallel for
    for (int j=0;j<tmpRank_ex;++j) {
        for (int i=0;i<nrowclu;++i) {
            rk_vec_x_e_ex->U.a[i+nrowclu*j] = coefA_ex*newa_vec_e_ex[i+nrowclu*j][0];
            rk_vec_y_e_ex->U.a[i+nrowclu*j] = coefA_ex*newa_vec_e_ex[i+nrowclu*j][1];
            rk_vec_z_e_ex->U.a[i+nrowclu*j] = coefA_ex*newa_vec_e_ex[i+nrowclu*j][2];
            rk_sca_e_ex->U.a[i+nrowclu*j] *= coefV_ex;
            rk_vec_x_m_ex->U.a[i+nrowclu*j] = newa_vec_m_ex[i+nrowclu*j][0];
            rk_vec_y_m_ex->U.a[i+nrowclu*j] = newa_vec_m_ex[i+nrowclu*j][1];
            rk_vec_z_m_ex->U.a[i+nrowclu*j] = newa_vec_m_ex[i+nrowclu*j][2];
        }
        for (int i=0;i<ncolclu;++i) {
            rk_vec_x_e_ex->V.a[i+ncolclu*j] = conj(newb_vec_e_ex[i+ncolclu*j][0]);
            rk_vec_y_e_ex->V.a[i+ncolclu*j] = conj(newb_vec_e_ex[i+ncolclu*j][1]);
            rk_vec_z_e_ex->V.a[i+ncolclu*j] = conj(newb_vec_e_ex[i+ncolclu*j][2]);
            rk_sca_e_ex->V.a[i+ncolclu*j] = conj(rk_sca_e_ex->V.a[i+ncolclu*j]);
            rk_vec_x_m_ex->V.a[i+ncolclu*j] = conj(newb_vec_m_ex[i+ncolclu*j][0]);
            rk_vec_y_m_ex->V.a[i+ncolclu*j] = conj(newb_vec_m_ex[i+ncolclu*j][1]);
            rk_vec_z_m_ex->V.a[i+ncolclu*j] = conj(newb_vec_m_ex[i+ncolclu*j][2]);
        }
    }

    /// move the values into rk matrices, for external
    #pragma omp parallel for
    for (int j=0;j<tmpRank_in;++j) {
        for (int i=0;i<nrowclu;++i) {
            rk_vec_x_e_in->U.a[i+nrowclu*j] = coefA_in*newa_vec_e_in[i+nrowclu*j][0];
            rk_vec_y_e_in->U.a[i+nrowclu*j] = coefA_in*newa_vec_e_in[i+nrowclu*j][1];
            rk_vec_z_e_in->U.a[i+nrowclu*j] = coefA_in*newa_vec_e_in[i+nrowclu*j][2];
            rk_sca_e_in->U.a[i+nrowclu*j] *= coefV_in;
            rk_vec_x_m_in->U.a[i+nrowclu*j] = newa_vec_m_in[i+nrowclu*j][0];
            rk_vec_y_m_in->U.a[i+nrowclu*j] = newa_vec_m_in[i+nrowclu*j][1];
            rk_vec_z_m_in->U.a[i+nrowclu*j] = newa_vec_m_in[i+nrowclu*j][2];
        }
        for (int i=0;i<ncolclu;++i) {
            rk_vec_x_e_in->V.a[i+ncolclu*j] = conj(newb_vec_e_in[i+ncolclu*j][0]);
            rk_vec_y_e_in->V.a[i+ncolclu*j] = conj(newb_vec_e_in[i+ncolclu*j][1]);
            rk_vec_z_e_in->V.a[i+ncolclu*j] = conj(newb_vec_e_in[i+ncolclu*j][2]);
            rk_sca_e_in->V.a[i+ncolclu*j] = conj(rk_sca_e_in->V.a[i+ncolclu*j]);
            rk_vec_x_m_in->V.a[i+ncolclu*j] = conj(newb_vec_m_in[i+ncolclu*j][0]);
            rk_vec_y_m_in->V.a[i+ncolclu*j] = conj(newb_vec_m_in[i+ncolclu*j][1]);
            rk_vec_z_m_in->V.a[i+ncolclu*j] = conj(newb_vec_m_in[i+ncolclu*j][2]);
        }
    }

    /// L operator external
    prkmatrix rk_vec_nx_e_ex = new_rkmatrix(nrowclu, ncolclu, tmpRank_ex);
    prkmatrix rk_vec_ny_e_ex = new_rkmatrix(nrowclu, ncolclu, tmpRank_ex);
    prkmatrix rk_vec_nz_e_ex = new_rkmatrix(nrowclu, ncolclu, tmpRank_ex);
    prkmatrix rk_sca_n_e_ex = new_rkmatrix(nrowclu, ncolclu, tmpRank_ex);
    copy_rkmatrix(false, rk_vec_x_e_ex, rk_vec_nx_e_ex);
    copy_rkmatrix(false, rk_vec_y_e_ex, rk_vec_ny_e_ex);
    copy_rkmatrix(false, rk_vec_z_e_ex, rk_vec_nz_e_ex);
    copy_rkmatrix(false, rk_sca_e_ex,   rk_sca_n_e_ex);


    if (nrowclu<=ncolclu) {
        mmp_z(rk_vec_x_e_ex->U.a, CblasNoTrans, gs_mat_ex, CblasNoTrans, rk_vec_nx_e_ex->U.a, nrowclu, tmpRank_ex, tmpRank_ex);
        mmp_z(rk_vec_y_e_ex->U.a, CblasNoTrans, gs_mat_ex, CblasNoTrans, rk_vec_ny_e_ex->U.a, nrowclu, tmpRank_ex, tmpRank_ex);
        mmp_z(rk_vec_z_e_ex->U.a, CblasNoTrans, gs_mat_ex, CblasNoTrans, rk_vec_nz_e_ex->U.a, nrowclu, tmpRank_ex, tmpRank_ex);
        mmp_z(rk_sca_e_ex->U.a,   CblasNoTrans, gs_mat_ex, CblasNoTrans, rk_sca_n_e_ex->U.a,  nrowclu, tmpRank_ex, tmpRank_ex);
    }
    else {
        mmp_z(rk_vec_x_e_ex->V.a, CblasNoTrans, gs_mat_ex, CblasConjTrans, rk_vec_nx_e_ex->V.a, ncolclu, tmpRank_ex, tmpRank_ex);
        mmp_z(rk_vec_y_e_ex->V.a, CblasNoTrans, gs_mat_ex, CblasConjTrans, rk_vec_ny_e_ex->V.a, ncolclu, tmpRank_ex, tmpRank_ex);
        mmp_z(rk_vec_z_e_ex->V.a, CblasNoTrans, gs_mat_ex, CblasConjTrans, rk_vec_nz_e_ex->V.a, ncolclu, tmpRank_ex, tmpRank_ex);
        mmp_z(rk_sca_e_ex->V.a,   CblasNoTrans, gs_mat_ex, CblasConjTrans, rk_sca_n_e_ex->V.a,  ncolclu, tmpRank_ex, tmpRank_ex);
    }

    copy_rkmatrix(false, rk_vec_nx_e_ex, tmp_rk_e_ex);
    ptruncmode tm = new_truncmode();
    add_rkmatrix(1.0, rk_vec_ny_e_ex, tm, svd_threshold, tmp_rk_e_ex);
    add_rkmatrix(1.0, rk_vec_nz_e_ex, tm, svd_threshold, tmp_rk_e_ex);
    add_rkmatrix(1.0, rk_sca_n_e_ex,  tm, svd_threshold, tmp_rk_e_ex);
    del_rkmatrix(rk_vec_nx_e_ex);
    del_rkmatrix(rk_vec_ny_e_ex);
    del_rkmatrix(rk_vec_nz_e_ex);
    del_rkmatrix(rk_sca_n_e_ex);

    /// L operator internal
    prkmatrix rk_vec_nx_e_in = new_rkmatrix(nrowclu, ncolclu, tmpRank_in);
    prkmatrix rk_vec_ny_e_in = new_rkmatrix(nrowclu, ncolclu, tmpRank_in);
    prkmatrix rk_vec_nz_e_in = new_rkmatrix(nrowclu, ncolclu, tmpRank_in);
    prkmatrix rk_sca_n_e_in = new_rkmatrix(nrowclu, ncolclu, tmpRank_in);
    copy_rkmatrix(false, rk_vec_x_e_in, rk_vec_nx_e_in);
    copy_rkmatrix(false, rk_vec_y_e_in, rk_vec_ny_e_in);
    copy_rkmatrix(false, rk_vec_z_e_in, rk_vec_nz_e_in);
    copy_rkmatrix(false, rk_sca_e_in,   rk_sca_n_e_in);

    if (nrowclu<=ncolclu) {
        mmp_z(rk_vec_x_e_in->U.a, CblasNoTrans, gs_mat_in, CblasNoTrans, rk_vec_nx_e_in->U.a, nrowclu, tmpRank_in, tmpRank_in);
        mmp_z(rk_vec_y_e_in->U.a, CblasNoTrans, gs_mat_in, CblasNoTrans, rk_vec_ny_e_in->U.a, nrowclu, tmpRank_in, tmpRank_in);
        mmp_z(rk_vec_z_e_in->U.a, CblasNoTrans, gs_mat_in, CblasNoTrans, rk_vec_nz_e_in->U.a, nrowclu, tmpRank_in, tmpRank_in);
        mmp_z(rk_sca_e_in->U.a,   CblasNoTrans, gs_mat_in, CblasNoTrans, rk_sca_n_e_in->U.a,  nrowclu, tmpRank_in, tmpRank_in);
    }
    else {
        mmp_z(rk_vec_x_e_in->V.a, CblasNoTrans, gs_mat_in, CblasConjTrans, rk_vec_nx_e_in->V.a, ncolclu, tmpRank_in, tmpRank_in);
        mmp_z(rk_vec_y_e_in->V.a, CblasNoTrans, gs_mat_in, CblasConjTrans, rk_vec_ny_e_in->V.a, ncolclu, tmpRank_in, tmpRank_in);
        mmp_z(rk_vec_z_e_in->V.a, CblasNoTrans, gs_mat_in, CblasConjTrans, rk_vec_nz_e_in->V.a, ncolclu, tmpRank_in, tmpRank_in);
        mmp_z(rk_sca_e_in->V.a,   CblasNoTrans, gs_mat_in, CblasConjTrans, rk_sca_n_e_in->V.a,  ncolclu, tmpRank_in, tmpRank_in);
    }

    copy_rkmatrix(false, rk_vec_nx_e_in, tmp_rk_e_in);
    add_rkmatrix(1.0, rk_vec_ny_e_in, tm, svd_threshold, tmp_rk_e_in);
    add_rkmatrix(1.0, rk_vec_nz_e_in, tm, svd_threshold, tmp_rk_e_in);
    add_rkmatrix(1.0, rk_sca_n_e_in,  tm, svd_threshold, tmp_rk_e_in);
    del_rkmatrix(rk_vec_nx_e_in);
    del_rkmatrix(rk_vec_ny_e_in);
    del_rkmatrix(rk_vec_nz_e_in);
    del_rkmatrix(rk_sca_n_e_in);

    /// K operator external
    prkmatrix rk_vec_nx_m_ex = new_rkmatrix(nrowclu, ncolclu, tmpRank_ex);
    prkmatrix rk_vec_ny_m_ex = new_rkmatrix(nrowclu, ncolclu, tmpRank_ex);
    prkmatrix rk_vec_nz_m_ex = new_rkmatrix(nrowclu, ncolclu, tmpRank_ex);
    copy_rkmatrix(false, rk_vec_x_m_ex, rk_vec_nx_m_ex);
    copy_rkmatrix(false, rk_vec_y_m_ex, rk_vec_ny_m_ex);
    copy_rkmatrix(false, rk_vec_z_m_ex, rk_vec_nz_m_ex);

    if (nrowclu<=ncolclu) {
        mmp_z(rk_vec_x_m_ex->U.a, CblasNoTrans, gs_mat_ex, CblasNoTrans, rk_vec_nx_m_ex->U.a, nrowclu, tmpRank_ex, tmpRank_ex);
        mmp_z(rk_vec_y_m_ex->U.a, CblasNoTrans, gs_mat_ex, CblasNoTrans, rk_vec_ny_m_ex->U.a, nrowclu, tmpRank_ex, tmpRank_ex);
        mmp_z(rk_vec_z_m_ex->U.a, CblasNoTrans, gs_mat_ex, CblasNoTrans, rk_vec_nz_m_ex->U.a, nrowclu, tmpRank_ex, tmpRank_ex);
    }
    else {
        mmp_z(rk_vec_x_m_ex->V.a, CblasNoTrans, gs_mat_ex, CblasConjTrans, rk_vec_nx_m_ex->V.a, ncolclu, tmpRank_ex, tmpRank_ex);
        mmp_z(rk_vec_y_m_ex->V.a, CblasNoTrans, gs_mat_ex, CblasConjTrans, rk_vec_ny_m_ex->V.a, ncolclu, tmpRank_ex, tmpRank_ex);
        mmp_z(rk_vec_z_m_ex->V.a, CblasNoTrans, gs_mat_ex, CblasConjTrans, rk_vec_nz_m_ex->V.a, ncolclu, tmpRank_ex, tmpRank_ex);
    }

    copy_rkmatrix(false, rk_vec_nx_m_ex, tmp_rk_m_ex);
    add_rkmatrix(1.0, rk_vec_ny_m_ex, tm, svd_threshold, tmp_rk_m_ex);
    add_rkmatrix(1.0, rk_vec_nz_m_ex, tm, svd_threshold, tmp_rk_m_ex);
    del_rkmatrix(rk_vec_nx_m_ex);
    del_rkmatrix(rk_vec_ny_m_ex);
    del_rkmatrix(rk_vec_nz_m_ex);

    /// K operator internal
    prkmatrix rk_vec_nx_m_in = new_rkmatrix(nrowclu, ncolclu, tmpRank_in);
    prkmatrix rk_vec_ny_m_in = new_rkmatrix(nrowclu, ncolclu, tmpRank_in);
    prkmatrix rk_vec_nz_m_in = new_rkmatrix(nrowclu, ncolclu, tmpRank_in);
    copy_rkmatrix(false, rk_vec_x_m_in, rk_vec_nx_m_in);
    copy_rkmatrix(false, rk_vec_y_m_in, rk_vec_ny_m_in);
    copy_rkmatrix(false, rk_vec_z_m_in, rk_vec_nz_m_in);

    if (nrowclu<=ncolclu) {
        mmp_z(rk_vec_x_m_in->U.a, CblasNoTrans, gs_mat_in, CblasNoTrans, rk_vec_nx_m_in->U.a, nrowclu, tmpRank_in, tmpRank_in);
        mmp_z(rk_vec_y_m_in->U.a, CblasNoTrans, gs_mat_in, CblasNoTrans, rk_vec_ny_m_in->U.a, nrowclu, tmpRank_in, tmpRank_in);
        mmp_z(rk_vec_z_m_in->U.a, CblasNoTrans, gs_mat_in, CblasNoTrans, rk_vec_nz_m_in->U.a, nrowclu, tmpRank_in, tmpRank_in);
    }
    else {
        mmp_z(rk_vec_x_m_in->V.a, CblasNoTrans, gs_mat_in, CblasConjTrans, rk_vec_nx_m_in->V.a, ncolclu, tmpRank_in, tmpRank_in);
        mmp_z(rk_vec_y_m_in->V.a, CblasNoTrans, gs_mat_in, CblasConjTrans, rk_vec_ny_m_in->V.a, ncolclu, tmpRank_in, tmpRank_in);
        mmp_z(rk_vec_z_m_in->V.a, CblasNoTrans, gs_mat_in, CblasConjTrans, rk_vec_nz_m_in->V.a, ncolclu, tmpRank_in, tmpRank_in);
    }

    copy_rkmatrix(false, rk_vec_nx_m_in, tmp_rk_m_in);
    add_rkmatrix(1.0, rk_vec_ny_m_in, tm, svd_threshold, tmp_rk_m_in);
    add_rkmatrix(1.0, rk_vec_nz_m_in, tm, svd_threshold, tmp_rk_m_in);
    del_rkmatrix(rk_vec_nx_m_in);
    del_rkmatrix(rk_vec_ny_m_in);
    del_rkmatrix(rk_vec_nz_m_in);

    /*****************************************************************/
    /**                        assembly rkmat                       **/
    /*****************************************************************/

    /// for hm0
    copy_rkmatrix(false, tmp_rk_e_ex, hm0->r);
    add_rkmatrix(1.0, tmp_rk_e_in, tm, svd_threshold, hm0->r);

    cout << hm0->r->k << endl;

    /// for hm1
    copy_rkmatrix(false, tmp_rk_m_ex, hm1->r);
    add_rkmatrix(1.0, tmp_rk_m_in, tm, svd_threshold, hm1->r);
//    cout << hm2->r->k << endl;

    /// for hm2
    copy_rkmatrix(false, hm1->r, hm2->r);
    for (int i=0;i<hm2->r->U.cols*hm2->r->U.rows;++i)   hm2->r->U.a[i] *= -1.0;
    cout << hm1->r->k << endl;

    /// for hm3
    const DComplex eta_2(ETA0*ETA0);
    const DComplex coef_ex = epsilon_ex / eta_2;
    const DComplex coef_in = epsilon_in / eta_2;
    copy_rkmatrix(false, tmp_rk_e_ex, hm3->r);
    for (int i=0;i<hm3->r->U.cols*hm3->r->U.rows;++i)   hm3->r->U.a[i] *= coef_ex;
    add_rkmatrix(coef_in, tmp_rk_e_in, tm, svd_threshold, hm3->r);
    cout << hm3->r->k << endl;

    /*****************************************************************/
    /**                        assembly rkmat                       **/
    /*****************************************************************/

    delete [] gs_mat_ex;
    delete [] gs_mat_in;

    delete [] newa_vec_e_ex;
    delete [] newb_vec_e_ex;
    del_rkmatrix(rk_vec_x_e_ex);
    del_rkmatrix(rk_vec_y_e_ex);
    del_rkmatrix(rk_vec_z_e_ex);
    del_rkmatrix(rk_sca_e_ex);
    del_rkmatrix(tmp_rk_e_ex);

    delete [] newa_vec_e_in;
    delete [] newb_vec_e_in;
    del_rkmatrix(rk_vec_x_e_in);
    del_rkmatrix(rk_vec_y_e_in);
    del_rkmatrix(rk_vec_z_e_in);
    del_rkmatrix(rk_sca_e_in);
    del_rkmatrix(tmp_rk_e_in);

    delete [] newa_vec_m_ex;
    delete [] newb_vec_m_ex;
    del_rkmatrix(rk_vec_x_m_ex);
    del_rkmatrix(rk_vec_y_m_ex);
    del_rkmatrix(rk_vec_z_m_ex);
    del_rkmatrix(tmp_rk_m_ex);

    delete [] newa_vec_m_in;
    delete [] newb_vec_m_in;
    del_rkmatrix(rk_vec_x_m_in);
    del_rkmatrix(rk_vec_y_m_in);
    del_rkmatrix(rk_vec_z_m_in);
    del_rkmatrix(tmp_rk_m_in);
}

void
FastBem::hca(phmatrix hm, const Integration &intg, const char c)
{
//    cout << hca_times++ << ": use HCA" << endl;

    /** parameter for integral **/
    const double freq(intg.get_freq());
    const double omega(2*PI*freq);
    const DComplex epsilon_r(1.0);
    const DComplex wavenum = omega*(1.0/Constants::SPEED_OF_LIGHT)*sqrt(epsilon_r);
    const double lambda(2*Constants::PI/wavenum.real());
    const DComplex coefA(CJ*omega*MU0), coefV(-CJ/(epsilon_r*EPSILON0*omega));

    /** the bounding box is crucial for hca **/
    const double *rbox_min = hm->rc->bmin;
    const double *rbox_max = hm->rc->bmax;
    const double *cbox_min = hm->cc->bmin;
    const double *cbox_max = hm->cc->bmax;

    /** fixed points number may be too small or too large
     ** use the space resolution method to select the size adaptively **/
    const double rbox_xlen = rbox_max[0] - rbox_min[0];
    const double rbox_ylen = rbox_max[1] - rbox_min[1];
    const double rbox_zlen = rbox_max[2] - rbox_min[2];
    const int r_nx = get_inter_num(rbox_xlen, lambda);
    const int r_ny = get_inter_num(rbox_ylen, lambda);
    const int r_nz = get_inter_num(rbox_zlen, lambda);
//    cout << r_nx << ", " << r_ny << ", " << r_nz << endl;

    /** build the point set for row cluster **/
    vector<Point3D> rbox_nodes_set;
    const int rbox_nodes_num(r_nx*r_ny*r_nz);
    ChebItp3D row_cheb(r_nx, r_ny, r_nz);
    double rxa = rbox_min[0], rxb = rbox_max[0];
    double rya = rbox_min[1], ryb = rbox_max[1];
    double rza = rbox_min[2], rzb = rbox_max[2];
    for (int i=0;i<r_nx;++i) {
        double px = row_cheb.get_cheb_point_x(i, rxa, rxb);
        for (int j=0;j<r_ny;++j) {
            double py = row_cheb.get_cheb_point_y(j, rya, ryb);
            for (int k=0;k<r_nz;++k) {
                double pz = row_cheb.get_cheb_point_z(k, rza, rzb);
                rbox_nodes_set.push_back(Point3D(px, py, pz));
            }
        }
    }

    /** build the point set for col cluster **/
    const double cbox_xlen = cbox_max[0] - cbox_min[0];
    const double cbox_ylen = cbox_max[1] - cbox_min[1];
    const double cbox_zlen = cbox_max[2] - cbox_min[2];
    const int c_nx = get_inter_num(cbox_xlen, lambda);
    const int c_ny = get_inter_num(cbox_ylen, lambda);
    const int c_nz = get_inter_num(cbox_zlen, lambda);
//    cout << c_nx << ", " << c_ny << ", " << c_nz << endl;

    vector<Point3D> cbox_nodes_set;
    const int cbox_nodes_num(c_nx*c_ny*c_nz);
    ChebItp3D col_cheb(c_nx, c_ny, c_nz);
    double cxa = cbox_min[0], cxb = cbox_max[0];
    double cya = cbox_min[1], cyb = cbox_max[1];
    double cza = cbox_min[2], czb = cbox_max[2];
    for (int i=0;i<c_nx;++i) {
        double px = col_cheb.get_cheb_point_x(i, cxa, cxb);
        for (int j=0;j<c_ny;++j) {
            double py = col_cheb.get_cheb_point_y(j, cya, cyb);
            for (int k=0;k<c_nz;++k) {
                double pz = col_cheb.get_cheb_point_z(k, cza, czb);
                cbox_nodes_set.push_back(Point3D(px, py, pz));
            }
        }
    }

    /** **************************************** **/
    /** interpolation points generation finished **/
    /** **************************************** **/

    vector<int> pre_row_idx; /// 我们有row & col Flags
    vector<int> pre_col_idx; /// So, 不需要额外加变量了? No.

    /** parameters for ACA **/
    const char norm = 'F';
    int nrowclu(rbox_nodes_num);
    int ncolclu(cbox_nodes_num);
    int tmpRank(0);
    double frbNormZ, frbNormZ_2, frbNormU, frbNormV;

    /** record the current row and col index **/
    int curRowNum, curColNum;

    /**  allocate temp space  **/
    vector<bool> flagRow(nrowclu, true);
    vector<bool> flagCol(ncolclu, true);
    vector<DComplex> resRowVal(ncolclu, 0);
    vector<DComplex> resColVal(nrowclu, 0);

    DComplex *table = new DComplex [nrowclu*ncolclu];

    /** the auxiliary rank-k matrix**/
    /** Attention: it is for the Green's Function matrix, not the final rank-k hmatrix **/
    const prkmatrix rkmat = new_rkmatrix(rbox_nodes_num, cbox_nodes_num, 0);

    /**Start from the 0th row**/
    curRowNum = 0;
    flagRow[curRowNum] = false;
    #pragma omp parallel for
    for (int j=0;j<ncolclu;j++) {
        table[curRowNum+j*nrowclu] = resRowVal[j] = Green0(rbox_nodes_set[curRowNum], cbox_nodes_set[j], wavenum);
    }

    curColNum = get_next_idx_aca(resRowVal, flagCol);
    if (curColNum == -1) return ;
    #pragma omp parallel for
    for (int i=0;i<nrowclu;i++) {
        table[i+curColNum*nrowclu] = resColVal[i] = Green0(rbox_nodes_set[i], cbox_nodes_set[curColNum], wavenum);
    }
    if ( abs( resRowVal[curColNum] ) < 1e-20 )  return ; // too close to zero is regarded as zero matrix
    DComplex temp_max(resRowVal[curColNum]);

    pre_row_idx.push_back(curRowNum);   // record which row is selected in ACA
    pre_col_idx.push_back(curColNum);   // record which col is selected in ACA

    rkmat->V.a = new DComplex [ncolclu];
    rkmat->U.a = new DComplex [nrowclu];

    for (int j=0;j<ncolclu;j++) rkmat->V.a[j] = resRowVal[j]/temp_max;
    for (int i=0;i<nrowclu;i++) rkmat->U.a[i] = resColVal[i];

    frbNormV = LAPACKE_zlange(LAPACK_COL_MAJOR, norm, ncolclu, 1, (MKL_Complex16 *)rkmat->V.a, ncolclu);
    frbNormU = LAPACKE_zlange(LAPACK_COL_MAJOR, norm, nrowclu, 1, (MKL_Complex16 *)rkmat->U.a, nrowclu);
    frbNormZ = frbNormU * frbNormV;
    frbNormZ_2 = frbNormZ * frbNormZ;

    while ( frbNormU*frbNormV > aca_threshold*frbNormZ ) {

        tmpRank++;

        if((curRowNum = get_next_idx_aca(resColVal, flagRow)) == -1) {
            tmpRank--;
            break;
        }

        #pragma omp parallel for
        for (int j=0;j<ncolclu;j++) {
            resRowVal[j] = Green0(rbox_nodes_set[curRowNum], cbox_nodes_set[j], wavenum);
            table[curRowNum+j*nrowclu] = resRowVal[j];
            for (int i=0;i<tmpRank;i++) {
                resRowVal[j] -= rkmat->U.a[curRowNum+i*nrowclu]*rkmat->V.a[j+i*ncolclu];
            }
        }

        if ((curColNum = get_next_idx_aca(resRowVal, flagCol)) == -1) {
            tmpRank--;
            break;
        }

        #pragma omp parallel for
        for (int i=0;i<nrowclu;i++) {
            resColVal[i] = Green0(rbox_nodes_set[i], cbox_nodes_set[curColNum], wavenum);
            table[i+curColNum*nrowclu] = resColVal[i];
            for (int j=0;j<tmpRank;j++) {
                resColVal[i] -= rkmat->U.a[i+j*nrowclu] * rkmat->V.a[curColNum+j*ncolclu];
            }
        }

        if ( abs( resRowVal[curColNum] ) < 1e-20 ) {
            tmpRank--;
            break;
        }
        temp_max = resRowVal[curColNum];

        rkmat->V.a = (DComplex *) realloc(rkmat->V.a, (tmpRank+1)*ncolclu*sizeof(DComplex));
        rkmat->U.a = (DComplex *) realloc(rkmat->U.a, (tmpRank+1)*nrowclu*sizeof(DComplex));

        pre_row_idx.push_back(curRowNum);
        pre_col_idx.push_back(curColNum);

        #pragma omp parallel for
        for (int j=0;j<ncolclu;j++) rkmat->V.a[j+tmpRank*ncolclu] = resRowVal[j]/temp_max;
        #pragma omp parallel for
        for (int i=0;i<nrowclu;i++) rkmat->U.a[i+tmpRank*nrowclu] = resColVal[i];

        frbNormV = LAPACKE_zlange(LAPACK_COL_MAJOR, norm, ncolclu, 1, (MKL_Complex16 *)rkmat->V.a+tmpRank*ncolclu, ncolclu);
        frbNormU = LAPACKE_zlange(LAPACK_COL_MAJOR, norm, nrowclu, 1, (MKL_Complex16 *)rkmat->U.a+tmpRank*nrowclu, nrowclu);
        for (int i=0;i<tmpRank;i++) {
            DComplex c1, c2;
            cblas_zdotc_sub(nrowclu, rkmat->U.a+i*nrowclu, 1, rkmat->U.a+tmpRank*nrowclu, 1, &c1);
            cblas_zdotc_sub(ncolclu, rkmat->V.a+i*ncolclu, 1, rkmat->V.a+tmpRank*ncolclu, 1, &c2);
            frbNormZ_2 += 2*abs( c1 )*abs( c2 );
        }
        frbNormZ_2 += frbNormU*frbNormU*frbNormV*frbNormV;
        frbNormZ = sqrt(frbNormZ_2);
    }

    tmpRank ++;

    cout << "m=" << nrowclu << ", " << "n=" << ncolclu << ": " << tmpRank << endl;

    nrowclu = hm->rc->size;
    ncolclu = hm->cc->size;

    /** actually, the table can be avoided **/
    DComplex *gs_mat = new DComplex [tmpRank*tmpRank];
    #pragma omp parallel for
    for (int i=0;i<tmpRank;++i) {
        for (int j=0;j<tmpRank;++j) {
            gs_mat[i+j*tmpRank] = table[pre_row_idx[i]+pre_col_idx[j]*rbox_nodes_num];
        }
    }/// This S matrix is still low rank, so direct inverse is not feasible.
    delete [] table;


    vector<int> row_idx, col_idx;
    int curRank = tmpRank;
    DComplex *ks_mat = gs_mat;
    for (int i=0;i<tmpRank;++i) {
        row_idx.push_back(pre_row_idx[i]);
        col_idx.push_back(pre_col_idx[i]);
    }
    matInv_z(ks_mat, curRank);


    if ('e' == c) {

        Complex3D *newa_vec = new Complex3D [curRank*nrowclu];
        Complex3D *newb_vec = new Complex3D [curRank*ncolclu];
        prkmatrix rk_vec_x = new_rkmatrix(nrowclu, ncolclu, curRank);
        prkmatrix rk_vec_y = new_rkmatrix(nrowclu, ncolclu, curRank);
        prkmatrix rk_vec_z = new_rkmatrix(nrowclu, ncolclu, curRank);
        prkmatrix rk_sca = new_rkmatrix(nrowclu, ncolclu, curRank);

        #pragma omp parallel for
        for (int j=0;j<curRank;++j) {
            for (int i=0;i<nrowclu;++i) {
                intg.single_layer_intg(hm->rc->idx[i], cbox_nodes_set[col_idx[j]],
                                    newa_vec[i+nrowclu*j], rk_sca->U.a[i+nrowclu*j]);
            }
            for (int i=0;i<ncolclu;++i) {
                intg.single_layer_intg(hm->cc->idx[i], rbox_nodes_set[row_idx[j]],
                                    newb_vec[i+ncolclu*j], rk_sca->V.a[i+ncolclu*j]);
            }
        }

        #pragma omp parallel for
        for (int j=0;j<curRank;++j) {
            for (int i=0;i<nrowclu;++i) {
                rk_vec_x->U.a[i+nrowclu*j] = coefA*newa_vec[i+nrowclu*j][0];
                rk_vec_y->U.a[i+nrowclu*j] = coefA*newa_vec[i+nrowclu*j][1];
                rk_vec_z->U.a[i+nrowclu*j] = coefA*newa_vec[i+nrowclu*j][2];
                rk_sca->U.a[i+nrowclu*j] *= coefV;
            }
            for (int i=0;i<ncolclu;++i) {
                rk_vec_x->V.a[i+ncolclu*j] = conj(newb_vec[i+ncolclu*j][0]);
                rk_vec_y->V.a[i+ncolclu*j] = conj(newb_vec[i+ncolclu*j][1]);
                rk_vec_z->V.a[i+ncolclu*j] = conj(newb_vec[i+ncolclu*j][2]);
                rk_sca->V.a[i+ncolclu*j] = conj(rk_sca->V.a[i+ncolclu*j]);
            }
        }

        prkmatrix rk_vec_nx = new_rkmatrix(nrowclu, ncolclu, curRank);
        prkmatrix rk_vec_ny = new_rkmatrix(nrowclu, ncolclu, curRank);
        prkmatrix rk_vec_nz = new_rkmatrix(nrowclu, ncolclu, curRank);
        prkmatrix rk_sca_n = new_rkmatrix(nrowclu, ncolclu, curRank);

        copy_rkmatrix(0, rk_vec_x, rk_vec_nx);
        copy_rkmatrix(0, rk_vec_y, rk_vec_ny);
        copy_rkmatrix(0, rk_vec_z, rk_vec_nz);
        copy_rkmatrix(0, rk_sca,   rk_sca_n);

        if (nrowclu<=ncolclu) {
            mmp_z(rk_vec_x->U.a, CblasNoTrans, ks_mat, CblasNoTrans, rk_vec_nx->U.a, nrowclu, curRank, curRank);
            mmp_z(rk_vec_y->U.a, CblasNoTrans, ks_mat, CblasNoTrans, rk_vec_ny->U.a, nrowclu, curRank, curRank);
            mmp_z(rk_vec_z->U.a, CblasNoTrans, ks_mat, CblasNoTrans, rk_vec_nz->U.a, nrowclu, curRank, curRank);
            mmp_z(rk_sca->U.a,   CblasNoTrans, ks_mat, CblasNoTrans, rk_sca_n->U.a,  nrowclu, curRank, curRank);
        }
        else {
            mmp_z(rk_vec_x->V.a, CblasNoTrans, ks_mat, CblasConjTrans, rk_vec_nx->V.a, ncolclu, curRank, curRank);
            mmp_z(rk_vec_y->V.a, CblasNoTrans, ks_mat, CblasConjTrans, rk_vec_ny->V.a, ncolclu, curRank, curRank);
            mmp_z(rk_vec_z->V.a, CblasNoTrans, ks_mat, CblasConjTrans, rk_vec_nz->V.a, ncolclu, curRank, curRank);
            mmp_z(rk_sca->V.a,   CblasNoTrans, ks_mat, CblasConjTrans, rk_sca_n->V.a,  ncolclu, curRank, curRank);
        }

        copy_rkmatrix(false, rk_vec_nx, hm->r);

        #if 1
        ptruncmode tm = new_truncmode();
        add_rkmatrix(1.0, rk_vec_ny, tm, svd_threshold, hm->r);
        add_rkmatrix(1.0, rk_vec_nz, tm, svd_threshold, hm->r);
//        cout << hm->r->k << " ";
        add_rkmatrix(1.0, rk_sca_n,  tm, svd_threshold, hm->r);
        #else
        add_rkmatrix_keep_rank(1.0, rk_vec_ny, hm->r);
        add_rkmatrix_keep_rank(1.0, rk_vec_nz, hm->r);
            #if 0
            add_rkmatrix_keep_rank(1.0, rk_sca_n, hm->r);
            #endif // 0
        hm->r->U.a = (DComplex *) realloc(hm->r->U.a, (hm->r->k+rk_sca_n->k)*nrowclu*sizeof(DComplex));
        hm->r->V.a = (DComplex *) realloc(hm->r->V.a, (hm->r->k+rk_sca_n->k)*ncolclu*sizeof(DComplex));
        hm->r->U.cols = hm->r->k+rk_sca_n->k;
        hm->r->V.cols = hm->r->k+rk_sca_n->k;
        #pragma omp parallel for
        for (int i=0;i<rk_sca_n->k;++i) {
            for (int j=0;j<nrowclu;++j) {
                hm->r->U.a[j+(i+hm->r->k)*nrowclu] = rk_sca_n->U.a[j+i*nrowclu];
            }
            for (int j=0;j<ncolclu;++j) {
                hm->r->V.a[j+(i+hm->r->k)*ncolclu] = rk_sca_n->V.a[j+i*ncolclu];
            }
        }
        hm->r->k +=rk_sca_n->k;
        #endif // 0
        del_rkmatrix(rk_vec_nx);
        del_rkmatrix(rk_vec_ny);
        del_rkmatrix(rk_vec_nz);
        del_rkmatrix(rk_sca_n);

        delete [] ks_mat;
        delete [] newa_vec;
        delete [] newb_vec;
        del_rkmatrix(rk_vec_x);
        del_rkmatrix(rk_vec_y);
        del_rkmatrix(rk_vec_z);
        del_rkmatrix(rk_sca);
        del_rkmatrix(rkmat);
    }
    else if ('m' == c) {
        Complex3D *newa_vec = new Complex3D [curRank*nrowclu];
        Complex3D *newb_vec = new Complex3D [curRank*ncolclu];
        prkmatrix rk_vec_x = new_rkmatrix(nrowclu, ncolclu, curRank);
        prkmatrix rk_vec_y = new_rkmatrix(nrowclu, ncolclu, curRank);
        prkmatrix rk_vec_z = new_rkmatrix(nrowclu, ncolclu, curRank);

        #pragma omp parallel for
        for (int j=0;j<curRank;++j) {
            for (int i=0;i<nrowclu;++i) {
                newa_vec[i+nrowclu*j] = intg.single_layer_intg_m1(hm->rc->idx[i], cbox_nodes_set[col_idx[j]]);
            }
            for (int i=0;i<ncolclu;++i) {
                newb_vec[i+ncolclu*j] = intg.single_layer_intg_m2(hm->cc->idx[i], rbox_nodes_set[row_idx[j]]);
            }
        }

        #pragma omp parallel for
        for (int j=0;j<curRank;++j) {
            for (int i=0;i<nrowclu;++i) {
                rk_vec_x->U.a[i+nrowclu*j] = newa_vec[i+nrowclu*j][0];
                rk_vec_y->U.a[i+nrowclu*j] = newa_vec[i+nrowclu*j][1];
                rk_vec_z->U.a[i+nrowclu*j] = newa_vec[i+nrowclu*j][2];
            }
            for (int i=0;i<ncolclu;++i) {
                rk_vec_x->V.a[i+ncolclu*j] = conj(newb_vec[i+ncolclu*j][0]);
                rk_vec_y->V.a[i+ncolclu*j] = conj(newb_vec[i+ncolclu*j][1]);
                rk_vec_z->V.a[i+ncolclu*j] = conj(newb_vec[i+ncolclu*j][2]);
            }
        }

        hm->r->U.a = new DComplex [curRank*nrowclu];
        hm->r->V.a = new DComplex [curRank*ncolclu];

        prkmatrix rk_vec_nx = new_rkmatrix(nrowclu, ncolclu, curRank);
        prkmatrix rk_vec_ny = new_rkmatrix(nrowclu, ncolclu, curRank);
        prkmatrix rk_vec_nz = new_rkmatrix(nrowclu, ncolclu, curRank);
        copy_rkmatrix(false, rk_vec_x, rk_vec_nx);
        copy_rkmatrix(false, rk_vec_y, rk_vec_ny);
        copy_rkmatrix(false, rk_vec_z, rk_vec_nz);

        if (nrowclu<=ncolclu) {
            mmp_z(rk_vec_x->U.a, CblasNoTrans, ks_mat, CblasNoTrans, rk_vec_nx->U.a, nrowclu, curRank, curRank);
            mmp_z(rk_vec_y->U.a, CblasNoTrans, ks_mat, CblasNoTrans, rk_vec_ny->U.a, nrowclu, curRank, curRank);
            mmp_z(rk_vec_z->U.a, CblasNoTrans, ks_mat, CblasNoTrans, rk_vec_nz->U.a, nrowclu, curRank, curRank);
        }
        else {
            mmp_z(rk_vec_x->V.a, CblasNoTrans, ks_mat, CblasConjTrans, rk_vec_nx->V.a, ncolclu, curRank, curRank);
            mmp_z(rk_vec_y->V.a, CblasNoTrans, ks_mat, CblasConjTrans, rk_vec_ny->V.a, ncolclu, curRank, curRank);
            mmp_z(rk_vec_z->V.a, CblasNoTrans, ks_mat, CblasConjTrans, rk_vec_nz->V.a, ncolclu, curRank, curRank);
        }

        copy_rkmatrix(false, rk_vec_nx, hm->r);

        #if 0
        add_rkmatrix_keep_rank(1.0, rk_vec_ny, hm->r);
        add_rkmatrix_keep_rank(1.0, rk_vec_nz, hm->r);
        #else
        ptruncmode tm = new_truncmode();
        add_rkmatrix(1.0, rk_vec_ny, tm, svd_threshold, hm->r);
        add_rkmatrix(1.0, rk_vec_nz, tm, svd_threshold, hm->r);
        #endif // 0

        del_rkmatrix(rk_vec_nx);
        del_rkmatrix(rk_vec_ny);
        del_rkmatrix(rk_vec_nz);

        delete [] ks_mat;
        delete [] newa_vec;
        delete [] newb_vec;
        del_rkmatrix(rk_vec_x);
        del_rkmatrix(rk_vec_y);
        del_rkmatrix(rk_vec_z);
    }
    else if ('c' == c) {
        Complex3D *newa_vec_e = new Complex3D [curRank*nrowclu];
        Complex3D *newb_vec_e = new Complex3D [curRank*ncolclu];
        prkmatrix rk_vec_x_e = new_rkmatrix(nrowclu, ncolclu, curRank);
        prkmatrix rk_vec_y_e = new_rkmatrix(nrowclu, ncolclu, curRank);
        prkmatrix rk_vec_z_e = new_rkmatrix(nrowclu, ncolclu, curRank);
        prkmatrix rk_sca_e = new_rkmatrix(nrowclu, ncolclu, curRank);
        Complex3D *newa_vec_m = new Complex3D [curRank*nrowclu];
        Complex3D *newb_vec_m = new Complex3D [curRank*ncolclu];
        prkmatrix rk_vec_x_m = new_rkmatrix(nrowclu, ncolclu, curRank);
        prkmatrix rk_vec_y_m = new_rkmatrix(nrowclu, ncolclu, curRank);
        prkmatrix rk_vec_z_m = new_rkmatrix(nrowclu, ncolclu, curRank);

        const DComplex eta(ETA0);

        #pragma omp parallel for
        for (int j=0;j<curRank;++j) {
            for (int i=0;i<nrowclu;++i) {
                intg.single_layer_intg(hm->rc->idx[i], cbox_nodes_set[col_idx[j]],
                                    newa_vec_e[i+nrowclu*j], rk_sca_e->U.a[i+nrowclu*j]);
                newa_vec_m[i+nrowclu*j] = eta * intg.single_layer_intg_m1(hm->rc->idx[i], cbox_nodes_set[col_idx[j]]);
            }
            for (int i=0;i<ncolclu;++i) {
                intg.single_layer_intg(hm->cc->idx[i], rbox_nodes_set[row_idx[j]],
                                    newb_vec_e[i+ncolclu*j], rk_sca_e->V.a[i+ncolclu*j]);
                newb_vec_m[i+ncolclu*j] = intg.single_layer_intg_m2(hm->cc->idx[i], rbox_nodes_set[row_idx[j]]);
            }
        }

        #pragma omp parallel for
        for (int j=0;j<curRank;++j) {
            for (int i=0;i<nrowclu;++i) {
                rk_vec_x_e->U.a[i+nrowclu*j] = coefA*newa_vec_e[i+nrowclu*j][0];
                rk_vec_y_e->U.a[i+nrowclu*j] = coefA*newa_vec_e[i+nrowclu*j][1];
                rk_vec_z_e->U.a[i+nrowclu*j] = coefA*newa_vec_e[i+nrowclu*j][2];
                rk_sca_e->U.a[i+nrowclu*j] *= coefV;
                rk_vec_x_m->U.a[i+nrowclu*j] = newa_vec_m[i+nrowclu*j][0];
                rk_vec_y_m->U.a[i+nrowclu*j] = newa_vec_m[i+nrowclu*j][1];
                rk_vec_z_m->U.a[i+nrowclu*j] = newa_vec_m[i+nrowclu*j][2];
            }
            for (int i=0;i<ncolclu;++i) {
                rk_vec_x_e->V.a[i+ncolclu*j] = conj(newb_vec_e[i+ncolclu*j][0]);
                rk_vec_y_e->V.a[i+ncolclu*j] = conj(newb_vec_e[i+ncolclu*j][1]);
                rk_vec_z_e->V.a[i+ncolclu*j] = conj(newb_vec_e[i+ncolclu*j][2]);
                rk_sca_e->V.a[i+ncolclu*j] = conj(rk_sca_e->V.a[i+ncolclu*j]);
                rk_vec_x_m->V.a[i+ncolclu*j] = conj(newb_vec_m[i+ncolclu*j][0]);
                rk_vec_y_m->V.a[i+ncolclu*j] = conj(newb_vec_m[i+ncolclu*j][1]);
                rk_vec_z_m->V.a[i+ncolclu*j] = conj(newb_vec_m[i+ncolclu*j][2]);
            }
        }

        prkmatrix rk_vec_nx_e = new_rkmatrix(nrowclu, ncolclu, curRank);
        prkmatrix rk_vec_ny_e = new_rkmatrix(nrowclu, ncolclu, curRank);
        prkmatrix rk_vec_nz_e = new_rkmatrix(nrowclu, ncolclu, curRank);
        prkmatrix rk_sca_n_e = new_rkmatrix(nrowclu, ncolclu, curRank);
        copy_rkmatrix(false, rk_vec_x_e, rk_vec_nx_e);
        copy_rkmatrix(false, rk_vec_y_e, rk_vec_ny_e);
        copy_rkmatrix(false, rk_vec_z_e, rk_vec_nz_e);
        copy_rkmatrix(false, rk_sca_e,   rk_sca_n_e);
        prkmatrix rk_vec_nx_m = new_rkmatrix(nrowclu, ncolclu, curRank);
        prkmatrix rk_vec_ny_m = new_rkmatrix(nrowclu, ncolclu, curRank);
        prkmatrix rk_vec_nz_m = new_rkmatrix(nrowclu, ncolclu, curRank);
        copy_rkmatrix(false, rk_vec_x_m, rk_vec_nx_m);
        copy_rkmatrix(false, rk_vec_y_m, rk_vec_ny_m);
        copy_rkmatrix(false, rk_vec_z_m, rk_vec_nz_m);

        if (nrowclu<=ncolclu) {
            mmp_z(rk_vec_x_e->U.a, CblasNoTrans, ks_mat, CblasNoTrans, rk_vec_nx_e->U.a, nrowclu, curRank, curRank);
            mmp_z(rk_vec_y_e->U.a, CblasNoTrans, ks_mat, CblasNoTrans, rk_vec_ny_e->U.a, nrowclu, curRank, curRank);
            mmp_z(rk_vec_z_e->U.a, CblasNoTrans, ks_mat, CblasNoTrans, rk_vec_nz_e->U.a, nrowclu, curRank, curRank);
            mmp_z(rk_sca_e->U.a,   CblasNoTrans, ks_mat, CblasNoTrans, rk_sca_n_e->U.a,  nrowclu, curRank, curRank);
            mmp_z(rk_vec_x_m->U.a, CblasNoTrans, ks_mat, CblasNoTrans, rk_vec_nx_m->U.a, nrowclu, curRank, curRank);
            mmp_z(rk_vec_y_m->U.a, CblasNoTrans, ks_mat, CblasNoTrans, rk_vec_ny_m->U.a, nrowclu, curRank, curRank);
            mmp_z(rk_vec_z_m->U.a, CblasNoTrans, ks_mat, CblasNoTrans, rk_vec_nz_m->U.a, nrowclu, curRank, curRank);
        }
        else {
            mmp_z(rk_vec_x_e->V.a, CblasNoTrans, ks_mat, CblasConjTrans, rk_vec_nx_e->V.a, ncolclu, curRank, curRank);
            mmp_z(rk_vec_y_e->V.a, CblasNoTrans, ks_mat, CblasConjTrans, rk_vec_ny_e->V.a, ncolclu, curRank, curRank);
            mmp_z(rk_vec_z_e->V.a, CblasNoTrans, ks_mat, CblasConjTrans, rk_vec_nz_e->V.a, ncolclu, curRank, curRank);
            mmp_z(rk_sca_e->V.a,   CblasNoTrans, ks_mat, CblasConjTrans, rk_sca_n_e->V.a,  ncolclu, curRank, curRank);
            mmp_z(rk_vec_x_m->V.a, CblasNoTrans, ks_mat, CblasConjTrans, rk_vec_nx_m->V.a, ncolclu, curRank, curRank);
            mmp_z(rk_vec_y_m->V.a, CblasNoTrans, ks_mat, CblasConjTrans, rk_vec_ny_m->V.a, ncolclu, curRank, curRank);
            mmp_z(rk_vec_z_m->V.a, CblasNoTrans, ks_mat, CblasConjTrans, rk_vec_nz_m->V.a, ncolclu, curRank, curRank);
        }

        ptruncmode tm = new_truncmode();
        add_rkmatrix(1.0, rk_vec_ny_e, tm, svd_threshold, rk_vec_nx_e);
        add_rkmatrix(1.0, rk_vec_nz_e, tm, svd_threshold, rk_vec_nx_e);
        add_rkmatrix(1.0, rk_sca_n_e,  tm, svd_threshold, rk_vec_nx_e);
        add_rkmatrix(1.0, rk_vec_ny_m,  tm, svd_threshold, rk_vec_nx_m);
        add_rkmatrix(1.0, rk_vec_nz_m,  tm, svd_threshold, rk_vec_nx_m);

        copy_rkmatrix(false, rk_vec_nx_e, hm->r);
        add_rkmatrix(1.0, rk_vec_nx_m, tm, svd_threshold, hm->r);

        del_rkmatrix(rk_vec_nx_e);
        del_rkmatrix(rk_vec_ny_e);
        del_rkmatrix(rk_vec_nz_e);
        del_rkmatrix(rk_sca_n_e);
        del_rkmatrix(rk_vec_nx_m);
        del_rkmatrix(rk_vec_ny_m);
        del_rkmatrix(rk_vec_nz_m);
        delete [] ks_mat;
        delete [] newa_vec_e;
        delete [] newb_vec_e;
        delete [] newa_vec_m;
        delete [] newb_vec_m;
        del_rkmatrix(rk_vec_x_e);
        del_rkmatrix(rk_vec_y_e);
        del_rkmatrix(rk_vec_z_e);
        del_rkmatrix(rk_sca_e);
        del_rkmatrix(rk_vec_x_m);
        del_rkmatrix(rk_vec_y_m);
        del_rkmatrix(rk_vec_z_m);
        del_rkmatrix(rkmat);
    }
    else if ('p' == c) {
    }
    else {
        cerr << "Unrecognized parameter for HCA:" << c << " ?" << endl;
        exit(1);
    }

    cout << "m = " << nrowclu << ", n = " << ncolclu << ", rank = " << hm->r->k << endl;

    #if TST_HMATRIX_ERROR
    /** Full SVD to verify the performance of ACA **/
    /** Or just to compute the error matrix ||Z-Z^hat||2 directly  **/
    DComplex *Z = new DComplex [nrowclu*ncolclu];
    DComplex *aca_Z = new DComplex [nrowclu*ncolclu];

    mmp_z(hm->r->U.a, CblasNoTrans, hm->r->V.a, CblasConjTrans, aca_Z, nrowclu, hm->r->k, ncolclu);
    #pragma omp parallel for
    for (int i=0;i<nrowclu;++i) {
        for (int j=0;j<ncolclu;++j) {
            Z[i+j*nrowclu] = intg.xfie(hm->rc->idx[i], hm->cc->idx[j], c);
        }
    }
    for (int i=0;i<nrowclu;++i) {
        for (int j=0;j<ncolclu;++j) {
            DComplex tmp = Z[i+j*nrowclu];
            frob_norm_2_total += abs(tmp*tmp);
            tmp -= aca_Z[i+j*nrowclu];
            frob_norm_error_2_abs += abs(tmp*tmp);
        }
    }

    #if 0
    double *sigma = new double [max(nrowclu, ncolclu)];
    double *work = new double [max(nrowclu, ncolclu)];
    MKL_Complex16 *u, *vt;
    LAPACKE_zgesvd(h2_col_major, 'N', 'N', nrowclu, ncolclu, (MKL_Complex16 *)Z, nrowclu, sigma, u, 1, vt, 1, work);
    int real_rank = min(nrowclu, ncolclu);
    for (int i=0;i<real_rank;++i)
        cout << i<< ": " << sigma[i] << endl;
    delete [] sigma;
    delete [] work;
    #endif // 0

    delete [] Z;
    delete [] aca_Z;

    #endif // 0
}

void
FastBem::aca(phmatrix hm, const Integration &intg, const char c)
{
//    cout << aca_times++ << ": use ACA" << endl;

    /** ACA required parameters **/
    char norm= 'F';
    const int nrowclu(hm->rc->size);
    const int ncolclu(hm->cc->size);

    /** auxiliary variables **/
    int curRank(0);
    double frbNormZ, frbNormZ_2, frbNormU, frbNormV;
    vector<bool> flagRow(nrowclu, 1);
    vector<bool> flagCol(ncolclu, 1);
    vector<unsigned> row_indices_set;
    vector<unsigned> col_indices_set;

    /** record the info of the current row and col **/
    int curRowNum, curColNum;
    vector<DComplex> resRowVal(ncolclu, 0);
    vector<DComplex> resColVal(nrowclu, 0);

    /** Let's BEGIN!
        Start from the main rk-matrix **/
    const prkmatrix rkmat = hm->r;

    /** Start from the 0th row **/
    /** Revise to random version ( incomplete ) **/
    curRowNum = 0;
    flagRow[curRowNum] = 0;
    #pragma omp parallel for
    for (int j=0;j<ncolclu;j++) resRowVal[j] = intg.xfie(hm->rc->idx[curRowNum], hm->cc->idx[j], c);

    curColNum = get_next_idx_aca(resRowVal, flagCol);
    if (curColNum == -1) return ;   /** Return value -1 means game over ;) **/

    #pragma omp parallel for
    for (int i=0;i<nrowclu;i++) resColVal[i] = intg.xfie(hm->rc->idx[i], hm->cc->idx[curColNum], c);
    DComplex temp_max = resRowVal[curColNum];
    if ( abs( temp_max ) < 1e-20 )  return ; /** too close to zero, relative value is perferable **/

    row_indices_set.push_back(hm->rc->idx[curRowNum]);
    col_indices_set.push_back(hm->cc->idx[curColNum]);

    rkmat->V.a = new DComplex [ncolclu]; /** now V is still Vt **/
    rkmat->U.a = new DComplex [nrowclu]; /** rank 1 mat now is allocated **/

    #pragma omp parallel for
    for (int j=0;j<ncolclu;j++) rkmat->V.a[j] = resRowVal[j]/temp_max;
    #pragma omp parallel for
    for (int i=0;i<nrowclu;i++) rkmat->U.a[i] = resColVal[i];

    frbNormV = LAPACKE_zlange(LAPACK_COL_MAJOR, norm, ncolclu, 1, (MKL_Complex16 *)rkmat->V.a, ncolclu);
    frbNormU = LAPACKE_zlange(LAPACK_COL_MAJOR, norm, nrowclu, 1, (MKL_Complex16 *)rkmat->U.a, nrowclu);
    frbNormZ = frbNormU * frbNormV;
    frbNormZ_2 = frbNormZ * frbNormZ;
    /// up to now, the U1, V1, and residual matrix have been set up
    /// frobenius norms have been initilized by the rank 1 matrix

    while ( frbNormU*frbNormV > aca_threshold*frbNormZ ) { /// Stop criterion may have the potential to be improved.

        curRank++;

        if((curRowNum = get_next_idx_aca(resColVal, flagRow)) == -1) {
            curRank--;
            break;
        }

        #pragma omp parallel for
        for (int j=0;j<ncolclu;j++) {
            /// Here the original matrix values are discarded.
            resRowVal[j] = intg.xfie(hm->rc->idx[curRowNum], hm->cc->idx[j], c);
            for (int i=0;i<curRank;i++) {
                resRowVal[j] -= rkmat->U.a[curRowNum+i*nrowclu]*rkmat->V.a[j+i*ncolclu];
            }
            /// the Row Vector of the residual matrix is updated by the rank-k matrix
        }

        if ((curColNum = get_next_idx_aca(resRowVal, flagCol)) == -1) {
            curRank--;
            break;
        }

        #pragma omp parallel for
        for (int i=0;i<nrowclu;i++) {
            resColVal[i] = intg.xfie(hm->rc->idx[i], hm->cc->idx[curColNum], c);
            for (int j=0;j<curRank;j++) {
                resColVal[i] -= rkmat->U.a[i+j*nrowclu] * rkmat->V.a[curColNum+j*ncolclu];
            }
            /// the Col Vector of the residual matrix is updated by the rank-k matrix
        }

        temp_max = resRowVal[curColNum];

        if (abs(temp_max) < 1e-20) {
            curRank--;
            break;
        }

        row_indices_set.push_back(hm->rc->idx[curRowNum]);
        col_indices_set.push_back(hm->cc->idx[curColNum]);

        // allocated one more row and column space for rank-k mat to rank-(k+1) mat
        rkmat->V.a = (DComplex *) realloc(rkmat->V.a, (curRank+1)*ncolclu*sizeof(DComplex));
        rkmat->U.a = (DComplex *) realloc(rkmat->U.a, (curRank+1)*nrowclu*sizeof(DComplex));

        #pragma omp parallel for
        for (int j=0;j<ncolclu;j++) rkmat->V.a[j+curRank*ncolclu] = resRowVal[j]/temp_max;
        #pragma omp parallel for
        for (int i=0;i<nrowclu;i++) rkmat->U.a[i+curRank*nrowclu] = resColVal[i];

        // update the newly formed rank-k matrix, without implicitly computation
        frbNormV = LAPACKE_zlange(LAPACK_COL_MAJOR, norm, ncolclu, 1, (MKL_Complex16 *)rkmat->V.a+curRank*ncolclu, ncolclu);
        frbNormU = LAPACKE_zlange(LAPACK_COL_MAJOR, norm, nrowclu, 1, (MKL_Complex16 *)rkmat->U.a+curRank*nrowclu, nrowclu);

//        #pragma omp parallel for
//        add with something in parallel for will lead to errors
        for (int i=0;i<curRank;i++) {
            DComplex c1, c2;
            cblas_zdotc_sub(nrowclu, rkmat->U.a+i*nrowclu, 1, rkmat->U.a+curRank*nrowclu, 1, &c1);
            cblas_zdotc_sub(ncolclu, rkmat->V.a+i*ncolclu, 1, rkmat->V.a+curRank*ncolclu, 1, &c2);
            frbNormZ_2 += 2*abs( c1 )*abs( c2 );
        }
        frbNormZ_2 += frbNormU*frbNormU*frbNormV*frbNormV;
        frbNormZ = sqrt(frbNormZ_2);
    }

    /// set up parameters for the final rank-k matrix
    rkmat->k = curRank+1;
    rkmat->U.rows = nrowclu;
    rkmat->U.cols = rkmat->k;
    rkmat->V.rows = ncolclu;
	rkmat->V.cols = rkmat->k;

	conj_amatrix(&(rkmat->V));  // convert U*V^T to U*V^H
//	cout << "m=" << nrowclu << ", n=" << ncolclu << ", rank=" << rkmat->k << endl;


    #if TST_HMATRIX_ERROR
    /** Full SVD to verify the performance of ACA **/
    /** Or just to compute the error matrix ||Z-Z^hat||2 directly  **/
    DComplex *Z = new DComplex [nrowclu*ncolclu];
    DComplex *aca_Z = new DComplex [nrowclu*ncolclu];

    mmp_z(rkmat->U.a, CblasNoTrans, rkmat->V.a, CblasConjTrans, aca_Z, rkmat->U.rows, rkmat->k, rkmat->V.rows);
    #pragma omp parallel for
    for (int i=0;i<nrowclu;++i) {
        for (int j=0;j<ncolclu;++j) {
            Z[i+j*nrowclu] = intg.xfie(hm->rc->idx[i], hm->cc->idx[j], c);
        }
    }
    for (int i=0;i<nrowclu;++i) {
        for (int j=0;j<ncolclu;++j) {
            DComplex tmp = Z[i+j*nrowclu];
            frob_norm_2_total += abs(tmp*tmp);
            tmp -= aca_Z[i+j*nrowclu];
            frob_norm_error_2_abs += abs(tmp*tmp);
        }
    }

    #if 0
    double *sigma = new double [max(nrowclu, ncolclu)];
    double *work = new double [max(nrowclu, ncolclu)];
    MKL_Complex16 *u, *vt;
    LAPACKE_zgesvd(h2_col_major, 'N', 'N', nrowclu, ncolclu, (MKL_Complex16 *)Z, nrowclu, sigma, u, 1, vt, 1, work);
    int real_rank = min(nrowclu, ncolclu);
    for (int i=0;i<real_rank;++i)
        cout << i<< ": " << sigma[i] << endl;
    delete [] sigma;
    delete [] work;
    #endif // 0

    delete [] Z;
    delete [] aca_Z;

    #endif // 0
}

void
FastBem::aca_plus(phmatrix hm, const Integration &intg, const char c)
{
    cout << aca_plus_times++ << ": use ACA+" << endl;

    /** ACA needed parameters **/
    char norm = 'F';
    double aca_plus_threshold = aca_threshold;

    const int m(hm->rc->size);
    const int n(hm->cc->size);

    int k(0);
    double frbNormZ, frbNormZ_2, frbNormU, frbNormV;

    /** record the current row and col index, reference row and col index **/
    int curRowNum, curColNum;
    int ref_col_idx, ref_row_idx;

    /** examine whetehr reference row and col index conflict with current selected row and col index **/
    bool use_ref_col, use_ref_row;

    /**  allocate temp space for ACA+  **/
    vector<bool> flagRow(m, true);
    vector<bool> flagCol(n, true);

    vector<int> row_indices_set;
    vector<int> col_indices_set;

    vector<DComplex> resRowVal(n, 0); // residual matrix row vector
    vector<DComplex> resColVal(m, 0); // residual matrix col vector

    vector<DComplex> refRowVal(n, 0); // reference row vector
    vector<DComplex> refColVal(m, 0); // reference col vector

    /** the main rank-k matrix **/
    const prkmatrix rkmat = hm->r;

    #if 0
    DComplex *DirectZ = new DComplex [m*n];
    for (int i=0;i<m;++i) {
        for (int j=0;j<n;++j) {
            DirectZ[i+j*m] = intg.xfie(hm->rc->idx[i], hm->cc->idx[j], c);
        }
    }
    #endif // 1

    /// initilize the residua matrix and rank-1 matrix (U1, V1)
    ref_col_idx = get_ref_idx_rand(flagCol, k); // The reference column is selected randomly
    for (int i=0;i<m;++i) refColVal[i] = intg.xfie(hm->rc->idx[i], hm->cc->idx[ref_col_idx], c);

    ref_row_idx = get_ref_idx_min(refColVal, flagRow);  // the reference row is selected by the minimum value of the referernce column
    for (int j=0;j<n;j++) refRowVal[j] = intg.xfie(hm->rc->idx[ref_row_idx], hm->cc->idx[j], c);

    curRowNum = get_next_idx_aca_plus(refColVal, flagRow);
    curColNum = get_next_idx_aca_plus(refRowVal, flagCol);

    if (curColNum >= 0) {
        if (curRowNum>=0) {
            if ( abs(refColVal[curRowNum]) > abs(refRowVal[curColNum]) ) {
                use_ref_col = true;
                use_ref_row = false;
            }
            else {
                use_ref_col = false;
                use_ref_row = true;
            }
        }
        else {
            use_ref_col = false;
            use_ref_row = true;
        }
    }
    else {
        if (curRowNum>=0) {
            use_ref_col = true;
            use_ref_row = false;
        }
        else {
            use_ref_col = use_ref_row = false;
            cout << "Full Zero Matrix!" << endl;
            #if 0
            double frobNorm = LAPACKE_zlange(LAPACK_COL_MAJOR, norm, m, n, (MKL_Complex16 *)DirectZ, m);
            total_diff += frobNorm * frobNorm;
            #endif // 0
            return ;
        }
    }

    if ( use_ref_col ) {
        for (int j=0;j<n;++j) resRowVal[j] = intg.xfie(hm->rc->idx[curRowNum], hm->cc->idx[j], c);
        flagRow[curRowNum] = 0;

        curColNum = get_next_idx_aca(resRowVal, flagCol);
        if ( curColNum==-1 ) {
            cout << "Ref Col: Full Zeroes!" << endl;
            assert(0);
        }
        else {
            for (int i=0;i<m;++i) resColVal[i] = intg.xfie(hm->rc->idx[i], hm->cc->idx[curColNum], c);
        }
    }
    if ( use_ref_row ) {
        for (int i=0;i<m;++i) resColVal[i] = intg.xfie(hm->rc->idx[i], hm->cc->idx[curColNum], c);
        flagCol[curColNum] = 0;

        curRowNum = get_next_idx_aca(resColVal, flagRow);
        if (curRowNum == -1) {
            cout << "Ref Row: Full Zeroes!" << endl;
            assert(0);
        }
        else {
            for (int j=0;j<n;++j) resRowVal[j] = intg.xfie(hm->rc->idx[curRowNum], hm->cc->idx[j], c);
        }
    }

    if ( abs( resRowVal[curColNum] ) < 1e-20 ) {
        #if 0
        double frobNorm = LAPACKE_zlange(LAPACK_COL_MAJOR, norm, m, n, (MKL_Complex16 *)DirectZ, m);
        total_diff += frobNorm * frobNorm;
        #endif // 0
        return ;
    }

    row_indices_set.push_back(hm->rc->idx[curRowNum]);
    col_indices_set.push_back(hm->cc->idx[curColNum]);

    rkmat->V.a = new DComplex [n];
    rkmat->U.a = new DComplex [m];

    for (int j=0;j<n;j++) rkmat->V.a[j] = resRowVal[j]/resRowVal[curColNum];
    for (int i=0;i<m;i++) rkmat->U.a[i] = resColVal[i];

    frbNormV = LAPACKE_zlange(LAPACK_COL_MAJOR, norm, n, 1, (MKL_Complex16 *)rkmat->V.a, n);
    frbNormU = LAPACKE_zlange(LAPACK_COL_MAJOR, norm, m, 1, (MKL_Complex16 *)rkmat->U.a, m);
    frbNormZ = frbNormU*frbNormV;
    frbNormZ_2 = frbNormZ*frbNormZ;

/// init R U1 V1, done!

    while ( frbNormU*frbNormV > aca_plus_threshold*frbNormZ ) {
        k++;

        // solve the conflict problem
        if (curColNum==ref_col_idx) {
            if (curRowNum==ref_row_idx) {
                ref_col_idx = get_ref_idx_rand(flagCol, k);
                if ( ref_col_idx!=-1 ) {
                    #pragma omp parallel for
                    for (int i=0;i<m;++i) {
                        refColVal[i] = intg.xfie(hm->rc->idx[i], hm->cc->idx[ref_col_idx], c);
                        for (int j=0;j<k;j++) refColVal[i] -= rkmat->U.a[i+j*m] * rkmat->V.a[ref_col_idx+j*n];
                    }
                }
                else {
                    k--;
                    break;
                }
                ref_row_idx = get_ref_idx_min(refColVal, flagRow);
                if ( ref_row_idx!=-1 ) {
                    #pragma omp parallel for
                    for (int j=0;j<n;++j) {
                        refRowVal[j] = intg.xfie(hm->rc->idx[ref_row_idx], hm->cc->idx[j], c);
                        for (int i=0;i<k;i++) refRowVal[j] -= rkmat->U.a[ref_row_idx+i*m]*rkmat->V.a[j+i*n];
                    }
                }
                else {
                    k--;
                    break;
                }
            }
            else {
                ref_col_idx = get_ref_idx_min(refRowVal, flagCol);
                if ( ref_col_idx!=-1 ) {
                    #pragma omp parallel for
                    for (int i=0;i<m;++i) {
                        refColVal[i] = intg.xfie(hm->rc->idx[i], hm->cc->idx[ref_col_idx], c);
                        for (int j=0;j<k;j++) refColVal[i] -= rkmat->U.a[i+j*m] * rkmat->V.a[ref_col_idx+j*n];
                    }
                }
                else {
                    k--;
                    break;
                }
                #pragma omp parallel for
                for (int j=0;j<n;++j) {
                    for (int i=0;i<k;i++) refRowVal[j] -= rkmat->U.a[ref_row_idx+i*m]*rkmat->V.a[j+i*n];
                }
            }
        }
        else {
            if (curRowNum==ref_row_idx) {
                #pragma omp parallel for
                for (int i=0;i<m;++i) {
                    for (int j=0;j<k;j++) refColVal[i] -= rkmat->U.a[i+j*m]*rkmat->V.a[ref_col_idx+j*n];
                }
                ref_row_idx = get_ref_idx_min(refColVal, flagRow);
                if ( ref_row_idx!=-1 ) {
                    #pragma omp parallel for
                    for (int j=0;j<n;++j) {
                        refRowVal[j] = intg.xfie(hm->rc->idx[ref_row_idx], hm->cc->idx[j], c);
                        for (int i=0;i<k;i++) refRowVal[j] -= rkmat->U.a[ref_row_idx+i*m]*rkmat->V.a[j+i*n];
                    }
                }
                else {
                    k--;
                    break;
                }
            }
            else {
                #pragma omp parallel for
                for (int i=0;i<m;++i) {
                    for (int j=0;j<k;j++) refColVal[i] -= rkmat->U.a[i+j*m]*rkmat->V.a[ref_col_idx+j*n];
                }
                #pragma omp parallel for
                for (int j=0;j<n;++j) {
                    for (int i=0;i<k;i++) refRowVal[j] -= rkmat->U.a[ref_row_idx+i*m]*rkmat->V.a[j+i*n];
                }
            }
        }

        curRowNum = get_next_idx_aca_plus(refColVal, flagRow);  // select current row index
        curColNum = get_next_idx_aca_plus(refRowVal, flagCol);  // select current col index

        /** judge which pivot element to use **/
        if (curColNum >= 0) {
            if (curRowNum >= 0) {
                if ( abs(refColVal[curRowNum]) > abs(refRowVal[curColNum]) ) {
                    use_ref_col = true;
                    use_ref_row = false;
                }
                else {
                    use_ref_col = false;
                    use_ref_row = true;
                }
            }
            else {
                use_ref_col = false;
                use_ref_row = true;
            }
        }
        else {
            if (curRowNum>=0) {
                use_ref_col = true;
                use_ref_row = false;
            }
            else {
                k--;
                cout << "Full Zero Matrix!" << endl;
                break;
            }
        }

        if ( use_ref_col ) {
            #pragma omp parallel for
            for (int j=0;j<n;++j) {
                resRowVal[j] = intg.xfie(hm->rc->idx[curRowNum], hm->cc->idx[j], c);
                for (int i=0;i<k;i++) resRowVal[j] -= rkmat->U.a[curRowNum+i*m]*rkmat->V.a[j+i*n];
            }
            flagRow[curRowNum] = 0;

            curColNum = get_next_idx_aca(resRowVal, flagCol);
            if (curColNum == -1) {
                k--;
                break;
            }
            else {
                #pragma omp parallel for
                for (int i=0;i<m;++i) {
                    resColVal[i] = intg.xfie(hm->rc->idx[i], hm->cc->idx[curColNum], c);
                    for (int j=0;j<k;j++) resColVal[i] -= rkmat->U.a[i+j*m] * rkmat->V.a[curColNum+j*n];
                }
            }
        }
        if ( use_ref_row ) {
            #pragma omp parallel for
            for (int i=0;i<m;++i) {
                resColVal[i] = intg.xfie(hm->rc->idx[i], hm->cc->idx[curColNum], c);
                for (int j=0;j<k;j++) resColVal[i] -= rkmat->U.a[i+j*m] * rkmat->V.a[curColNum+j*n];
            }
            flagCol[curColNum] = 0;

            curRowNum = get_next_idx_aca(resColVal, flagRow);
            if ( curRowNum==-1 ) {
                k--;
                break;
            }
            else {
                #pragma omp parallel for
                for (int j=0;j<n;++j) {
                    resRowVal[j] = intg.xfie(hm->rc->idx[curRowNum], hm->cc->idx[j], c);
                    for (int i=0;i<k;i++) resRowVal[j] -= rkmat->U.a[curRowNum+i*m]*rkmat->V.a[j+i*n];
                }
            }
        }

        if (abs(resRowVal[curColNum]) < 1e-20) {
            k--;
            break;
        }

    row_indices_set.push_back(hm->rc->idx[curRowNum]);
    col_indices_set.push_back(hm->cc->idx[curColNum]);

//        assert(abs(resRowVal[curColNum]) > 1e-20);

        rkmat->V.a = (DComplex *) realloc(rkmat->V.a, (k+1)*n*sizeof(DComplex));
        rkmat->U.a = (DComplex *) realloc(rkmat->U.a, (k+1)*m*sizeof(DComplex));

        /** update the residual matrix **/
        #pragma omp parallel for
        for (int j=0;j<n;j++) rkmat->V.a[j+k*n] = resRowVal[j]/resRowVal[curColNum];
        #pragma omp parallel for
        for (int i=0;i<m;i++) rkmat->U.a[i+k*m] = resColVal[i];

        /** update the frobenius norm **/
        frbNormV = LAPACKE_zlange(LAPACK_COL_MAJOR, norm, n, 1, (MKL_Complex16 *)rkmat->V.a+k*n, n);
        frbNormU = LAPACKE_zlange(LAPACK_COL_MAJOR, norm, m, 1, (MKL_Complex16 *)rkmat->U.a+k*m, m);
        #pragma omp parallel for
        for (int i=0;i<k;i++) {
            DComplex c1, c2;
            cblas_zdotc_sub(m, rkmat->U.a+i*m, 1, rkmat->U.a+k*m, 1, &c1);
            cblas_zdotc_sub(n, rkmat->V.a+i*n, 1, rkmat->V.a+k*n, 1, &c2);
            frbNormZ_2 += 2 * abs( c1 ) * abs( c2 );
        }
        frbNormZ_2 += frbNormU*frbNormU * frbNormV*frbNormV;
        frbNormZ = sqrt(frbNormZ_2);
    }

    rkmat->k = k+1;
    rkmat->U.rows = m;
    rkmat->U.cols = rkmat->k;
    rkmat->V.rows = n;
	rkmat->V.cols = rkmat->k;

	conj_amatrix(&(rkmat->V));
	cout << "m=" << m << ", n=" << n << ", rank=" << rkmat->k << endl;

#if TST_FROB_NORM
    DComplex *ACAZ = new DComplex [m*n];
    mmp_z(rkmat->U.a, CblasNoTrans, rkmat->V.a, CblasConjTrans, ACAZ, m, k, n);

    for (int j=0;j<n;++j) {
        for (int i=0;i<m;++i) {
            ACAZ[i+j*m] = DirectZ[i+j*m] - ACAZ[i+j*m];
        }
    }

    double frobNorm = LAPACKE_zlange(LAPACK_COL_MAJOR, norm, m, n, (MKL_Complex16 *)ACAZ, m);
    total_diff += frobNorm*frobNorm;
    if (c == 'e')
        cout << sqrt(total_diff) << endl << endl;
    else
        cout << sqrt(total_diff)*ETA0 << endl << endl;

        /// 那些早早就 return 的，岂不是差距更大。
    delete [] DirectZ;
    delete [] ACAZ;
#endif // TST_FROB_NORM
    assert(row_indices_set.size() == col_indices_set.size());
    sort(row_indices_set.begin(), row_indices_set.end());
    for (unsigned i=0;i<row_indices_set.size();++i) {
        cout << row_indices_set[i] << " ";
    }
    cout << endl;
    sort(col_indices_set.begin(), col_indices_set.end());
    for (unsigned i=0;i<col_indices_set.size();++i) {
        cout << col_indices_set[i] << " ";
    }
    cout << endl;
}

void
FastBem::cost_reduction(const double eps)
{
    max_rank_in_rk = 0;
    cost_reduction_hmat(g_hmatrix, eps);
}

void
FastBem::cost_reduction_hmat(phmatrix hm, const double eps)
{
    if (hm->son) {
        for (int j=0;j<hm->csons;++j) {
            for (int i=0;i<hm->rsons;++i) {
                cost_reduction_hmat(hm->son[i+j*hm->rsons], eps);
            }
        }
    }
    else if (hm->r) {
        ptruncmode tm = new_truncmode();
        trunc_rkmatrix(tm, eps, hm->r);
        if (max_rank_in_rk < hm->r->k) max_rank_in_rk = hm->r->k;
    }
}

int
full_col_pivot_aca(const DComplex *A, vector<int> &row_flag, vector<int> &col_flag)
{
    int k(0);
    char norm_type('F');
    int nrow = row_flag.size();
    int ncol = col_flag.size();

    DComplex tempA[nrow*ncol];

    for (int i =0;i<nrow;++i)
        for (int j=0;j<ncol;++j)
            tempA[i+nrow*j] = A[i+nrow*j];

    DComplex temp_row_val[ncol];
    DComplex temp_col_val[nrow];

    double norm = LAPACKE_zlange(LAPACK_COL_MAJOR, norm_type, nrow, ncol, (MKL_Complex16 *)tempA, nrow);

    while (norm > 1e-10) {

        DComplex temp_max(0);
        int i_max(-1), j_max(-1);
        for (int i=0;i<nrow;++i) {
            for (int j=0;j<ncol;++j) {
                if ( abs(tempA[i+nrow*j])>abs(temp_max) ) {
                    temp_max = tempA[i+nrow*j];
                    i_max = i;
                    j_max = j;
                }
            }
        }
        if (i_max == -1 || j_max == -1) break;
        else {
            k++;
            row_flag[i_max] = 0;
            col_flag[j_max] = 0;
            for (int i=0;i<nrow;++i) temp_col_val[i] = tempA[i+nrow*j_max];
            for (int j=0;j<ncol;++j) temp_row_val[j] = tempA[i_max+nrow*j] / temp_max;

            for (int i=0;i<nrow;++i) {
                for (int j=0;j<ncol;++j) {
                    tempA[i+nrow*j] -= temp_col_val[i] * temp_row_val[j];
                }
            }
            norm = LAPACKE_zlange(LAPACK_COL_MAJOR, norm_type, nrow, ncol, (MKL_Complex16 *)tempA, nrow);
        }
    }

    return k;
}


int
get_next_idx_aca(const vector<DComplex> &vals, vector<bool> &flag)
{
    assert(vals.size()==flag.size());
    unsigned len( flag.size() );  // store length of the array
    unsigned max_val_idx(0);        // store the index of the maximum value
    for (max_val_idx=0;max_val_idx<len;++max_val_idx) {
        if ( true == flag[max_val_idx] ) {
            break;          // find the first available index
        }
    }
    if (max_val_idx == len) return -1;  // all the elements have been labeled
    for (unsigned i=max_val_idx+1;i<len;++i) {
        if (true == flag[i]) {
            if (abs(vals[i]) > abs(vals[max_val_idx]))
                max_val_idx = i;        // find the index of the maximum value
        }
    }
    if (max_val_idx == len) return -1;  // failed to find the maximum value, impossible to excute
    else {
        flag[max_val_idx] = false;      // label the current maximum value in the flag array
        return max_val_idx;             // and return its index
    }
}

int
get_next_idx_aca_plus(const vector<DComplex> &vals, const vector<bool> &flag)
{
    assert(vals.size()==flag.size());
    int len = flag.size();
    int max_val_idx;
    for (max_val_idx=0;max_val_idx<len;++max_val_idx) {
        if ( flag[max_val_idx]==1 ) {
            break;
        }
    }
    if (max_val_idx == len) return -1;
    for (int i=max_val_idx+1;i<len;++i) {
        if (flag[i]) {
            if (abs(vals[i]) > abs(vals[max_val_idx]))
                max_val_idx = i;
        }
    }
    if (max_val_idx == len) return -1;
    else {
        /** flag[max_val_idx] = 0; **/
        // the only difference with the get_next_idx_aca is the set of flag array
        // because this returned index may not be selected in the end
        return max_val_idx;
    }
}

int
get_ref_idx_rand(const vector<bool> &flag, const int k)
{
    int final_idx(0);
    int len(flag.size());
#if 1
    srand((unsigned)time(0));

    if (len == k) return -1;

    int idx_count = rand() % (len-k);

    for (int i=0;i<idx_count;++i) {
        final_idx++;
        while (!flag[final_idx])
            final_idx++;
    }

    assert(final_idx<len);
#else
    for (final_col_num=0;final_col_num<n;++final_col_num) {
        if (flagCol[final_col_num]) {
            break;
        }
    }
#endif // 0
    return final_idx;
}

int
get_ref_idx_min(const vector<DComplex> &vals, const vector<bool> &flag)
{
    int min_val_idx(0);
    int len(flag.size());

    for (min_val_idx=0;min_val_idx<len;++min_val_idx) {
        if (flag[min_val_idx]) {
            break;
        }
    }
    if (min_val_idx == len) return -1;
    for (int i=min_val_idx+1;i<len;++i) {
        if (flag[i]) {
            if ( abs(vals[i]) < abs(vals[min_val_idx]) )
                min_val_idx = i;
        }
    }
    if (min_val_idx == len) return -1;
    return min_val_idx;
}

int
get_inter_num(const double length, const double lambda)
{
    int num(0);
    double electric_length(length/lambda);
    if (electric_length < 0.5) {
        double resolution = lambda * 0.05;
        num = (int)(length/resolution) + 1;
    }
    else if (electric_length < 1) {
        double resolution = lambda * 0.1;
        num = (int)(length/resolution) + 1;
    }
    else if (electric_length < 1.5) {
        double resolution = lambda * 0.15;
        num = (int)(length/resolution) + 1;
    }
    else if (electric_length < 2) {
        double resolution = lambda * 0.2;
        num = (int)(length/resolution) + 1;
    }
    else if (electric_length < 3) {
        double resolution = lambda * 0.3;
        num = (int)(length/resolution) + 1;
    }
    else if (electric_length < 4) {
        double resolution = lambda * 0.4;
        num = (int)(length/resolution) + 1;
    }
    else {
        double resolution = lambda * 0.5;
        num = (int)(length/resolution) + 1;
    }
    return num;
}

void
bfs_hmatrix(phmatrix hm, vector< vector<phmatrix> > &info)
{
    queue<phmatrix> hm_queue;
    hm_queue.push(hm);
    int level = 0;
    while (!hm_queue.empty()) {
        vector<phmatrix> &cur_level = info[level-1];
        int n = hm_queue.size();
        for (int i=0;i<n;++i) {
            phmatrix tmp = hm_queue.front();
            if (tmp->r)  {
                if (level > 0) cur_level.push_back(tmp);
            }
            for (int j=0;j<tmp->csons;++j) {
                for (int k=0;k<tmp->rsons;++k) {
                    hm_queue.push(tmp->son[k+j*tmp->rsons]);
                }
            }
            hm_queue.pop();
        }
        level++;
    }
}
