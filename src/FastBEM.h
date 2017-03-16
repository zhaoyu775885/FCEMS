#ifndef FASTBEM_H
#define FASTBEM_H

#include "settings.h"
#include "h2lib/cluster.h"
#include "h2lib/block.h"
#include "h2lib/hmatrix.h"
#include "h2lib/harith.h"
#include "h2lib/truncation.h"
#include "Intg.h"
#include "LinAlg.h"

#include "h2lib/matrixplot.h"

typedef struct {
    int max_leaf_size;
    double eta;
    double threshold;
} FastBemConf;

class FastBem
{
public:
    FastBem(const BasisFunc &func_space, const FastBemConf &config);

    int get_nrow() const {return nrow;}
    int get_nrhs() const {return nrhs;}
    const DComplex *export_data() const {return uvec->v;}
    void set_cluster(pccluster pc);

    void build_equation(const double freq, const char c='e');
    void build_system_matrix(const double freq, const char c='e');
    void build_excitation(const double freq, const char c='e');

    void build_main_hmatrix();
    phmatrix get_hmatrix() const {return g_hmatrix;}

    void assemble(const Integration &intg, const char c);
    void direct_solve_vec();
    void iterative_solve_vec();
    void print_zmat(const string filename) const;
    void print_umat(const string filename) const;

    phmatrix assemble_hmatrix(phmatrix hm, const Integration &intg, const char c);
    void assemble_hmatrix(phmatrix hm0, phmatrix hm1, phmatrix hm2, phmatrix hm3, const Integration &intg);
    phmatrix assemble_hmatrix_fmatrix(phmatrix hm, const Integration &intg, const char c);
    phmatrix assemble_hmatrix_rkmatrix(phmatrix hm, const Integration &intg, const char c);

    phmatrix assemble_rkmatrix(phmatrix hm, const Integration &intg, const char c);
    void assemble_rkmatrix(phmatrix hm0, phmatrix hm1, phmatrix hm2, phmatrix hm3, const Integration &intg);
    phmatrix assemble_fmatrix(phmatrix hm, const Integration &intg, const char c);
    void assemble_fmatrix(phmatrix hm0, phmatrix hm1, phmatrix hm2, phmatrix hm3, const Integration &intg);
    void aca_plus(phmatrix hm, const Integration &intg, const char c);
    void aca(phmatrix hm, const Integration &intg, const char c);
    void hca(phmatrix hm, const Integration &intg, const char c);
    void hca(phmatrix hm0, phmatrix hm1, phmatrix hm2, phmatrix hm3, const Integration &intg);

    void hca_naive(phmatrix hm, const Integration &intg, const char c);
    void hca_std(phmatrix hm, const Integration &intg, const char c);
    void hca_m(phmatrix hm, const Integration &intg, const char c);

    void cost_reduction(const double eps);
    void cost_reduction_hmat(phmatrix hm, const double eps);
    double get_rel_error() const {return frob_norm_error_2_rel;}
    int get_max_rank() const {return max_rank_in_rk;}
    void precondition();

private:
    const BasisFunc &func_space;

    pcluster g_cluster;
    phmatrix g_hmatrix;
    pavector uvec;
    int *index_tree;
    int *aux_index_tree;
    pccluster aux_cluster;

    int nrow;
    int nrhs;
    int leaf_size_threshold;
    double kesi;

    int aca_plus_times=0;
    int aca_times=0;
    int hca_times=0;
    double aca_threshold;
    double svd_threshold;

    int max_rank_in_rk;
    double frob_norm_2_total;
    double frob_norm_error_2_abs;
    double frob_norm_error_2_rel;
};


void
bfs_hmatrix(phmatrix hm, vector< vector<phmatrix> > &info);

int
get_next_idx_aca(const vector<DComplex> &vals, vector<bool> &flag);

int
get_next_idx_aca_plus(const vector<DComplex> &vals, const vector<bool> &flag);

int
get_ref_idx_rand(const vector<bool> &flag, const int cur_rank);

int
get_ref_idx_min(const vector<DComplex> &vals, const vector<bool> &flag);

int
full_col_pivot_aca(const DComplex *A, vector<int> &row_flag, vector<int> &col_flag);

int
get_inter_num(const double length, const double lambda);

#endif // FASTBEM_H
