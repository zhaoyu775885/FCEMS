#ifndef POSTPROCESS_H
#define POSTPROCESS_H

#include "settings.h"
#include "BasisFunction.h"
#include "LinAlg.h"

class PostProcess
{
public:
    PostProcess(const BasisFunc &func_space, const DComplex *src,
                const int r, const int c, const double f, const char ie_type) :
        basis_space(func_space), data(src), rows(r), cols(c),
        freq(f), deg_vec(0), rcs_vec(0), type(ie_type) {};

    void print_data() const { for(int i=0;i<rows;++i) for(int j=0;j<cols;++j) cout << data[i+j*cols] << endl; }
    void gen_rcs(const double beg, const double end, const int count);

    void write_rcs(const string &filename);
    void write_mat(const string &filename, const char c) const;

private:
    const BasisFunc &basis_space;
    const DComplex *data;
    /**
    rows, cols: The original data form of the data.
    For further application, one should make use of
    part or the whole of it to obtain the specific
    parameters, RCS, S-mat, etc.
     **/
    int rows;
    int cols;
    double freq;
    vector<double> deg_vec;
    vector<double> rcs_vec;

    const char type;

};

#endif // POSTPROCESS_H
