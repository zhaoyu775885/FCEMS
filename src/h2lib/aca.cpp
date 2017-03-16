#include "aca.h"
#include <ctime>
#include <fstream>

int aca_times(0);

#if 0
void
mom_aca_plus(const BasisFunc &rbf, phmatrix hm, double wavenum, Intg intg)
{
    cout << aca_times++ << ": use ACA+..." << endl;
    int m = hm->rc->size;
    int n = hm->cc->size;
    prkmatrix rkmat = hm->r;
#if 0
    if (aca_times == 31) {
        ofstream of("./Z_tmp.txt");
        cout << m << ", " << n << endl;
        DComplex *Z = new DComplex [m*n];
        for (int i=0;i<m;++i) {
            for (int j=0;j<n;++j) {
                Z[i+j*m] = intg(rbf[hm->rc->idx[i]], rbf[hm->cc->idx[j]], wavenum);
                of << Z[i+j*m] << " ";
            }
            of << endl;
        }
    }
#endif // 0
    DComplex *R;
    R = new field [m*n];

    double frbNormZ, frbNormZ_2, frbNormU, frbNormV;
    int curRowNum, curColNum;
    char norm = 'F';
    double threshold = 1e-3;

    int flagRow[m], flagCol[n];
    DComplex ref_row_val[n], ref_col_val[m];
    for (int i=0;i<m;i++) {flagRow[i] = 1;} // 1 means have not been touched. 0 means have been touched.
    for (int j=0;j<n;j++) {flagCol[j] = 1;}

/// init R, U1, V1
    int k = 0;

    int ref_col_idx = getRefCol(flagCol, n, k), use_ref_col;
    for (int i=0;i<m;++i) {
        ref_col_val[i] = intg(rbf[hm->rc->idx[i]], rbf[hm->cc->idx[ref_col_idx]], wavenum);
    }
    int ref_row_idx = getRefRow(ref_col_val, flagRow, m, k), use_ref_row;
    for (int j=0;j<n;j++) {
        ref_row_val[j] = intg(rbf[hm->rc->idx[ref_row_idx]], rbf[hm->cc->idx[j]], wavenum);
    }

    curRowNum = get_rowidx_acaplus(ref_col_val, m, flagRow);
    curColNum = get_colidx_acaplus(ref_row_val, n, flagCol);

    if (curColNum >= 0) {
        if (curRowNum>=0) {
            if ( abs(R[curRowNum+ref_col_idx*m]) > abs(R[ref_row_idx+curColNum*m]) ) {
                use_ref_col = 1;
                use_ref_row = 0;
            }
            else {
                use_ref_col = 0;
                use_ref_row = 1;
            }
        }
        else {
            use_ref_col = 0;
            use_ref_row = 1;
        }
    }
    else {
        if (curRowNum>=0) {
            use_ref_col = 1;
            use_ref_row = 0;
        }
        else {
            use_ref_col = use_ref_row = 0;
            cout << "Full Zeros?!" << endl;
            return ;
        }
    }

    rkmat->V.a = new DComplex [n];
    rkmat->U.a = new DComplex [m];

    if ( use_ref_col ) {
        for (int j=0;j<n;++j) {
            R[curRowNum+j*m] = intg(rbf[hm->rc->idx[curRowNum]], rbf[hm->cc->idx[j]], wavenum);
        }
        flagRow[curRowNum] = 0;
        curColNum = getColACA(R, m, n, curRowNum, flagCol);
        if ( curColNum==-1 ) {
            cout << "Ref Col: Full Zeroes!" << endl;
            assert(0);
        }
        else {
            for (int i=0;i<m;++i) {
                R[i+curColNum*m] = intg(rbf[hm->rc->idx[i]], rbf[hm->cc->idx[curColNum]], wavenum);
            }
        }
    }
    if ( use_ref_row ) {
        for (int i=0;i<m;++i) {
            R[i+curColNum*m] = intg(rbf[hm->rc->idx[i]], rbf[hm->cc->idx[curColNum]], wavenum);
        }
        flagCol[curColNum] = 0;
        curRowNum = getRowACA(R, m, n, curColNum, flagRow);
        if (curRowNum == -1) {
            cout << "Ref Row: Full Zeroes!" << endl;
            assert(0);
        }
        else {
            for (int j=0;j<n;++j) {
                R[curRowNum+j*m] = intg(rbf[hm->rc->idx[curRowNum]], rbf[hm->cc->idx[j]], wavenum);
            }
        }
    }

    assert( abs(R[curRowNum+curColNum*m])>1e-20 );
    for (int j=0;j<n;j++) rkmat->V.a[j] = R[curRowNum+j*m]/R[curRowNum+curColNum*m];
    for (int i=0;i<m;i++) rkmat->U.a[i] = R[i+curColNum*m];

    frbNormV = LAPACKE_zlange(LAPACK_COL_MAJOR, norm, n, 1, (MKL_Complex16 *)rkmat->V.a, n);
    frbNormU = LAPACKE_zlange(LAPACK_COL_MAJOR, norm, m, 1, (MKL_Complex16 *)rkmat->U.a, m);
    frbNormZ = frbNormU*frbNormV;
    frbNormZ_2 = frbNormZ*frbNormZ;

/// init R U1 V1, done!

    while ( frbNormU*frbNormV > threshold*frbNormZ ) {
        k++;

        if (curColNum == ref_col_idx) {
            ref_col_idx = getRefCol(flagCol, n, k);
            if ( ref_col_idx!=-1 ) {
                for (int i=0;i<m;++i) {
                    ref_col_val[i] = intg(rbf[hm->rc->idx[i]], rbf[hm->cc->idx[ref_col_idx]], wavenum);
                    for (int j=0;j<k;j++) ref_col_val[i] -= rkmat->U.a[i+j*m] * rkmat->V.a[ref_col_idx+j*n];
                }
            }
        }
        else {
            for (int i=0;i<m;++i) {
                for (int j=0;j<k;j++) ref_col_val[i] -= rkmat->U.a[i+j*m]*rkmat->V.a[ref_col_idx+j*n];
            }
        }

        if (curRowNum == ref_row_idx) {
            ref_row_idx = getRefRow(ref_col_val, flagRow, m, k);
            if ( ref_row_idx!=-1 ) {
                for (int j=0;j<n;++j) {
                    ref_row_val[j] = intg(rbf[hm->rc->idx[ref_row_idx]], rbf[hm->cc->idx[j]], wavenum);
                    for (int i=0;i<k;i++) ref_row_val[j] -= rkmat->U.a[ref_row_idx+i*m]*rkmat->V.a[j+i*n];
                }
            }
        }
        else {
            for (int j=0;j<n;++j) {
                for (int i=0;i<k;i++) ref_row_val[j] -= rkmat->U.a[ref_row_idx+i*m]*rkmat->V.a[j+i*n];
            }
        }

        curRowNum = get_rowidx_acaplus(ref_col_val, m, flagRow);
        curColNum = get_colidx_acaplus(ref_row_val, n, flagCol);

        if (curColNum >= 0) {
            if (curRowNum>=0) {
                if ( abs(R[curRowNum+ref_col_idx*m]) > abs(R[ref_row_idx+curColNum*m]) ) {
                    use_ref_col = 1;
                    use_ref_row = 0;
                }
                else {
                    use_ref_col = 0;
                    use_ref_row = 1;
                }
            }
            else {
                use_ref_col = 0;
                use_ref_row = 1;
            }
        }
        else {
            if (curRowNum>=0) {
                use_ref_col = 1;
                use_ref_row = 0;
            }
            else {
                cout << "Full Zeros?!" << endl;
                return;
            }
        }

        if ( use_ref_col ) {
            #pragma omp parallel for
            for (int j=0;j<n;++j) {
                R[curRowNum+j*m] = intg(rbf[hm->rc->idx[curRowNum]], rbf[hm->cc->idx[j]], wavenum);
                for (int i=0;i<k;i++) R[curRowNum+j*m] -= rkmat->U.a[curRowNum+i*m]*rkmat->V.a[j+i*n];
            }
            flagRow[curRowNum] = 0;
            curColNum = getColACA(R, m, n, curRowNum, flagCol);
            if (curColNum == -1) {
                k--;
                break;
            }
            else {
                #pragma omp parallel for
                for (int i=0;i<m;++i) {
                    R[i+curColNum*m] = intg(rbf[hm->rc->idx[i]], rbf[hm->cc->idx[curColNum]], wavenum);
                    for (int j=0;j<k;j++) R[i+curColNum*m] -= rkmat->U.a[i+j*m] * rkmat->V.a[curColNum+j*n];
                }
            }
        }
        if ( use_ref_row ) {
            #pragma omp parallel for
            for (int i=0;i<m;++i) {
                R[i+curColNum*m] = intg(rbf[hm->rc->idx[i]], rbf[hm->cc->idx[curColNum]], wavenum);
                for (int j=0;j<k;j++) R[i+curColNum*m] -= rkmat->U.a[i+j*m] * rkmat->V.a[curColNum+j*n];
            }
            flagCol[curColNum] = 0;
            curRowNum = getRowACA(R, m, n, curColNum, flagRow);
            if ( curRowNum==-1 ) {
                k--;
                break;
            }
            else {
                #pragma omp parallel for
                for (int j=0;j<n;++j) {
                    R[curRowNum+j*m] = intg(rbf[hm->rc->idx[curRowNum]], rbf[hm->cc->idx[j]], wavenum);
                    for (int i=0;i<k;i++) R[curRowNum+j*m] -= rkmat->U.a[curRowNum+i*m]*rkmat->V.a[j+i*n];
                }
            }
        }

//        cout << ref_col_idx << ", " << ref_row_idx << endl;
//        cout << curColNum << ", " << curRowNum << endl;
//        cout << k << endl;

        assert(abs(R[curRowNum+curColNum*m]) > 1e-20);

        rkmat->V.a = (DComplex *) realloc(rkmat->V.a, (k+1)*n*sizeof(DComplex));
        rkmat->U.a = (DComplex *) realloc(rkmat->U.a, (k+1)*m*sizeof(DComplex));
        for (int j=0;j<n;j++) rkmat->V.a[j+k*n] = R[curRowNum+j*m]/R[curRowNum+curColNum*m];
        for (int i=0;i<m;i++) rkmat->U.a[i+k*m] = R[i+curColNum*m];

        frbNormV = LAPACKE_zlange(LAPACK_COL_MAJOR, norm, n, 1, (MKL_Complex16 *)rkmat->V.a+k*n, n);
        frbNormU = LAPACKE_zlange(LAPACK_COL_MAJOR, norm, m, 1, (MKL_Complex16 *)rkmat->U.a+k*m, m);

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
//	cout << "m=" << m << ", n=" << n << ", rank=" << rkmat->k << endl;
    free(R);
}

void
mom_aca(const BasisFunc &rbf, phmatrix hm, double wavenum, Intg intg)
{
    cout << aca_times++ << ": use ACA" << endl;
    int m = hm->rc->size;
    int n = hm->cc->size;

    DComplex *R, *Z;

    int *flagRow, *flagCol;
    double frbNormZ, frbNormZ_2, frbNormU, frbNormV;
    int curRowNum, curColNum;

    char norm = 'F';
    double threshold = 1e-3;
    prkmatrix rkmat = hm->r;

    R = new field [m*n];
    Z = new field [m*n];
    flagRow = new int[m];
    flagCol = new int[n];

    for (int i=0;i<m;i++) {flagRow[i] = 1;} // 1 means have not been touched. 0 means have been touched.
    for (int j=0;j<n;j++) {flagCol[j] = 1;}

/// init R, U1, V1
    int k = 0;

    curRowNum = 0;
    flagRow[curRowNum] = 0;
    for (int j=0;j<n;j++) {
        Z[curRowNum+j*m] = intg(rbf[hm->rc->idx[curRowNum]], rbf[hm->cc->idx[j]], wavenum);
        R[curRowNum+j*m]= Z[curRowNum+j*m];
    }
    curColNum = getColACA(R, m, n, curRowNum, flagCol);
    /// The selected column may have invalid value, like -1.
    for (int i=0;i<m;i++) {
        Z[i+curColNum*m] = intg(rbf[hm->rc->idx[i]], rbf[hm->cc->idx[curColNum]], wavenum);
        R[i+curColNum*m] = Z[i+curColNum*m];
    }

    assert( abs(R[curRowNum+curColNum*m])>1e-20 );
    for (int j=0;j<n;j++) rkmat->V.a[j] = R[curRowNum+j*m]/R[curRowNum+curColNum*m];
    for (int i=0;i<m;i++) rkmat->U.a[i] = R[i+curColNum*m];

    frbNormV = LAPACKE_zlange(LAPACK_COL_MAJOR, norm, n, 1, (MKL_Complex16 *)rkmat->V.a, n);
    frbNormU = LAPACKE_zlange(LAPACK_COL_MAJOR, norm, m, 1, (MKL_Complex16 *)rkmat->U.a, m);
    frbNormZ = frbNormU * frbNormV;
    frbNormZ_2 = frbNormZ * frbNormZ;
/// init R U1 V1, done!

    while ( frbNormU*frbNormV > threshold*frbNormZ ) {
        k++;

        if((curRowNum = getRowACA(R, m, n, curColNum, flagRow)) == -1) {
            k--;
            break;
        }
        #pragma omp parallel for
        for (int j=0;j<n;j++) {
            Z[curRowNum+j*m] = intg(rbf[hm->rc->idx[curRowNum]], rbf[hm->cc->idx[j]], wavenum);
            R[curRowNum+j*m] = Z[curRowNum+j*m];
            for (int i=0;i<k;i++) {
                R[curRowNum+j*m] -= rkmat->U.a[curRowNum+i*m]*rkmat->V.a[j+i*n];
            }
        }
        if ((curColNum = getColACA(R, m, n, curRowNum, flagCol)) == -1) {
            k--;
            break;
        }
        #pragma omp parallel for
        for (int i=0;i<m;i++) {
            Z[i+curColNum*m] = intg(rbf[hm->rc->idx[i]], rbf[hm->cc->idx[curColNum]], wavenum);
            R[i+curColNum*m] = Z[i+curColNum*m];
            for (int j=0;j<k;j++) {
                R[i+curColNum*m] -= rkmat->U.a[i+j*m] * rkmat->V.a[curColNum+j*n];
            }
        }

        rkmat->V.a = (DComplex *) realloc(rkmat->V.a, (k+1)*n*sizeof(DComplex));
        for (int j=0;j<n;j++) {
            rkmat->V.a[j+k*n] = R[curRowNum+j*m]/R[curRowNum+curColNum*m];
        }
        rkmat->U.a = (DComplex *) realloc(rkmat->U.a, (k+1)*m*sizeof(DComplex));
        for (int i=0;i<m;i++) {
            rkmat->U.a[i+k*m] = R[i+curColNum*m];
        }

        frbNormV = LAPACKE_zlange(LAPACK_COL_MAJOR, norm, n, 1, (MKL_Complex16 *)rkmat->V.a+k*n, n);
        frbNormU = LAPACKE_zlange(LAPACK_COL_MAJOR, norm, m, 1, (MKL_Complex16 *)rkmat->U.a+k*m, m);

        for (int i=0;i<k;i++) {
            DComplex c1, c2;
            cblas_zdotc_sub(m, rkmat->U.a+i*m, 1, rkmat->U.a+k*m, 1, &c1);
            cblas_zdotc_sub(n, rkmat->V.a+i*n, 1, rkmat->V.a+k*n, 1, &c2);
            frbNormZ_2 += 2*abs( c1 )*abs( c2 );
        }
        frbNormZ_2 += frbNormU*frbNormU*frbNormV*frbNormV;
        frbNormZ = sqrt(frbNormZ_2);
    }

    rkmat->k = k+1;
    rkmat->U.rows = m;
    rkmat->U.cols = rkmat->k;
    rkmat->V.rows = n;
	rkmat->V.cols = rkmat->k;

	conj_amatrix(&(rkmat->V));
//	cout << "m=" << m << ", n=" << n << ", rank=" << rkmat->k << endl;

	delete [] Z;
    delete [] R;
    delete [] flagRow;
    delete [] flagCol;
}

#if 0
void
mom_aca_error(const BasisFunc &rbf, phmatrix hm, double wavenum, Intg intg)
{
    cout << "use ACA..." << endl;
    int m = hm->rc->size;
    int n = hm->cc->size;

    DComplex *R, *Z;

    int *flagRow, *flagCol;

    double frbNormZ, frbNormZ_2, frbNormU, frbNormV;
    int curRowNum, curColNum;

    char norm = 'F';

    double threshold = 1e-3;

    //hm->r = new_rkmatrix(m, n, 1);
    prkmatrix rkmat = hm->r;


    R = new field [m*n];
    Z = new field [m*n];

    flagRow = new int[m];
    flagCol = new int[n];

    for (int i=0;i<m;i++) {
        flagRow[i] = 1; // 1 means have not been touched. 0 means have been touched.
    }
    for (int j=0;j<n;j++) {
        flagCol[j] = 1;
    }

/// init R, U1, V1
    int k = 0;
    curRowNum = 0;
    for (int j=0;j<n;j++) {
        // [curRowNum, j]
        Z[curRowNum*n+j] = intg(rbf[hm->rc->idx[curRowNum]], rbf[hm->cc->idx[j]], wavenum);
        R[curRowNum*n+j]= Z[curRowNum*n+j];
    }
    curColNum = getColACA(R, m, n, curRowNum, flagCol);

    for (int j=0;j<n;j++) {
        rkmat->V.a[j] = R[curRowNum*n+j]/R[curRowNum*n+curColNum];
    }
    for (int i=0;i<m;i++) {
        // [i, curColNum]
        Z[i*n+curColNum] = intg(rbf[hm->rc->idx[i]], rbf[hm->cc->idx[curColNum]], wavenum);
        R[i*n+curColNum] = Z[i*n+curColNum];
    }

    for (int i=0;i<m;i++) {
        rkmat->U.a[i] = R[i*n+curColNum];
    }

    frbNormV = LAPACKE_zlange(LAPACK_COL_MAJOR, norm, n, 1, (MKL_Complex16 *)rkmat->V.a, n);
    frbNormU = LAPACKE_zlange(LAPACK_COL_MAJOR, norm, m, 1, (MKL_Complex16 *)rkmat->U.a, m);

    frbNormZ_2 = frbNormU*frbNormU * frbNormV*frbNormV;
    frbNormZ = sqrt(frbNormZ_2);

/// init R U1 V1, done!

    while ( frbNormU*frbNormV > threshold*frbNormZ ) {

        k++;

        if((curRowNum = getRowACA(R, m, n, curColNum, flagRow)) == -1) {
            k--;
//            printf("quit now, k = %d\n", k);
            break;
        }

        #pragma omp parallel for
        for (int j=0;j<n;j++) {
            Z[curRowNum*n+j] = intg(rbf[hm->rc->idx[curRowNum]], rbf[hm->cc->idx[j]], wavenum);
            R[curRowNum*n+j] = Z[curRowNum*n+j];
            for (int i=0;i<k;i++) {
                R[curRowNum*n+j] -= rkmat->U.a[i*m+curRowNum]*rkmat->V.a[i*n+j];
            }
        }
        if ((curColNum = getColACA(R, m, n, curRowNum, flagCol)) == -1) {
            k--;
//            printf("quit now, k = %d\n", k);
            break;
        }

        rkmat->V.a = (DComplex *) realloc(rkmat->V.a, (k+1)*n*sizeof(DComplex));
        for (int j=0;j<n;j++) {
            rkmat->V.a[j+k*n] = R[curRowNum*n+j]/R[curRowNum*n+curColNum];
        }

        #pragma omp parallel for
        for (int i=0;i<m;i++) {
            Z[i*n+curColNum] = intg(rbf[hm->rc->idx[i]], rbf[hm->cc->idx[curColNum]], wavenum);
            R[i*n+curColNum] = Z[i*n+curColNum];
            for (int j=0;j<k;j++) {
                R[i*n+curColNum] -= rkmat->V.a[j*n+curColNum]*rkmat->U.a[j*m+i];
            }
        }

        rkmat->U.a = (DComplex *) realloc(rkmat->U.a, (k+1)*m*sizeof(DComplex));
        for (int i=0;i<m;i++) {
            rkmat->U.a[k*m+i] = R[i*n+curColNum];
        }

        frbNormV = LAPACKE_zlange(LAPACK_COL_MAJOR, norm, n, 1, (MKL_Complex16 *)rkmat->V.a+k*n, n);
        frbNormU = LAPACKE_zlange(LAPACK_COL_MAJOR, norm, m, 1, (MKL_Complex16 *)rkmat->U.a+k*m, m);

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

    free(Z);
    free(R);
    free(flagRow);
    free(flagCol);
}

prkmatrix
aca(DComplex *Z, int m, int n)
{
    int curRowNum, curColNum;
    DComplex *R;
    int *flagRow, *flagCol;
    double frbNormZ, frbNormZ_2, frbNormU, frbNormV;
    char norm = 'f';
    double threshold = 1e-4;

    prkmatrix rkmat = new_rkmatrix(m, n, 1);


    R = (DComplex *)malloc(m*n*sizeof(DComplex));
    flagRow = (int *) malloc(m*sizeof(int));
    flagCol = (int *) malloc(n*sizeof(int));

    for (int i=0;i<m;i++) {
        flagRow[i] = 1; // 1 means have not been touched. 0 means have been touched.
    }
    for (int j=0;j<n;j++) {
        flagCol[j] = 1;
    }


/// init R, U1, V1
    int k = 0;
    curRowNum = 0;
    for (int j=0;j<n;j++) {
        // [curRowNum, j]
        R[curRowNum*n+j]= Z[curRowNum*n+j];
    }
    curColNum = getColACA(R, m, n, curRowNum, flagCol);

//    rkmat->V.a = (datatype *) malloc(n*sizeof(datatype));
    for (int j=0;j<n;j++) {
//
        rkmat->V.a[j] = R[curRowNum*n+j]/R[curRowNum*n+curColNum];
    }
    for (int i=0;i<m;i++) {
        // [i, curColNum]
        R[i*n+curColNum] = Z[i*n+curColNum];
    }

//    rkmat->U.a = (datatype *) malloc(m*sizeof(datatype));
    for (int i=0;i<m;i++) {
//
        rkmat->U.a[i] = R[i*n+curColNum];
    }

//    U = (DComplex *) realloc(U, m*sizeof(DComplex ));
//    printf("ok?\n");

//    frbNormV = LAPACKE_dlange(LAPACK_ROW_MAJOR, norm, 1, n, V, n);
    frbNormV = LAPACKE_zlange(LAPACK_ROW_MAJOR, norm, 1, n, (MKL_Complex16 *)rkmat->V.a, n);
//    frbNormU = LAPACKE_dlange(LAPACK_COL_MAJOR, norm, m, 1, U, m);
    frbNormU = LAPACKE_zlange(LAPACK_COL_MAJOR, norm, m, 1, (MKL_Complex16 *)rkmat->U.a, m);

    frbNormZ_2 = frbNormU*frbNormU * frbNormV*frbNormV;
    frbNormZ = sqrt(frbNormZ_2);

/// init R U1 V1, done!

    while ( frbNormU*frbNormV > threshold*frbNormZ ) {

        k++;

        if((curRowNum = getRowACA(R, m, n, curColNum, flagRow)) == -1) {
            k--;
//            printf("quit now, k = %d\n", k);
            break;
        }

        for (int j=0;j<n;j++) {
            R[curRowNum*n+j] = Z[curRowNum*n+j];
            for (int i=0;i<k;i++) {
                R[curRowNum*n+j] -= rkmat->U.a[i*m+curRowNum]*rkmat->V.a[i*n+j];
            }
        }
        if ((curColNum = getColACA(R, m, n, curRowNum, flagCol)) == -1) {
            k--;
//            printf("quit now, k = %d\n", k);
            break;
        }

        rkmat->V.a = (DComplex *) realloc(rkmat->V.a, (k+1)*n*sizeof(DComplex ));
        for (int j=0;j<n;j++) {
//
            rkmat->V.a[k*n+j] = R[curRowNum*n+j]/R[curRowNum*n+curColNum];
        }
        for (int i=0;i<m;i++) {
            R[i*n+curColNum] = Z[i*n+curColNum];
            for (int j=0;j<k;j++) {
                R[i*n+curColNum] -= rkmat->V.a[j*n+curColNum]*rkmat->U.a[j*m+i];
            }
        }

        rkmat->U.a = (DComplex *) realloc(rkmat->U.a, (k+1)*m*sizeof(DComplex ));
        for (int i=0;i<m;i++) {
//
            rkmat->U.a[k*m+i] = R[i*n+curColNum];
        }

//        frbNormV = LAPACKE_dlange(LAPACK_ROW_MAJOR, norm, 1, n, V+k*n, n);
//        frbNormU = LAPACKE_dlange(LAPACK_COL_MAJOR, norm, m, 1, U+k*m, m);
        frbNormV = LAPACKE_zlange(LAPACK_ROW_MAJOR, norm, 1, n, (MKL_Complex16 *)rkmat->V.a+k*n, n);
        frbNormU = LAPACKE_zlange(LAPACK_COL_MAJOR, norm, m, 1, (MKL_Complex16 *)rkmat->U.a+k*m, m);

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
    rkmat->V.rows = rkmat->k;
    rkmat->V.cols = n;

    free(R);
    free(flagRow);
    free(flagCol);

    return rkmat;
}
#endif // 0
#endif // 0

