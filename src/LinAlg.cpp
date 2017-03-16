#include "LinAlg.h"

#if 0
void leqSolve_z(DComplex *a, DComplex *b, const int n, const int nb)
{
    int ipiv[n];
    char trans = 'N';
    int info1, info2;
    info1 = LAPACKE_zgetrf(LAPACK_ROW_MAJOR, n, n, (MKL_Complex16 *)a, n, ipiv);
    info2 = LAPACKE_zgetrs(LAPACK_ROW_MAJOR, trans, n, nb, (MKL_Complex16 *)a, n, ipiv, (MKL_Complex16 *)b, nb);
    if (!info1 * !info2 == 0) {
        cout << "linear equation solved error!\n" << endl;
        cout << info1 << endl;
        cout << info2 << endl;
    }
}
#endif // 0
void leqSolve_z(DComplex *a, DComplex *b, int n, int nb)
{
    int ipiv[n];
    char trans = 'N';
    int info1, info2;
    info1 = LAPACKE_zgetrf(LAPACK_COL_MAJOR, n, n, (MKL_Complex16 *)a, n, ipiv);
    info2 = LAPACKE_zgetrs(LAPACK_COL_MAJOR, trans, n, nb, (MKL_Complex16 *)a, n, ipiv, (MKL_Complex16 *)b, n);
    if (!info1 * !info2 == 0) {
        cout << "linear equation solved error!\n" << endl;
        cout << info1 << endl;
        cout << info2 << endl;
    }
}

void mmp_z(DComplex *a, CBLAS_TRANSPOSE atrans, DComplex *b, CBLAS_TRANSPOSE btrans, DComplex *c, int m, int k, int n)
{
	DComplex alpha = 1.0, beta = 0.0;
	int lda = m, ldb = 0, ldc = m;
	if (btrans == CblasNoTrans) ldb = k;
	else ldb = n;
	cblas_zgemm(CblasColMajor,atrans,btrans,m,n,k,&alpha,(MKL_Complex16 *)a,lda,(MKL_Complex16 *)b,ldb,&beta,(MKL_Complex16 *)c,ldc);
}

void matInv_z(DComplex *a, int n)
{
    int ipiv[n];
    int info1, info2;
    info1 = LAPACKE_zgetrf(LAPACK_COL_MAJOR, n, n, (MKL_Complex16 *)a, n, ipiv);
    info2 = LAPACKE_zgetri(LAPACK_COL_MAJOR, n, (MKL_Complex16 *)a, n, ipiv);
    if (!info1 * !info2 == 0) {
        cout << "matrix inverse solved error!\n" << endl;
        cout << n << endl;
        cout << info1 << endl;
        cout << info2 << endl;
    }
}
