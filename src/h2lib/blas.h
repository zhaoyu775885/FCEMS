/* ------------------------------------------------------------
 This is the file "blas.h" of the H2Lib package.
 All rights reserved, Sven Christophersen 2015
 ------------------------------------------------------------ */

/**
 * @file blas.h
 * @author Sven Christophersen
 * @date 2015
 */

#ifndef BLAS_H_
#define BLAS_H_


/* C STD LIBRARY */
/* CORE 0 */
#include "h2settings.h"
/* CORE 1 */
/* CORE 2 */
/* CORE 3 */
/* SIMPLE */
/* PARTICLES */
/* BEM */

#define h2_dot(n, x, incx, y, incy, res) cblas_zdotc_sub(n, x, incx, y, incy, res)


#define h2_rdot(n, x, incx, y, incy) cblas_ddot(n, x, incx, y, incy)


#define h2_axpy(n, alpha, x, incx, y, incy) cblas_zaxpy(n, alpha, x, incx, y, incy)

#define h2_raxpy(n, alpha, x, incx, y, incy) cblas_zaxpy(n, alpha, x, incx, y, incy)

#define h2_scal(n, alpha, x, incx) cblas_zscal(n, alpha, x, incx)


#define h2_rscal(n, alpha, x, incx) cblas_zdscal(n, alpha, x, incx)

#define h2_nrm2(n, x, incx) cblas_dznrm2(n, x, incx)

#define h2_swap(n, x, incx, y, incy) cblas_zswap(n, x, incx, y, incy)


/****************************************************
 * BLAS level 2
 ****************************************************/

/****************************************************
 * _GEMV
 ****************************************************/

#define h2_gemv(layout, trans, m, n, alpha, a, lda, x, incx, beta, y, incy) cblas_zgemv(layout, trans, m, n, alpha, a, lda, x, incx, beta, y, incy)

#define h2_trmv(layout, uplo, trans, diag, n, a, lda, x, incx) cblas_ztrmv(layout, uplo, trans, diag, n, a, lda, x, incx)

/****************************************************
 * _GER
 ****************************************************/

#define h2_ger(layout, m, n, alpha, x, incx, y, incy, a, lda) cblas_zgerc(layout, m, n, alpha, x, incx, y, incy, a, lda)

#define h2_gerc(layout, m, n, alpha, x, incx, y, incy, a, lda) cblas_zgerc(layout, m, n, alpha, x, incx, y, incy, a, lda)

#define h2_geru(layout, m, n, alpha, x, incx, y, incy, a, lda) cblas_zgeru(layout, m, n, alpha, x, incx, y, incy, a, lda)


/****************************************************
 * _SYR
 ****************************************************/

#define h2_syr(layout, uplo, n, alpha, x, incx, a, lda) cblas_zher(layout, uplo, n, alpha, x, incx, a, lda)


/****************************************************
 * BLAS level 3
 ****************************************************/

/****************************************************
 * _GEMM
 ****************************************************/

#define h2_gemm(layouta, transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc) cblas_zgemm(layouta, transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)


/****************************************************
 * _TRMM
 ****************************************************/

#define h2_trmm(layout, side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb) cblas_ztrmm(layout, side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)

/****************************************************
 * _TRSM
 ****************************************************/

#define h2_trsm(layout, side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb) cblas_ztrsm(layout, side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)


/****************************************************
 * LAPACK
 ****************************************************/

/****************************************************
 * _LACGV
 ****************************************************/

#define h2_lacgv(n, x, incx) cblas_zlacgv(n, x, incx)


/****************************************************
 * _POTRF
 ****************************************************/

#define h2_potrf(layout, uplo, n, a, lda) LAPACKE_zpotrf(layout, uplo, n, a, lda)


/****************************************************
 * _LARF / _LARFG
 ****************************************************/

#define h2_larf(side, m, n, v, incv, tau, c, ldc, work) zlarf(side, m, n, v, incv, tau, c, ldc, work)

#define h2_larfg(n, alpha, x, incx, tau) LAPACKE_zlarfg(n, alpha, x, incx, tau)


/****************************************************
 * _GEQRF
 ****************************************************/

#define h2_geqrf(layout, m, n, a, lda, tau) LAPACKE_zgeqrf(layout, m, n, a, lda, tau)


/****************************************************
 * _ORMQR
 ****************************************************/

#define h2_ormqr(layout, side, trans, m, n, k, a, lda, tau, c, ldc) LAPACKE_zunmqr(layout, side, trans, m, n, k, a, lda, tau, c, ldc)


/****************************************************
 * _ORGQR
 ****************************************************/

#define h2_orgqr(layout, m, n, k, a, lda, tau) LAPACKE_zungqr(layout, m, n, k, (MKL_Complex16 *)a, lda, tau)


/****************************************************
 * _STEQR
 ****************************************************/

#define h2_steqr(layout, compz, n, d, e, z, ldz) LAPACKE_zsteqr(layout, compz, n, d, e, z, ldz)


/****************************************************
 * _STEV
 ****************************************************/

#define h2_stev(layout, jobz, n, d, e, z, ldz) LAPACKE_dstev(layout, jobz, n, d, e, z, ldz)

/****************************************************
 * _GESVD
 ****************************************************/

#define h2_gesvd(layout, jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, rwork) LAPACKE_zgesvd(layout, jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, rwork)


/****************************************************
 * _DBSQR
 ****************************************************/

#define h2_bdsqr(layout, uplo, n, ncvt, nru, ncc, d, e, vt, ldvt, u, ldu, c, ldc) LAPACKE_zbdsqr(layout, uplo, n, ncvt, nru, ncc, d, e, vt, ldvt, u, ldu, c, ldc)

/****************************************************
 * _SYEV / _HEEV
 ****************************************************/

#define h2_heev(layout, jobz, uplo, n, a, lda, w) LAPACKE_zheev(layout, jobz, uplo, n, a, lda, w)


#endif

