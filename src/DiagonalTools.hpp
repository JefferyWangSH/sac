/**
 *  This source file includes some diagonalizing tools with C++/Eigen interface
 *  for diagonalizing real matrices using mkl and lapack.
 *  including:
 *    1. generalized SVD decomposition for arbitrary M * N matrices,
 *    2. optimized diagonalizing mechanism for N * N real symmetric matrix
 *  and calculation accuracy is guaranteed.
 */


#include <iostream>

#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#include <Eigen/Core>
#include "mkl_lapacke.h"


/**
 * SVD decomposition of arbitrary M * N real matrix, using MKL_LAPACK:
 *      A  ->  U * S * V^T
 * Remind that V is returned in this subroutine, not V transpose.
 *
 * @param m -> number of rows.
 * @param n -> number of cols.
 * @param a -> arbitrary M*N real matrix to be solved.
 * @param u -> u matrix in Eigen::Matrix, M * M.
 * @param s -> eigenvalues s in Eigen::Vector, descending sorted.
 * @param v -> v matrix in Eigen::Matrix, N * N.
 */
void mkl_lapack_dgesvd(const int &m, const int &n, const Eigen::MatrixXd &a, Eigen::MatrixXd &u, Eigen::VectorXd &s, Eigen::MatrixXd &v) {
    assert(m == a.rows());
    assert(n == a.cols());

    // Matrix size
    int matrix_layout = LAPACK_ROW_MAJOR;
    lapack_int info, lda = m, ldu = m, ldvt = n;

    // Local arrays
    double s_[ldu * ldu], u_[ldu * m], vt_[ldvt * n];
    double a_[lda * n];
    double superb[ldu * lda];
    for (int i = 0; i < lda * n; ++i) {
        a_[i] = a(i / lda, i % lda);
    }

    // Compute SVD
    info = LAPACKE_dgesvd( matrix_layout, 'A', 'A', m, n, a_, lda, s_, u_, ldu, vt_, ldvt, superb );

    // Check for convergence
    if( info > 0 ) {
        std::cerr << "The algorithm computing SVD failed to converge." << std::endl;
        exit( 1 );
    }

    // Convert results to Eigen
    u = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(u_, n, n);
    s = Eigen::Map<Eigen::Matrix<double, 1, Eigen::Dynamic, Eigen::RowMajor>>(s_, 1, n);
    v = Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>>(vt_, m, m);
}


/**
 * Calculate eigenvalues and eigenstates given an arbitrary N * N real symmetric matrix, using MKL_LAPACK
 *      A  ->  T^dagger * S * T
 * where T is rotation matrix, which is orthogonal;
 *       S is diagonal matrix with eigenvalues being diagonal elements.
 *
 * @param n -> number of rows/cols.
 * @param a -> arbitrary N * N real symmetric matrix to be solved.
 * @param s -> diagonal eigen matrix, descending sorted.
 * @param T -> rotation matrix, columns being eigenstates.
 */
void mkl_lapack_dsyev(const int &n, const Eigen::MatrixXd &a, Eigen::MatrixXd &s, Eigen::MatrixXd &T) {
//    LAPACKE_dsyev();
}