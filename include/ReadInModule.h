#ifndef SAC_READINMODULE_H
#define SAC_READINMODULE_H
#pragma once

/**
  *   Interface class for the read-in of QMC data
  *   Features:
  *     1. read-in of imaginary time sequence tau, dynamical correlations G(tau)
  *        and their QMC errors in terms of bin.
  *     2. generate and diagonalize covariance matrix, further compute its eigenvalue and eigenvector for SAC use.
  *        ( The SAC calculation is performed in eigen space of the covariance matrix. )
  */

/*
 * TODO：
 *   1. use pointers for memory management
 *   2. discard correlations with poor data quality (done)
 *   3. optimize diagonalizing process, using dsyev (done)
 */

#include <string>

#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#include <Eigen/Core>


namespace QMCData {
    class ReadInModule {
    public:
        int lt{};                       // number of imaginary time slices
        int nbin{};                     // number of bins (after rebin)
        int nbin_total{};               // number of bins (before rebin)
        int rebin_pace{};               // rebin factor (pace of rebin), 1 for all bins
        int num_bootstrap{};            // number of bootstrap samples
        double beta{};                  // inverse temperature beta
        double g0{};                    // static correlation at tau = 0, worked as normalization parameter

        int cov_mat_dim{};              // dimension of covariance matrix
        Eigen::MatrixXd cov_mat;        // covariance matrix
        Eigen::MatrixXd rotate_mat;     // orthogonal matrix, rotating covariance matrix to its eigen space
        Eigen::VectorXd cov_eig;        // eigenvalues of covariance matrix

        Eigen::VectorXd tau, corr_mean, corr_err;

        // helping matrices
        Eigen::MatrixXd corr_tau_bin;
        Eigen::MatrixXd sample_bootstrap;


    public:

        ReadInModule() = default;

        void deallocate_memory();

        /* set up module params */
        void set_params(int lt, double beta, int nbin, int rebin_pace, int num_bootstrap);

        /* read sequence of imaginary time tau from input file */
        void read_tau_from_file(const std::string &infile_tau_seq);

        /* read sequence of correlation */
        void read_corr_from_file(const std::string &infile_corr_bin);

        /* calculate means and corr-errors of input bin correlations */
        void analyse_corr();

        /* discard correlations with poor data quality,
         * rotate correlations to eigen space of covariance matrix */
        void discard_and_rotate();


    private:

        void allocate_memory();

        /* compute means of correlation */
        void compute_corr_means();

        /* compute corr-errors of correlations */
        void compute_corr_errs();

        /* compute covariance matrix C_ij and diagonalize */
        void compute_cov_matrix();

        /* discard correlations with poor data quality */
        void discard_poor_quality_data();
    };
}


#endif //SAC_READINMODULE_H
