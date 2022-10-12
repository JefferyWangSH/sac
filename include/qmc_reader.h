#ifndef SAC_INITIALIZER_QMC_READER_H
#define SAC_INITIALIZER_QMC_READER_H
#pragma once

/**
  *  This header file defines the `SAC::InitializerQmcReader` class, 
  *  which serves as a interface class for the input of the QMC data.
  *  Features:
  *    1. read the imaginary-time dynamic correlation functions G(t),
  *       in form of the raw data with the bin samples included.
  *    2. generate and diagonalize the covariance matrix, 
  *       further compute its eigenvalue and eigenvector for the SAC usage.
  *       ( The SAC calculation is performed in the eigen space of the covariance matrix. )
  */

#include <string>
#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#include <Eigen/Core>


namespace SAC::Initializer {

    // ---------------------------------------  SAC::Initializer::QmcReader class  --------------------------------------------
    class QmcReader {

        private:

            int m_time_num{};                 // number of imaginary-time points
            int m_bin_num{};                  // number of bins (after rebin)
            int m_bin_num_total{};            // number of bins (before rebin)
            int m_rebin_pace{};               // pace of rebin, 1 for all bins
            int m_bootstrap_num{};            // number of bootstrap samples
            double m_beta{};                  // inverse temperature
            double m_g0{};                    // static correlation at G(t=0), serving as a normalization factor

            int m_cov_mat_dim{};              // dimension of the covariance matrix
            Eigen::MatrixXd m_cov_mat{};      // covariance matrix
            Eigen::MatrixXd m_rotate_mat{};   // orthogonal matrix, which rotates the covariance matrix to its eigen space
            Eigen::VectorXd m_eig_vec{};      // eigenvalue vector of the covariance matrix

            // the processed QMC data useful for the SAC simulation
            Eigen::VectorXd m_tgrids_qmc{}, m_corr_mean_qmc{}, m_corr_stddev_qmc{};

            // intermediate matrices
            Eigen::MatrixXd m_qmc_samples;
            Eigen::MatrixXd m_bootstrap_samples;


        public:

            QmcReader() = default;

            void deallocate_memory();
            
            // set up reader params
            void set_params( int time_size, double beta, int bin_num, int rebin_pace, int bootstrap_num );

            // read the imaginary-time grids from the input file
            void read_tgrids_from_file( const std::string& tgrids_file );

            // read the correlation functions from the input file
            // also rebin the raw data if necessary
            void read_corr_from_file( const std::string& corr_file );

            // compute the mean value and standard deviation of the input correlation functions
            void analyse_corr();
            
            // discard correlation data with poor data quality
            // rotate the correlations to the diagonal space of the covariance matrix
            void filter_and_rotate();

            // interface memeber functions
            int time_num()          const;
            int bin_num()           const;
            int bin_num_total()     const;
            int bootstrap_num()     const;
            int rebin_pace()        const;
            double beta()           const;
            double scaling_factor() const;

            int cov_mat_dim() const;
            const Eigen::MatrixXd& cov_mat()    const;
            const Eigen::MatrixXd& rotate_mat() const;
            const Eigen::VectorXd& eig_vec()    const;

            const Eigen::VectorXd& tgrids_qmc()      const;
            const Eigen::VectorXd& corr_mean_qmc()   const;
            const Eigen::VectorXd& corr_stddev_qmc() const;

        private:
        
            void allocate_memory();

            // compute the mean value of the correlations
            void compute_corr_mean();

            // compute the standard deviation of the correlations
            void compute_corr_stddev();

            // generate the covariance matrix and diagonalize it
            void compute_cov_matrix();

            // discard correlation data with poor data quality
            void discard_poor_quality_data();

    };


} // namespace SAC::Initializer


#endif // SAC_INITIALIZER_QMC_READER_H
