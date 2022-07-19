#include "qmc_reader.h"
#include "random.h"
#include "utils/linear_algebra.hpp"

#include <fstream>
#include <iostream>
#include <boost/format.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/lexical_cast.hpp>


namespace SAC {

    // interface memeber functions
    int QmcReader::time_num() const { return this->m_time_num; }
    int QmcReader::bin_num() const { return this->m_bin_num; }
    int QmcReader::bin_num_total() const { return this->m_bin_num_total; }
    int QmcReader::bootstrap_num() const { return this->m_bootstrap_num; }
    int QmcReader::rebin_pace() const { return this->m_rebin_pace; }
    double QmcReader::beta() const { return this->m_beta; }
    double QmcReader::scaling_factor() const { return this->m_g0; }

    int QmcReader::cov_mat_dim() const { return this->m_cov_mat_dim; }
    const Eigen::MatrixXd& QmcReader::cov_mat() const { return this->m_cov_mat; }
    const Eigen::MatrixXd& QmcReader::rotate_mat() const { return this->m_rotate_mat; }
    const Eigen::VectorXd& QmcReader::eig_vec() const { return this->m_eig_vec; }

    const Eigen::VectorXd& QmcReader::tgrids_qmc() const { return this->m_tgrids_qmc; }
    const Eigen::VectorXd& QmcReader::corr_mean_qmc() const { return this->m_corr_mean_qmc; }
    const Eigen::VectorXd& QmcReader::corr_err_qmc() const { return this->m_corr_err_qmc; }


    void QmcReader::set_params( int time_size, double beta, int bin_num, int rebin_pace, int bootstrap_num ) 
    {
        assert( time_size > 0 && beta > 0 );
        assert( bin_num > 0 && bootstrap_num > 0 );
        assert( rebin_pace > 0 && rebin_pace <= bin_num );

        this->m_time_num = time_size;
        this->m_beta = beta;
        this->m_bin_num_total = bin_num;
        this->m_rebin_pace = rebin_pace;
        this->m_bin_num = bin_num / rebin_pace;
        this->m_bootstrap_num = bootstrap_num;

        this->allocate_memory();
    }


    void QmcReader::allocate_memory() 
    {
        this->m_tgrids_qmc.resize(this->m_time_num);
        this->m_corr_mean_qmc.resize(this->m_time_num);
        this->m_corr_err_qmc.resize(this->m_time_num);

        this->m_bin_data_qmc.resize(this->m_bin_num, this->m_time_num);
        this->m_bootstrap_samples.resize(this->m_bootstrap_num, this->m_time_num);
    }


    void QmcReader::deallocate_memory() 
    {
        // free useless objects which would otherwise take lots of memory
        this->m_bin_data_qmc.resize(0, 0);
        this->m_bootstrap_samples.resize(0, 0);

        // this->m_tgrids_qmc.resize(0);
        // this->m_corr_mean_qmc.resize(0);
        // this->m_corr_err_qmc.resize(0);
    }


    void QmcReader::read_tgrids_from_file( const std::string& tgrids_file ) 
    {
        std::ifstream infile(tgrids_file, std::ios::in);
        if ( !infile.is_open() ) {
            std::cerr << "SAC::Initializer::QmcReader::read_tgrids_from_file(): "
                      << "fail to open file \'" << tgrids_file << "\'." << std::endl;
            exit(1);
        }

        // intermediate variables
        std::string line;
        std::vector<std::string> data;

        // read the heading messages
        // number of the time points and the inverse temperature
        getline(infile, line);
        boost::split(data, line, boost::is_any_of(" "), boost::token_compress_on);
        data.erase(std::remove(std::begin(data), std::end(data), ""), std::end(data));
        if ( this->m_time_num != boost::lexical_cast<int>(data[0]) ) {
            std::cerr << "SAC::Initializer::QmcReader::read_tgrids_from_file(): "
                      << "inconsistence of the number of imaginary-time points "
                      << "between SAC settings and the input tgrids file." << std::endl;
            exit(1);
        }
        if ( this->m_beta != boost::lexical_cast<double>(data[1]) ) {
            std::cerr << "SAC::Initializer::QmcReader::read_tgrids_from_file(): "
                      << "inconsistence of the inverse temperature "
                      << "between SAC settings and the input tgrids file." << std::endl;
            exit(1);
        }

        for ( auto t = 0; t < this->m_time_num; ++t ) {
            getline(infile, line);
            boost::split(data, line, boost::is_any_of(" "), boost::token_compress_on);
            data.erase(std::remove(std::begin(data), std::end(data), ""), std::end(data));
            this->m_tgrids_qmc[t] = boost::lexical_cast<double>(data[1]);
        }
        infile.close();
    }


    void QmcReader::read_corr_from_file( const std::string& corr_file ) 
    {
        std::ifstream infile(corr_file, std::ios::in);
        if ( !infile.is_open() ) {
            std::cerr << "SAC::Initializer::QmcReader::read_corr_from_file(): "
                      << "fail to open file \'" << corr_file << "\'." << std::endl;
            exit(1);
        }
        
        // intermediate variables
        std::string line;
        std::vector<std::string> data;

        // read the heading message, total number of bins
        getline(infile, line);
        boost::split(data, line, boost::is_any_of(" "), boost::token_compress_on);
        data.erase(std::remove(std::begin(data), std::end(data), ""), std::end(data));
        if ( this->m_bin_num_total != boost::lexical_cast<int>(data[0]) ) {
            std::cerr << "SAC::Initializer::QmcReader::read_corr_from_file(): "
                      << "inconsistence of the total number of bins "
                      << "between SAC settings and the input corr file." << std::endl;
            exit(1);
        }
        if ( this->m_time_num != boost::lexical_cast<int>(data[1]) ) {
            std::cerr << "SAC::Initializer::QmcReader::read_corr_from_file(): "
                      << "inconsistence of the number of imaginary-time points "
                      << "between SAC settings and the input corr file." << std::endl;
            exit(1);
        }
        // clear previous data
        this->m_bin_data_qmc = Eigen::MatrixXd::Zero(this->m_bin_num, this->m_time_num);

        // read and rebin the input data
        for ( auto bin = 0; bin < this->m_bin_num; ++bin ) {
            for ( auto rebin = 0; rebin < this->m_rebin_pace; ++rebin ) {
                for ( auto t = 0; t < this->m_time_num; ++t ) {
                    getline(infile, line);
                    boost::split(data, line, boost::is_any_of(" "), boost::token_compress_on);
                    data.erase(std::remove(std::begin(data), std::end(data), ""), std::end(data));
                    this->m_bin_data_qmc(bin, t) += boost::lexical_cast<double>(data[2]);
                }
            }
        }
        this->m_bin_data_qmc /= this->m_rebin_pace;
        infile.close();
    }


    void QmcReader::compute_corr_means() 
    {
        // clear previous data
        this->m_corr_mean_qmc.setZero();

        // calculate the mean values of correlations
        for ( auto t = 0; t < this->m_time_num; ++t ) {
            this->m_corr_mean_qmc[t] = this->m_bin_data_qmc.col(t).sum() / this->m_bin_num;
        }

        // record the static correlation G(t=0), which serves as the scaling factor in SAC
        // we will rescale the input correaltion data such that G(t=0) = 1.0,
        // and this is exactly the normalization condition for the recovered spectral function.
        this->m_g0 = this->m_corr_mean_qmc[0];
    }


    void QmcReader::compute_corr_errs() 
    {
        // clear previous data
        this->m_corr_err_qmc.setZero();
        this->m_bootstrap_samples.setZero();

        // first generate the bootstrap samples to compute the covariance matrix
        // `nbin` random selections out of the `nbin` number of bins
        std::uniform_int_distribution<> rand_bin(0, this->m_bin_num-1);
        for ( auto i = 0; i < this->m_bootstrap_num; ++i ) {
            for ( auto bin = 0; bin < this->m_bin_num; bin++ ) {
                this->m_bootstrap_samples.row(i) += this->m_bin_data_qmc.row( rand_bin(Utils::Random::Engine) );
            }
        }
        this->m_bootstrap_samples /= this->m_bin_num;

        // compute correlated errors of the input correlation functions
        for ( auto t = 0; t < this->m_time_num; ++t ) {
            this->m_corr_err_qmc[t] = ( this->m_bootstrap_samples.col(t).array() - this->m_corr_mean_qmc[t] ).square().sum();
        }
        this->m_corr_err_qmc = ( this->m_corr_err_qmc / this->m_bootstrap_num ).array().sqrt().matrix();
    }


    void QmcReader::analyse_corr() 
    {
        this->compute_corr_means();
        this->compute_corr_errs();
    }


    void QmcReader::discard_poor_quality_data() 
    {
        // discard correlations with poor data quality
        // criteria: the data with good data quality should have a relative error smaller than 0.1.
        // intermediate variable which contains indices of the correlations with `good` data quality
        std::vector<int> good_tgrids;

        // exclude the static correlation G(t=0)
        // which serves as the normalization condition for the recovered spectrum
        for ( auto t = 1; t < this->m_time_num; ++t ) {
            if ( std::abs(this->m_corr_err_qmc[t] / this->m_corr_mean_qmc[t]) < 0.1) {
                good_tgrids.push_back(t);
            }
        }
        // determine the dimension of the covariance matrix
        this->m_cov_mat_dim = good_tgrids.size();

        // allocate memory
        this->m_cov_mat.resize(this->m_cov_mat_dim, this->m_cov_mat_dim);
        this->m_eig_vec.resize(this->m_cov_mat_dim);
        this->m_rotate_mat.resize(this->m_cov_mat_dim, this->m_cov_mat_dim);

        // reshape data and bootstrap samples
        // this requires Eigen version > 3.3.9
        const Eigen::VectorXd& tmp_tgrids = this->m_tgrids_qmc(good_tgrids);
        const Eigen::VectorXd& tmp_corr_mean = this->m_corr_mean_qmc(good_tgrids);
        const Eigen::VectorXd& tmp_corr_err = this->m_corr_err_qmc(good_tgrids);
        const Eigen::MatrixXd& tmp_bootstrap_samples = this->m_bootstrap_samples(Eigen::all, good_tgrids);
        this->m_tgrids_qmc = tmp_tgrids;
        this->m_corr_mean_qmc = tmp_corr_mean;
        this->m_corr_err_qmc = tmp_corr_err;
        this->m_bootstrap_samples = tmp_bootstrap_samples;
    }


    void QmcReader::compute_cov_matrix() 
    {
        // clear previous data
        this->m_cov_mat = Eigen::MatrixXd::Zero(this->m_cov_mat_dim, this->m_cov_mat_dim);

        // compute the covariance matrix
        // correlations data with `poor` quality should be discarded in advance
        for ( auto i = 0; i < this->m_cov_mat_dim; ++i ) {
            for ( auto j = 0; j < this->m_cov_mat_dim; ++j ) {
                this->m_cov_mat(i, j) = ( (this->m_bootstrap_samples.col(i).array() - this->m_corr_mean_qmc[i])
                                        * (this->m_bootstrap_samples.col(j).array() - this->m_corr_mean_qmc[j]) ).sum();
            }
        }
    }


    void QmcReader::filter_and_rotate() 
    {
        // discard correlations with poor data quality
        this->discard_poor_quality_data();

        // rescale the correlation data by G(t=0)
        this->m_corr_mean_qmc /= this->m_g0;
        this->m_corr_err_qmc /= this->m_g0;
        this->m_bootstrap_samples /= this->m_g0;

        // compute the covariance matrix
        this->compute_cov_matrix();

        // diagonalize the covariance matrix
        // for a real symmetric matrix, this exists an orthogonal transformation T satisfying
        //   T * C * T^dagger -> diagonal space
        // where T is the orthogonal rotation matrix.
        Utils::LinearAlgebra::mkl_lapack_dsyev(this->m_cov_mat_dim, this->m_cov_mat, this->m_eig_vec, this->m_rotate_mat);

        // alternative method with lower accuracy, using SVD
        // Eigen::MatrixXd u(this->m_cov_mat_dim, this->m_cov_mat_dim), v(this->m_cov_mat_dim, this->m_cov_mat_dim);
        // Utils::LinearAlgebra::mkl_lapack_dgesvd(this->m_cov_mat_dim, this->m_cov_mat_dim, this->m_cov_mat, u, this->m_eig_vec, v);
        // this->m_rotate_mat = u.transpose();
    }

} // namespace SAC
