#include "qmc_data_reader.h"
#include "matrix_util.hpp"
#include "random.h"

#include <fstream>
#include <iostream>
#include <boost/format.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/lexical_cast.hpp>


namespace DataReader {

    void QMCDataReader::set_params(int lt, double beta, int nbin, int rebin_pace, int boostrap_num) {
        assert( lt > 0 && beta > 0 );
        assert( nbin > 0 && boostrap_num > 0 );
        assert( rebin_pace > 0 && rebin_pace <= nbin );

        this->lt = lt;
        this->beta = beta;
        this->nbin_total = nbin;
        this->rebin_pace = rebin_pace;
        this->nbin = nbin / rebin_pace;
        this->bootstrap_num = boostrap_num;

        this->allocate_memory();
    }

    void QMCDataReader::allocate_memory() {
        this->tau_qmc.resize(this->lt);
        this->corr_mean_qmc.resize(this->lt);
        this->corr_err_qmc.resize(this->lt);

        this->bin_data_qmc.resize(this->nbin, this->lt);
        this->bootstrap_samples.resize(this->bootstrap_num, this->lt);
    }

    void QMCDataReader::deallocate_memory() {
        // free useless objects which cost large memory
        this->bin_data_qmc.resize(0, 0);
        this->bootstrap_samples.resize(0, 0);

        this->tau_qmc.resize(0);
        this->corr_mean_qmc.resize(0);
        this->corr_err_qmc.resize(0);
    }

    void QMCDataReader::read_tau_from_file(const std::string &tau_file_path) {
        std::ifstream infile(tau_file_path, std::ios::in);
        if (!infile.is_open()) {
            std::cerr << boost::format(" Fail to open file %s, check the input.\n") % tau_file_path << std::endl;
            exit(1);
        }
        // temporary params
        std::string line;
        std::vector<std::string> data;

        // heading message, reading lt and beta
        getline(infile, line);
        boost::split(data, line, boost::is_any_of(" "), boost::token_compress_on);
        data.erase(std::remove(std::begin(data), std::end(data), ""), std::end(data));
        if ( this->lt != boost::lexical_cast<int>(data[0]) || this->beta != boost::lexical_cast<double>(data[1]) ) {
            std::cerr << " Inconsistence between QMC params and input QMC file, check the input. \n" << std::endl;
            exit(1);
        }

        int t = 0;
        while(getline(infile, line)) {
            boost::split(data, line, boost::is_any_of(" "), boost::token_compress_on);
            data.erase(std::remove(std::begin(data), std::end(data), ""), std::end(data));
            this->tau_qmc[t] = boost::lexical_cast<double>(data[0]);
            t++;
        }
        // check the consistence between input data and model settings for second time
        assert( this->lt == t );
        infile.close();
    }

    void QMCDataReader::read_corr_from_file(const std::string &corr_file_path) {
        std::ifstream infile(corr_file_path, std::ios::in);
        if (!infile.is_open()) {
            std::cerr << boost::format(" Fail to open file %s, check the input.\n") % corr_file_path << std::endl;
            exit(1);
        }
        // temporary params
        std::string line;
        std::vector<std::string> data;

        // heading message, reading total number of bins
        getline(infile, line);
        boost::split(data, line, boost::is_any_of(" "), boost::token_compress_on);
        data.erase(std::remove(std::begin(data), std::end(data), ""), std::end(data));
        if ( this->nbin_total != boost::lexical_cast<int>(data[0]) ) {
            std::cerr << " Inconsistence between QMC params and input QMC file, check the input. \n" << std::endl;
            exit(1);
        }
        // clear previous data
        this->bin_data_qmc = Eigen::MatrixXd::Zero(this->nbin, this->lt);

        // rebin and read data
        for (int bin = 0; bin < this->nbin; ++bin) {
            for (int rebin = 0; rebin < this->rebin_pace; ++rebin) {
                getline(infile, line);
                for (int l = 0; l < this->lt; ++l) {
                    getline(infile, line);
                    boost::split(data, line, boost::is_any_of(" "), boost::token_compress_on);
                    data.erase(std::remove(std::begin(data), std::end(data), ""), std::end(data));
                    this->bin_data_qmc(bin, l) += boost::lexical_cast<double>(data[0]);
                }
            }
        }
        this->bin_data_qmc /= this->rebin_pace;
        infile.close();
    }

    void QMCDataReader::compute_corr_means() {
        // clear previous data
        this->corr_mean_qmc.setZero();

        // calculate means of correlations
        for (int l = 0; l < this->lt; ++l) {
            this->corr_mean_qmc[l] = this->bin_data_qmc.col(l).sum() / this->nbin;
        }

        // record static correlation G(0), rescale such that G(0) = 1.0
        this->g0 = this->corr_mean_qmc[0];
    }

    void QMCDataReader::compute_corr_errs() {
        // clear previous data
        this->corr_err_qmc.setZero();

        // firstly generate bootstrap samples
        // `nbin` random selections out of the `nbin` bins
        std::uniform_int_distribution<> rand_bin(0, this->nbin-1);
        for (int i = 0; i < this->bootstrap_num; ++i) {
            for (int bin = 0; bin < this->nbin; bin++) {
                this->bootstrap_samples.row(i) += this->bin_data_qmc.row(rand_bin(Random::Engine));
            }
        }
        this->bootstrap_samples /= this->nbin;

        // compute correlated errors of input correlation functions
        for (int l = 0; l < this->lt; ++l) {
            this->corr_err_qmc[l] = (this->bootstrap_samples.col(l).array() - this->corr_mean_qmc[l]).square().sum();
        }
        this->corr_err_qmc = (this->corr_err_qmc/this->bootstrap_num).array().sqrt().matrix();
    }

    void QMCDataReader::analyse_corr() {
        this->compute_corr_means();
        this->compute_corr_errs();
    }

    void QMCDataReader::discard_poor_quality_data() {
        // discard correlations with poor data quality
        // criteria: relative error less than 0.1
        // helping param containing index of which the correlation has `good` data quality
        std::vector<int> selected_tau;

        // static correlation G(0) excluded
        for (int l = 1; l < this->lt; ++l) {
            if ( abs(this->corr_err_qmc[l] / this->corr_mean_qmc[l]) < 0.1) {
                selected_tau.push_back(l);
            }
        }
        // determine dimension of covariance matrix
        this->cov_mat_dim = selected_tau.size();

        // allocate memory
        this->cov_mat.resize(this->cov_mat_dim, this->cov_mat_dim);
        this->cov_eig.resize(this->cov_mat_dim);
        this->rotate_mat.resize(this->cov_mat_dim, this->cov_mat_dim);

        // reshape data and bootstrap samples
        // prerequisite: Eigen version > 3.3.9
        const Eigen::VectorXd& tmp_tau = this->tau_qmc(selected_tau);
        const Eigen::VectorXd& tmp_corr_mean = this->corr_mean_qmc(selected_tau);
        const Eigen::VectorXd& tmp_corr_err = this->corr_err_qmc(selected_tau);
        const Eigen::MatrixXd& tmp_bootstrap_samples = this->bootstrap_samples(Eigen::all, selected_tau);
        this->tau_qmc = tmp_tau;
        this->corr_mean_qmc = tmp_corr_mean;
        this->corr_err_qmc = tmp_corr_err;
        this->bootstrap_samples = tmp_bootstrap_samples;
    }

    void QMCDataReader::compute_cov_matrix() {
        // clear previous data
        this->cov_mat = Eigen::MatrixXd::Zero(this->cov_mat_dim, this->cov_mat_dim);

        // compute covariance matrix
        // means and errors of correlations with `poor` quality are discarded ahead of time
        for (int i = 0; i < this->cov_mat_dim; ++i) {
            for (int j = 0; j < this->cov_mat_dim; ++j) {
                // fixed G(0), not included in the data set defining chi square
                this->cov_mat(i, j) = ( (this->bootstrap_samples.col(i).array() - this->corr_mean_qmc[i])
                                    * (this->bootstrap_samples.col(j).array() - this->corr_mean_qmc[j]) ).sum();
            }
        }
    }

    void QMCDataReader::discard_and_rotate() {
        // discard correlations with poor data quality
        this->discard_poor_quality_data();

        // rescaled correlations by G(0)
        this->corr_mean_qmc /= this->g0;
        this->corr_err_qmc /= this->g0;
        this->bootstrap_samples /= this->g0;

        // compute covariance matrix
        this->compute_cov_matrix();

        // diagonalize covariance matrix
        // for a real symmetric matrix, orthogonal transformation T satisfies
        //   T * C * T^dagger -> Eigen space
        // where T is rotation matrix, which is orthogonal.
        MatrixUtil::mkl_lapack_dsyev(this->cov_mat_dim, this->cov_mat, this->cov_eig, this->rotate_mat);

        // alternative method with lower accuracy, using SVD
        // Eigen::MatrixXd u(cov_mat_dim, cov_mat_dim), v(cov_mat_dim, cov_mat_dim);
        // mkl_lapack_dgesvd(cov_mat_dim, cov_mat_dim, cov_mat, u, cov_eig, v);
        // rotate_mat = u.transpose();
    }

} // namespace DataRead
