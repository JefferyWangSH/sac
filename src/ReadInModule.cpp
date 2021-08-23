#include "ReadInModule.h"
#include "DiagonalTools.hpp"

#include <fstream>
#include <iostream>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/lexical_cast.hpp>


void QMCData::ReadInModule::set_params(int _lt, double _beta, int _nbin, int _rebin_pace, int _num_boostrap) {
    assert( _lt > 0 && _beta > 0 );
    assert( _nbin > 0 && _rebin_pace > 0 && _num_boostrap > 0 );
    assert( _rebin_pace <= _nbin );

    this->lt = _lt;
    this->beta = _beta;
    this->nbin_total = _nbin;
    this->rebin_pace = _rebin_pace;
    this->nbin = _nbin / _rebin_pace;
    this->num_bootstrap = _num_boostrap;

    this->allocate_memory();
}

void QMCData::ReadInModule::allocate_memory() {
    this->tau.resize(lt);
    this->corr_mean.resize(lt);
    this->corr_err.resize(lt);

    this->corr_tau_bin.resize(nbin, lt);
    this->sample_bootstrap.resize(num_bootstrap, lt);
}

void QMCData::ReadInModule::deallocate_memory() {
    // free useless objects which cost large memory
    this->corr_tau_bin.resize(0, 0);
    this->sample_bootstrap.resize(0, 0);

    this->tau.resize(0);
    this->corr_mean.resize(0);
    this->corr_err.resize(0);
}

void QMCData::ReadInModule::read_tau_from_file(const std::string &infile_tau_seq) {
    std::ifstream infile;
    infile.open(infile_tau_seq, std::ios::in);

    if (!infile.is_open()) {
        std::cerr << "=====================================================================" << std::endl
                  << " fail to open file " + infile_tau_seq + ", please check the input." << std::endl
                  << "=====================================================================" << std::endl;
        exit(1);
    }
    // temporary params
    std::string line;
    std::vector<std::string> data;

    // heading message, reading lt and beta
    getline(infile, line);
    boost::split(data, line, boost::is_any_of(" "), boost::token_compress_on);
    data.erase(std::remove(std::begin(data), std::end(data), ""), std::end(data));
    assert( this->lt == boost::lexical_cast<int>(data[0]) );
    assert( this->beta == boost::lexical_cast<double>(data[1]) );

    int t = 0;
    while(getline(infile, line)) {
        boost::split(data, line, boost::is_any_of(" "), boost::token_compress_on);
        data.erase(std::remove(std::begin(data), std::end(data), ""), std::end(data));
        this->tau[t] = boost::lexical_cast<double>(data[0]);
        t++;
    }

    // check the consistence between data and model settings
    assert( this->lt == t );
    infile.close();
}

void QMCData::ReadInModule::read_corr_from_file(const std::string &infile_g_bin) {
    std::ifstream infile;
    infile.open(infile_g_bin, std::ios::in);

    if (!infile.is_open()) {
        std::cerr << "=====================================================================" << std::endl
                  << " fail to open file " + infile_g_bin + ", please check the input." << std::endl
                  << "=====================================================================" << std::endl;
        exit(1);
    }
    // temporary params
    std::string line;
    std::vector<std::string> data;

    // heading message, reading total number of bins
    getline(infile, line);
    boost::split(data, line, boost::is_any_of(" "), boost::token_compress_on);
    data.erase(std::remove(std::begin(data), std::end(data), ""), std::end(data));
    this->nbin_total = boost::lexical_cast<int>(data[0]);

    // clear previous data
    this->corr_tau_bin = Eigen::MatrixXd::Zero(this->nbin, this->lt);

    for (int bin = 0; bin < this->nbin; ++bin) {
        for (int rebin = 0; rebin < this->rebin_pace; ++rebin) {
            getline(infile, line);
            for (int l = 0; l < this->lt; ++l) {
                getline(infile, line);
                boost::split(data, line, boost::is_any_of(" "), boost::token_compress_on);
                data.erase(std::remove(std::begin(data), std::end(data), ""), std::end(data));
                this->corr_tau_bin(bin, l) += boost::lexical_cast<double>(data[0]);
            }
        }
    }
    this->corr_tau_bin /= this->rebin_pace;
    infile.close();
}

void QMCData::ReadInModule::compute_corr_means() {
    // clear previous data
    this->corr_mean = Eigen::VectorXd::Zero(this->lt);

    // calculate means of correlations
    for (int l = 0; l < lt; ++l) {
        this->corr_mean[l] = this->corr_tau_bin.col(l).sum() / this->nbin;
    }

    // record static correlation G(0), rescale such that G(0) = 1.0
    this->g0 = this->corr_mean[0];
}

void QMCData::ReadInModule::compute_corr_errs() {
    // clear previous data
    this->corr_err = Eigen::VectorXd::Zero(this->lt);

    // first generate bootstrap samples
    // `nbin` random selections out of the `nbin` bins
    std::uniform_int_distribution<> rand_bin(0, this->nbin - 1);
    for (int i = 0; i < this->num_bootstrap; ++i) {
        for (int bin = 0; bin < this->nbin; bin++) {
            this->sample_bootstrap.row(i) += this->corr_tau_bin.row(rand_bin(rand_engine_readin));
        }
    }
    this->sample_bootstrap /= this->nbin;

    // compute corr-errors of correlations
    for (int l = 0; l < lt; ++l) {
        this->corr_err[l] = (this->sample_bootstrap.col(l).array() - this->corr_mean[l]).square().sum();
        this->corr_err[l] = sqrt(this->corr_err[l] / this->num_bootstrap);
    }
}

void QMCData::ReadInModule::analyse_corr() {
    this->compute_corr_means();
    this->compute_corr_errs();
}

void QMCData::ReadInModule::discard_poor_quality_data() {
    // discard correlations with poor data quality
    // criteria: relative error less than 0.1
    // helping param which contain tau index where the correlation has `good` data quality
    std::vector<int> selected_tau;

    // static correlation G(0) excluded
    for (int l = 1; l < lt; ++l) {
        if ( abs(this->corr_err[l] / this->corr_mean[l]) < 0.1) {
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
    // prerequisite: Eigen version > 3.9.9
    this->tau = this->tau(selected_tau);
    this->corr_mean = this->corr_mean(selected_tau);
    this->corr_err = this->corr_err(selected_tau);

    const Eigen::VectorXi rows = Eigen::VectorXi::LinSpaced(this->num_bootstrap, 0, this->num_bootstrap);
    const Eigen::MatrixXd sample_tmp = this->sample_bootstrap(rows, selected_tau);
    this->sample_bootstrap = sample_tmp;
}

void QMCData::ReadInModule::compute_cov_matrix() {
    // clear previous data
    this->cov_mat = Eigen::MatrixXd::Zero(this->cov_mat_dim, this->cov_mat_dim);

    // compute covariance matrix
    // means and errors of correlations with `poor` quality are discarded ahead of time
    for (int i = 0; i < this->cov_mat_dim; ++i) {
        for (int j = 0; j < this->cov_mat_dim; ++j) {
            // fixed G(0), not included in the data set defining chi square
            this->cov_mat(i, j) = ( (this->sample_bootstrap.col(i).array() - this->corr_mean[i])
                                * (this->sample_bootstrap.col(j).array() - this->corr_mean[j]) ).sum();
        }
    }
}

void QMCData::ReadInModule::discard_and_rotate() {
    // discard correlations with poor data quality
    this->discard_poor_quality_data();

    // rescaled correlations by G(0)
    this->corr_mean /= this->g0;
    this->corr_err /= this->g0;
    this->sample_bootstrap /= this->g0;

    // compute covariance matrix
    this->compute_cov_matrix();

    // diagonalize covariance matrix
    // for a real symmetric matrix, orthogonal transformation T satisfies
    //   T * C * T^dagger -> Eigen space
    // where T is rotation matrix, which is orthogonal.
    mkl_lapack_dsyev(this->cov_mat_dim, this->cov_mat, this->cov_eig, this->rotate_mat);

    // alternative method with lower accuracy, using SVD
    // Eigen::MatrixXd u(cov_mat_dim, cov_mat_dim), v(cov_mat_dim, cov_mat_dim);
    // mkl_lapack_dgesvd(cov_mat_dim, cov_mat_dim, cov_mat, u, cov_eig, v);
    // rotate_mat = u.transpose();
}

