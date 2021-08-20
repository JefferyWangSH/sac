#include "ReadInModule.h"
#include "DiagonalTools.hpp"

#include <fstream>
#include <iostream>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/lexical_cast.hpp>


void ReadInModule::set_params(int _lt, double _beta, int _nbin, int _rebin_pace, int _num_boostrap) {
    assert( _lt > 0 && _beta > 0 );
    assert( _nbin > 0 && _rebin_pace > 0 && _num_boostrap > 0 );
    assert( _rebin_pace <= _nbin );

    this->lt = _lt;
    this->beta = _beta;
    this->nbin_total = _nbin;
    this->rebin_pace = _rebin_pace;
    this->nbin = _nbin / _rebin_pace;
    this->num_bootstrap = _num_boostrap;
    this->cov_mat_dim = _lt - 1;                // FIXME: discard data with poor quality

    allocate_memory();
}

void ReadInModule::allocate_memory() {

    this->tau_seq_raw.resize(lt);
    this->corr_mean_seq_raw.resize(lt);
    this->corr_err_seq_raw.resize(lt);

    this->tau_seq.resize(cov_mat_dim);
    this->corr_mean_seq.resize(cov_mat_dim);
    this->corr_err_seq.resize(cov_mat_dim);

    this->cov_mat.resize(cov_mat_dim, cov_mat_dim);
    this->rotate_mat.resize(cov_mat_dim, cov_mat_dim);
    this->cov_eig.resize(cov_mat_dim);

    this->corr_tau_bin.resize(nbin, lt);
    this->sample_bootstrap.resize(num_bootstrap, lt);
}

void ReadInModule::deallocate_memory() {
    // free useless objects which cost large memory
    this->corr_tau_bin.resize(0, 0);
    this->sample_bootstrap.resize(0, 0);
}

void ReadInModule::read_tau_from_file(const std::string &infile_tau_seq) {
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
        tau_seq_raw[t] = boost::lexical_cast<double>(data[0]);
        t++;
    }

    // check the consistence between data and model settings
    assert( this->lt == t );
    infile.close();
}

void ReadInModule::read_corr_from_file(const std::string &infile_g_bin) {
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
    corr_tau_bin = Eigen::MatrixXd::Zero(nbin, lt);

    for (int bin = 0; bin < nbin; ++bin) {
        for (int rebin = 0; rebin < rebin_pace; ++rebin) {
            getline(infile, line);
            for (int l = 0; l < lt; ++l) {
                getline(infile, line);
                boost::split(data, line, boost::is_any_of(" "), boost::token_compress_on);
                data.erase(std::remove(std::begin(data), std::end(data), ""), std::end(data));
                corr_tau_bin(bin, l) += boost::lexical_cast<double>(data[0]);
            }
        }
    }
    corr_tau_bin /= rebin_pace;
    infile.close();
}

void ReadInModule::compute_corr_means() {
    // clear previous data
    corr_mean_seq_raw = Eigen::VectorXd::Zero(lt);
    corr_err_seq_raw = Eigen::VectorXd::Zero(lt);

    // calculate means of correlations
    for (int l = 0; l < lt; ++l) {
        corr_mean_seq_raw[l] = corr_tau_bin.col(l).sum() / nbin;
    }

    // generate bootstrap samples
    for (int i = 0; i < num_bootstrap; ++i) {
        for (int l = 0; l < lt; ++l) {
            for(int bin = 0; bin < nbin; bin++) {
                std::uniform_int_distribution<> rand_bin(0, nbin - 1);
                sample_bootstrap(i, l) += corr_tau_bin(rand_bin(rand_engine_readin), l);
            }
        }
    }
    sample_bootstrap /= nbin;

    // compute corr-errors of correlations
    for (int l = 0; l < lt; ++l) {
        corr_err_seq_raw[l] = ( (sample_bootstrap.col(l).array() - corr_mean_seq_raw[l])
                              * (sample_bootstrap.col(l).array() - corr_mean_seq_raw[l]) ).sum();
        corr_err_seq_raw[l] = sqrt(corr_err_seq_raw[l] / num_bootstrap);
    }
}

void ReadInModule::compute_cov_matrix() {
    // record static correlation G(0), rescale such that G(0) = 1.0
    this->g0 = corr_mean_seq_raw[0];

    // FIXME: discard data with poor quality
    corr_mean_seq = corr_mean_seq_raw.tail(corr_mean_seq_raw.size() - 1) / g0;
    corr_err_seq = corr_err_seq_raw.tail(corr_err_seq_raw.size() - 1) / g0;
    sample_bootstrap /= g0;

    // clear previous data
    cov_mat = Eigen::MatrixXd::Zero(cov_mat_dim, cov_mat_dim);

    for (int i = 0; i < cov_mat_dim; ++i) {
        for (int j = 0; j < cov_mat_dim; ++j) {
            // fixed G(0), not included in the data set defining chi square
            cov_mat(i, j) = ( (sample_bootstrap.col(i+1).array() - corr_mean_seq[i])
                            * (sample_bootstrap.col(j+1).array() - corr_mean_seq[j]) ).sum();
        }
    }

    // diagonalize covariance matrix
    // for a real symmetric matrix, orthogonal transformation T satisfies
    //   T * C * T^dagger -> Eigen space
    // where T is rotation matrix, which is orthogonal.
    Eigen::MatrixXd u(cov_mat_dim, cov_mat_dim), v(cov_mat_dim, cov_mat_dim);
    mkl_lapack_dgesvd(cov_mat_dim, cov_mat_dim, cov_mat, u, cov_eig, v);
    rotate_mat = u.transpose();
}
