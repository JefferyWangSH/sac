#include <cassert>
#include "Measure.h"
#include "MonteCarloSAC.h"

Measure::Measure(int nbin, int nconfig) {
    this->nbin = nbin;
    this->nconfig = nconfig;
    prepare();
}

void Measure::resize(int _nbin, int _nconfig) {
    this->nbin = _nbin;
    this->nconfig = _nconfig;
    prepare();
}

void Measure::prepare() {
    bin_Hamilton.clear();
    bin_n_x.clear();
    bin_Hamilton.reserve(nbin);
    bin_n_x.reserve(nbin);
    for (int bin = 0; bin < nbin; ++bin) {
        bin_Hamilton.emplace_back(0.0);
        bin_n_x.emplace_back(nconfig, 0.0);
    }
    bin_Hamilton.shrink_to_fit();
    bin_n_x.shrink_to_fit();

    mean_Hamilton = 0.0;
    err_Hamilton = 0.0;
    tmp_Hamilton = 0.0;

    mean_n_x.clear();
    err_n_x.clear();
    tmp_n_x.clear();
    mean_n_x.reserve(nconfig);
    err_n_x.reserve(nconfig);
    tmp_n_x.reserve(nconfig);

    for (int i = 0; i < nconfig; ++i) {
        mean_n_x.emplace_back(0.0);
        err_n_x.emplace_back(0.0);
        tmp_n_x.emplace_back(0.0);
    }
    mean_n_x.shrink_to_fit();
    err_n_x.shrink_to_fit();
    tmp_n_x.shrink_to_fit();
}

void Measure::clear_tmp_stats() {
    assert( tmp_n_x.size() == nconfig );

    nn = 0;
    tmp_Hamilton = 0.0;
    for (int i = 0; i < nconfig; ++i) {
        tmp_n_x[i] = 0.0;
    }
}

void Measure::measure(MonteCarloSAC& sac) {
    assert( sac.nconfig == nconfig );
    assert( sac.nbin == nbin );

    double tmp_h = 0.0;
    sac.cal_Config_Hamiltonian(tmp_h);
    tmp_Hamilton += tmp_h;

    for (int i = 0; i < nconfig; ++i) {
        tmp_n_x[i] += sac.n_list[i];
    }
    nn++;
}

void Measure::normalize_stats() {
    tmp_Hamilton /= nn;
    for (int i = 0; i < nconfig; ++i) {
        tmp_n_x[i] /= nn;
    }
}

void Measure::write_data_to_bin(const int& bin) {
    assert( bin >= 0 && bin < nbin );

    bin_Hamilton[bin] = tmp_Hamilton;
    for (int i = 0; i < nconfig; ++i) {
        bin_n_x[bin][i] = tmp_n_x[i];
    }
}

void Measure::analyse_stats() {

    // clear previous data
    mean_Hamilton = 0.0;
    err_Hamilton = 0.0;
    for (int i = 0; i < nconfig; ++i) {
        mean_n_x[i] = 0.0;
        err_n_x[i] = 0.0;
    }

    // average over bins
    for (int bin = 0; bin < nbin; ++bin) {
        mean_Hamilton += bin_Hamilton[bin];
        err_Hamilton += bin_Hamilton[bin] * bin_Hamilton[bin];

        for (int i = 0; i < nconfig; ++i) {
            mean_n_x[i] += bin_n_x[bin][i];
            err_n_x[i] += bin_n_x[bin][i] * bin_n_x[bin][i];
        }
    }

    // normalize data
    mean_Hamilton /= nbin;
    err_Hamilton /= nbin;
    for (int i = 0; i < nconfig; ++i) {
        mean_n_x[i] /= nbin;
        err_n_x[i] /= nbin;
    }

    // calculate statistical error
    err_Hamilton = pow(err_Hamilton - pow(mean_Hamilton, 2), 0.5) / pow(nbin - 1, 0.5);
    for (int i = 0; i < nalpha; ++i) {
        err_n_x[i] = pow(err_n_x[i] - pow(mean_n_x[i], 2), 0.5) / pow(nbin - 1, 0.5);
    }
}

