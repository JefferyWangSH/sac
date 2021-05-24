#include <cassert>
#include "Measure.h"
#include "MonteCarloSAC.h"

Measure::Measure(int nbin, int nalpha, int nconfig) {
    this->nbin = nbin;
    this->nalpha = nalpha;
    this->nconfig = nconfig;

    prepare();
}

void Measure::resize(int _nbin, int _nalpha, int _nconfig) {
    this->nbin = _nbin;
    this->nalpha = _nalpha;
    this->nconfig = _nconfig;

    prepare();
}

void Measure::prepare() {
    bin_H_alpha.clear();
    bin_n_x_alpha.clear();
    bin_H_alpha.reserve(nbin);
    bin_n_x_alpha.reserve(nbin);
    for (int bin = 0; bin < nbin; ++bin) {
        bin_H_alpha.emplace_back(nalpha, 0.0);
        bin_n_x_alpha.emplace_back(nalpha);

        for (int n = 0; n < nalpha; ++n) {
            bin_n_x_alpha.back()[n].reserve(nconfig);
            for (int i = 0; i < nconfig; ++i) {
                bin_n_x_alpha.back()[n].emplace_back(0.0);
            }
        }
    }
    bin_H_alpha.shrink_to_fit();
    bin_n_x_alpha.shrink_to_fit();

    H_alpha.clear();
    err_H_alpha.clear();
    tmp_H_alpha.clear();
    n_x_alpha.clear();
    err_n_x_alpha.clear();
    tmp_n_x_alpha.clear();

    H_alpha.reserve(nalpha);
    err_H_alpha.reserve(nalpha);
    tmp_H_alpha.reserve(nalpha);
    n_x_alpha.reserve(nalpha);
    err_n_x_alpha.reserve(nalpha);
    tmp_n_x_alpha.reserve(nalpha);

    for (int n = 0; n < nalpha; ++n) {
        H_alpha.emplace_back(0.0);
        err_H_alpha.emplace_back(0.0);
        tmp_H_alpha.emplace_back(0.0);
        n_x_alpha.emplace_back(nconfig, 0.0);
        err_n_x_alpha.emplace_back(nconfig, 0.0);
        tmp_n_x_alpha.emplace_back(nconfig, 0.0);
    }
    H_alpha.shrink_to_fit();
    err_H_alpha.shrink_to_fit();
    tmp_H_alpha.shrink_to_fit();
    n_x_alpha.shrink_to_fit();
    err_n_x_alpha.shrink_to_fit();
    tmp_n_x_alpha.shrink_to_fit();
}

void Measure::clear_tmp_stats() {
    assert( tmp_H_alpha.size() == nalpha );
    assert( tmp_n_x_alpha.size() == nalpha );

    for (int n = 0; n < nalpha; ++n) {
        assert( tmp_n_x_alpha[n].size() == nconfig );
        tmp_H_alpha[n] = 0.0;
        for ( int i = 0; i < nconfig; ++i) {
            tmp_n_x_alpha[n][i] = 0.0;
        }
    }
    nn = 0;
}

void Measure::measure(MonteCarloSAC &sac) {
    assert( sac.nalpha == nalpha );
    assert( sac.nconfig == nconfig );
    assert( sac.nbin == nbin );

    for (int n = 0; n < nalpha; ++n) {
        double aux_H_alpha = 0.0;
        sac.sac_list[n].cal_Config_Hamiltonian(aux_H_alpha);
        tmp_H_alpha[n] += aux_H_alpha;

        for (int i = 0; i < nconfig; ++i) {
            tmp_n_x_alpha[n][i] += sac.sac_list[n].n_list[i];
        }
    }
    nn++;
}

void Measure::normalize_stats() {
    for (int n = 0; n < nalpha; ++n) {
        tmp_H_alpha[n] /= nn;
        for (int i = 0; i < nconfig; ++i) {
            tmp_n_x_alpha[n][i] /= nn;
        }
    }
}

void Measure::write_data_to_bin(const int &bin) {
    for (int n = 0; n < nalpha; ++n) {
        bin_H_alpha[bin][n] = tmp_H_alpha[n];
        for (int i = 0; i < nconfig; ++i) {
            bin_n_x_alpha[bin][n][i] = tmp_n_x_alpha[n][i];
        }
    }
}

void Measure::analyse_stats() {

    // clear previous data
    for (int n = 0; n < nalpha; ++n) {
        H_alpha[n] = 0.0;
        err_H_alpha[n] = 0.0;
        for (int i = 0; i < nconfig; ++i) {
            n_x_alpha[n][i] = 0.0;
            err_n_x_alpha[n][i] = 0.0;
        }
    }

    // average over bins
    for (int bin = 0; bin < nbin; ++bin) {
        for (int n = 0; n < nalpha; ++n) {
            H_alpha[n] += bin_H_alpha[bin][n];
            err_H_alpha[n] += bin_H_alpha[bin][n] * bin_H_alpha[bin][n];
            for (int i = 0; i < nconfig; ++i) {
                n_x_alpha[n][i] += bin_n_x_alpha[bin][n][i];
                err_n_x_alpha[n][i] += bin_n_x_alpha[bin][n][i] * bin_n_x_alpha[bin][n][i];
            }
        }
    }

    // normalize data
    for (int n = 0; n < nalpha; ++n) {
        H_alpha[n] /= nbin;
        err_H_alpha[n] /= nbin;
        for (int i = 0; i < nconfig; ++i) {
            n_x_alpha[n][i] /= nbin;
            err_n_x_alpha[n][i] /= nbin;
        }
    }

    // calculate statistical error
    for (int n = 0; n < nalpha; ++n) {
        err_H_alpha[n] = pow(err_H_alpha[n] - pow(H_alpha[n], 2), 0.5) / pow(nbin - 1, 0.5);
        for (int i = 0; i < nconfig; ++i) {
            err_n_x_alpha[n][i] = pow(err_n_x_alpha[n][i] - pow(n_x_alpha[n][i], 2), 0.5) / pow(nbin - 1, 0.5);
        }
    }
}

