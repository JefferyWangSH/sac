#include <cassert>
#include <iostream>

#include "MonteCarloSAC.h"


void MonteCarloSAC::set_SAC_params(int _lt, double _beta, int _n_config, double _omega_min, double _omega_max, int _n_moment) {
    assert( _lt > 0 && _beta > 0 );
    assert( _n_config > 0 );
    assert( _omega_min < _omega_max );

    this->lt = _lt;
    this->beta = _beta;
    this->n_config = _n_config;
    this->omega_min = _omega_min;
    this->omega_max = _omega_max;
    this->n_moment = _n_moment;
}

void MonteCarloSAC::set_tempering_profile(const int &_nalpha, const std::vector<double> &_alpha_list) {
    assert( _alpha_list.size() == _nalpha );

    this->nalpha = _nalpha;

    // initialize alpha profiles
    this->alpha_list.clear();
    this->alpha_list.reserve(_nalpha);
    for (int i = 0; i < _nalpha; ++i) {
        this->alpha_list.emplace_back(_alpha_list[i]);
    }
    this->alpha_list.shrink_to_fit();
}

void MonteCarloSAC::set_measure_params(int _nbin, int _nstep_1bin, int _step_between_bins, int _nwarm, int _n_swap_pace) {
    this->nbin = _nbin;
    this->nstep_1bin = _nstep_1bin;
    this->step_between_bins = _step_between_bins;
    this->nwarm = _nwarm;
    this->n_swap_pace = _n_swap_pace;
}

void MonteCarloSAC::set_input_file(const std::string &_infile_green, const std::string &_infile_A) {
    this->infile_green = _infile_green;
    this->infile_A = _infile_A;
}

void MonteCarloSAC::prepare() {
    /*
     *  preparations for parallel-tempering monte carlo procedure
     */

    // prepare for SACs
    // TODO: construct function by copy
    sac_list.clear();
    sac_list.reserve(nalpha);
    for (int n = 0; n < nalpha; ++n) {
        sac_list.emplace_back();
        sac_list[n].set_SAC_params(lt, beta, n_config, omega_min, omega_max, n_moment);
        sac_list[n].set_alpha(alpha_list[n]);
        sac_list[n].set_QMC_filename(infile_green);
        sac_list[n].prepare();
    }

    p_list.clear();
    p_list.reserve(nalpha - 1);
    for (int n = 0; n < nalpha - 1; ++n) {
        p_list.emplace_back();
    }
    p_list.shrink_to_fit();

    // prepare for measuring
    // TODO
}

void MonteCarloSAC::run_Monte_Carlo() {
    /*
     *  Metropolis Monte Carlo update in parallel,
     *  swap configs of adjacent layers every n_swap_pace steps.
     */

    // warm up process
    for (int warm = 0; warm < nwarm; ++warm) {

        // loop for different alpha slices
        #pragma omp parallel for num_threads(4) default(none)
        for (int n = 0; n < nalpha; ++n) {
            sac_list[n].Metropolis_update_1step();
        }

        // swap configs for different temperature layers
        if ((warm + 1) % n_swap_pace == 0) {
            swap_configs_between_layers();
        }
    }

    // measuring
    // TODO

}

void MonteCarloSAC::swap_configs_between_layers() {
    /*
     *  Swap the configurations of adjacent temperature layers with probability
     *      exp( ( alpha_{p} - alpha_{q} ) * ( H_{p} - H_{q} ) )
     *  if necessary, record the probability of swapping
     */

    // calculate hamiltonian for each layer
    std::vector<double> H_list(nalpha);
    for (int n = 0; n < nalpha; ++n) {
        sac_list[n].cal_Config_Hamiltonian(H_list[n]);
    }

    // loop for each pair of adjacent temperature layers and swap
    for (int n = 0; n < nalpha - 1; ++n) {

        const double p = exp( ( alpha_list[n] - alpha_list[n+1] ) * ( H_list[n] - H_list[n+1] ) );

        // FIXME:
        if ( n == 0 ) {
            std::cout << p << std::endl;
        }

        // FIXME: measurement part
        // calculate and record average accepted rate p
        if ( p_list[n].empty() ) {
            p_list[n].emplace_back(std::min(1.0, p));
        }
        else {
            const double tmp_p = ( (double)p_list[n].size() * p_list[n].back() + std::min(1.0, p) ) / ( (double)p_list[n].size() + 1 );
            p_list[n].emplace_back(tmp_p);
        }

        // if accepted
        if (std::bernoulli_distribution(std::min(1.0, p))(gen_MC_SAC)) {

            // swap configs n(x)
            sac_list[n].n_list.swap(sac_list[n+1].n_list);

            // swap hamiltonian density h(\tau)
            sac_list[n].h_tau.swap(sac_list[n+1].h_tau);

            // swap hamiltonian H
            std::swap(H_list[n], H_list[n+1]);
        }
    }
}