#include <cassert>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "MonteCarloSAC.h"
#include "ProgressBar.hpp"


void MonteCarloSAC::set_SAC_params(int _lt, double _beta, int _nconfig, double _omega_min, double _omega_max, int _nMoment) {
    assert( _lt > 0 && _beta > 0 );
    assert( _nconfig > 0 );
    assert( _omega_min < _omega_max );

    this->lt = _lt;
    this->beta = _beta;
    this->nconfig = _nconfig;
    this->omega_min = _omega_min;
    this->omega_max = _omega_max;
    this->nMoment = _nMoment;
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
        sac_list[n].set_SAC_params(lt, beta, nconfig, omega_min, omega_max, nMoment);
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
    measure.resize(nbin, nalpha, nconfig);
}

void MonteCarloSAC::run_Monte_Carlo() {
    /*
     *  Metropolis Monte Carlo update in parallel,
     *  swap configs of adjacent layers every n_swap_pace steps.
     */

    begin_t = std::chrono::steady_clock::now();

    // warm up process
    // progress bar
    // progresscpp::ProgressBar progress_bar_warm(nwarm, 40, '#', '-');

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

        //++progress_bar_warm;
        if ( warm % 10 == 0 ) {
            // std::cout << "Warm-up progress:   ";
            //progress_bar_warm.display();
        }
    }
    // std::cout << "Warm-up progress:   ";
    // progress_bar_warm.done();

    // clear measuring data previously
    measure.clear_tmp_stats();

    // measuring and loop for bins
    // progresscpp::ProgressBar progress_bar_measure(nbin * nstep_1bin, 40, '#', '-');

    for (int bin = 0; bin < nbin; ++bin) {

        for (int step = 0; step < nstep_1bin; ++step) {

            #pragma omp parallel for num_threads(4) default(none)
            for (int n = 0; n < nalpha; ++n) {
                sac_list[n].Metropolis_update_1step();
            }
            measure.measure(*this);

            if ( ( step + 1 ) % n_swap_pace == 0 ) {
                swap_configs_between_layers();

                // update for some steps shortly after swapping
                for (int i = 0; i < step_between_bins; ++i) {
                    #pragma omp parallel for num_threads(4) default(none)
                    for (int n = 0; n < nalpha; ++n) {
                        sac_list[n].Metropolis_update_1step();
                    }
                }
            }

            // ++progress_bar_measure;
            if ( step % 10 == 0 ) {
                // std::cout << "Measuring progress: ";
                // progress_bar_measure.display();
            }
        }

        measure.normalize_stats();

        measure.write_data_to_bin(bin);

        measure.clear_tmp_stats();

        // avoid correlation between bins
        for (int i = 0; i < step_between_bins; ++i) {
            #pragma omp parallel for num_threads(4) default(none)
            for (int n = 0; n < nalpha; ++n) {
                sac_list[n].Metropolis_update_1step();
            }
        }
    }
    // std::cout << "Measuring progress: ";
    // progress_bar_measure.done();

    measure.analyse_stats();

    end_t = std::chrono::steady_clock::now();
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

        if ( n == 14 ) {
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

void MonteCarloSAC::print_stats() const {

    // calculate cpu time
    auto time = std::chrono::duration_cast<std::chrono::milliseconds>(end_t - begin_t).count();
    const int minute = std::floor((double)time / 1000 / 60);
    const double sec = (double)time / 1000 - 60 * minute;

    // Todo: print data

    // print cpu time
    std::cout << "Time Cost:    " << minute << " min " << sec << " s" << std::endl;
}

void MonteCarloSAC::output_stats(const std::string &outfilename) {
    std::ofstream outfile;
    outfile.open(outfilename, std::ios::out | std::ios::trunc);

    for (int n = 0; n < nalpha; ++n) {
        outfile << std::setiosflags(std::ios::right)
                << std::setw(15) << n
                << std::setw(15) << log(measure.H_alpha[n])
                << std::setw(15) << measure.err_H_alpha[n] / measure.H_alpha[n]
                << std::endl;
    }
    outfile.close();
}
