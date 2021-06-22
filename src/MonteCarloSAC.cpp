#include <cassert>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "MonteCarloSAC.h"
#include "ProgressBar.hpp"


void MonteCarloSAC::set_measure_params(int _nbin, int _nstep_1bin, int _step_between_bins, int _nwarm) {
    this->nbin = _nbin;
    this->nstep_1bin = _nstep_1bin;
    this->step_between_bins = _step_between_bins;
    this->nwarm = _nwarm;
}

void MonteCarloSAC::prepare() {
    /*
     *  allocate memory for SAC and initialize temp vector and matrices for simulation use,
     */

    // clear previous data
    x_list.clear();
    n_list.clear();
    omega_list.clear();
    A_config.clear();
    tau_list.clear();
    g_tau.clear();
    sigma_tau.clear();
    h_tau.clear();
    kernel.clear();

    // reserve memory for vectors
    x_list.reserve(nconfig);
    n_list.reserve(nconfig);
    omega_list.reserve(nconfig);
    A_config.reserve(nconfig);

    for (int i = 0; i < nconfig; ++i){
        x_list.emplace_back(0.0);
        n_list.emplace_back(0.0);
        omega_list.emplace_back(0.0);
        A_config.emplace_back(0.0);
    }

    tau_list.reserve(lt);
    g_tau.reserve(lt);
    sigma_tau.reserve(lt);
    h_tau.reserve(lt);
    kernel.reserve(lt);

    for (int t = 0; t < lt; ++t) {
        tau_list.emplace_back(0.0);
        g_tau.emplace_back(0.0);
        sigma_tau.emplace_back(0.0);
        h_tau.emplace_back(0.0);
        kernel.emplace_back(nconfig, 0.0);
    }

    // shrink to fit
    x_list.shrink_to_fit();
    n_list.shrink_to_fit();
    omega_list.shrink_to_fit();
    A_config.shrink_to_fit();
    tau_list.shrink_to_fit();
    g_tau.shrink_to_fit();
    sigma_tau.shrink_to_fit();
    h_tau.shrink_to_fit();
    kernel.shrink_to_fit();

    // initialize tau
    for (int t = 0; t < lt; ++t) {
        tau_list[t] = (double)(t + 1) / lt * beta;
    }

    // initialize x and omega
    for (int i = 0; i < nconfig; ++i) {
        x_list[i] = delta_n_config * (i + 1);
        omega_list[i] = omega_min + ( omega_max - omega_min ) * x_list[i];
    }

    read_QMC_data(filename_greens);
    read_Configs_data(filename_configs);

    // calculate kernel matrix
    for (int t = 0; t < lt; ++t) {
        const double tau = tau_list[t];
        for (int i = 0; i < nconfig; ++i) {
            const double omega = omega_list[i];
            kernel[t][i] = exp(- tau * omega) / (1 + exp(- beta * omega));
        }
    }

    // calculate hamiltonian density h(\tau) for default configurations
    cal_Config_Hamiltonian_Density(h_tau);

    // prepare for measuring
    measure.resize(nbin, nconfig);
}

void MonteCarloSAC::run_Monte_Carlo() {
    /*
     *  Metropolis Monte Carlo update in parallel,
     *  swap configs of adjacent layers every n_swap_pace steps.
     */

    begin_t = std::chrono::steady_clock::now();

    // warm up process
    // progress bar
    progresscpp::ProgressBar progress_bar_warm(nwarm, 40, '#', '-');

    for (int warm = 0; warm < nwarm; ++warm) {

        // one Monte Carlo step
        Metropolis_update_1step();

        ++progress_bar_warm;
        if ( warm % 5 == 0 ) {
            std::cout << " Warm-up progress:   ";
            progress_bar_warm.display();
        }
    }
    std::cout << " Warm-up progress:   ";
    progress_bar_warm.done();

    // clear measuring data previously
    measure.clear_tmp_stats();

    // measuring and loop for bins
    progresscpp::ProgressBar progress_bar_measure(nbin * nstep_1bin, 40, '#', '-');

    for (int bin = 0; bin < nbin; ++bin) {

        for (int step = 0; step < nstep_1bin; ++step) {

            Metropolis_update_1step();

            measure.measure(*this);

            ++progress_bar_measure;
            if ( step % 5 == 0 ) {
                std::cout << " Measuring progress: ";
                progress_bar_measure.display();
            }
        }

        measure.normalize_stats();

        measure.write_data_to_bin(bin);

        measure.clear_tmp_stats();

        // avoid correlation between bins
        for (int i = 0; i < step_between_bins; ++i) {
            Metropolis_update_1step();
        }
    }
    std::cout << " Measuring progress: ";
    progress_bar_measure.done();

    measure.analyse_stats();

    end_t = std::chrono::steady_clock::now();
}

void MonteCarloSAC::print_stats() const {

    // calculate cpu time
    auto time = std::chrono::duration_cast<std::chrono::milliseconds>(end_t - begin_t).count();
    const int minute = std::floor((double)time / 1000 / 60);
    const double sec = (double)time / 1000 - 60 * minute;

    std::cout << std::setiosflags(std::ios::left)
              << "  ln ( \\alpha ): " << std::setw(8) << log(alpha)
              << "  Goodness of fit: " << std::setw(8) << log(measure.mean_Hamilton)
              << std::endl;

    // print cpu time
    std::cout << "  Time Cost:    " << minute << " min " << sec << " s" << std::endl;
}

void MonteCarloSAC::output_stats(const std::string& filename) const {
    std::ofstream outfile;
    outfile.open(filename, std::ios::out | std::ios::app);

    outfile << std::setiosflags(std::ios::right)
            << std::setw(15) << log(alpha)
            << std::setw(15) << measure.mean_Hamilton
            << std::setw(15) << measure.err_Hamilton
            << std::setw(15) << log(measure.mean_Hamilton)
            << std::setw(15) << measure.err_Hamilton / measure.mean_Hamilton
            << std::endl;
    outfile.close();
}

void MonteCarloSAC::output_Configs(const std::string& filename) const {
    std::ofstream outfile;
    outfile.open(filename, std::ios::out | std::ios::trunc);

    for (int i = 0; i < nconfig; ++i) {
        outfile << std::setiosflags(std::ios::right)
                << std::setw(15) << i + 1
                << std::setw(15) << measure.mean_n_x[i]
                << std::setw(15) << measure.err_n_x[i]
                << std::endl;
    }
    outfile.close();
    std::cout << "=====================================================================" << std::endl
              << " configs have been written into " + filename + "." << std::endl
              << "=====================================================================" << std::endl;
    std::cout << std::endl;
}

void MonteCarloSAC::output_averaged_spectrum(const std::string& filename) const {
    std::ofstream outfile;
    outfile.open(filename, std::ios::out | std::ios::trunc);

    for (int i = 0; i < nconfig; ++i) {
        outfile << std::setiosflags(std::ios::right)
                << std::setw(15) << omega_list[i]
                << std::setw(15) << measure.mean_n_x[i]
                << std::setw(15) << measure.err_n_x[i]
                << std::endl;
    }
    outfile.close();
}
