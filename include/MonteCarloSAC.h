#ifndef STOCHASTIC_ANALYTIC_CONTINUATION_MONTECARLOSAC_H
#define STOCHASTIC_ANALYTIC_CONTINUATION_MONTECARLOSAC_H
#pragma once

/*
 *  This head file includes MonteCarloSAC class
 *  This is the critical part for our SAC calculation.
 *  By method of parallel tempering and standard Metropolis Monte Carlo algorithm,
 *  we calculate SACs for large range of temperature profiles in parallel to accelerate simulation.
 */

#include <ctime>
#include <chrono>
#include <random>

// random engine
static std::default_random_engine gen_MC_SAC(time(nullptr));

#include "SAC.h"
#include "Measure.h"


class MonteCarloSAC {

public:

    // temperature ( alpha ) parameters
    int nalpha = 60;
    std::vector<double> alpha_list;

    // SACs in different temperature slices and SAC params
    std::vector<SAC> sac_list;
    int lt = 80;
    double beta = 4.0;
    int nconfig = 50;
    double omega_min = -5;
    double omega_max = 5;
    int nMoment = 1;

    // pace of swap
    int n_swap_pace = 500;

    // measuring parameters
    int nbin = 20;
    int nstep_1bin = 100;
    int step_between_bins = 10;

    int nwarm = (int) pow(10, 5);

    Measure measure;

    std::vector<std::vector<double>> p_list;

    // input filename
    std::string infile_green = "default.txt";
    std::string infile_A;

    // record cpu time
    std::chrono::steady_clock::time_point begin_t;
    std::chrono::steady_clock::time_point end_t;

    bool is_read_data = false;

    MonteCarloSAC() = default;

    /* set parameters */
    void set_SAC_params(int lt, double beta, int n_config, double omega_min, double omega_max, int n_moment);

    void set_measure_params(int nbin, int nstep_1bin, int step_between_bins, int nwarm, int n_swap_pace);

    void set_tempering_profile(const int& nalpha, const std::vector<double> &alpha_list);

    void set_input_file(const std::string &infile_green, const std::string &infile_A = "");

    /* prepare for measuring: reading data from input file */
    void prepare();

    /* Monte Carlo process */
    void run_Monte_Carlo();

    /* Swap configurations between adjacent temperature slices */
    void swap_configs_between_layers();

    void print_stats() const;

    void output_stats(const std::string &outfilename);

    void output_Config(const std::string &outfilename) const;


private:


};

#endif //STOCHASTIC_ANALYTIC_CONTINUATION_MONTECARLOSAC_H
