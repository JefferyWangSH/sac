#ifndef STOCHASTIC_ANALYTIC_CONTINUATION_MONTECARLOSAC_H
#define STOCHASTIC_ANALYTIC_CONTINUATION_MONTECARLOSAC_H
#pragma once

/*
 *  This head file includes MonteCarloSAC class
 *  to obtain spectrum functions of statistical significance
 *  by standard MonteCarlo method with Metropolis algorithm.
 */

#include <ctime>

#include "SAC.h"


class MonteCarloSAC {

public:
    // measuring parameters
    int nbin = 20;
    int nBetweenBins = 5;
    int nstep = 100;
    int nwarm = (int)pow(10, 5);

    // SAC class and sampling parameters
    SAC sac;
    double theta = exp(7);
    int nCst = 4;

    // input filename
    std::string infile_Green = "../results/benchmark_g.txt";
    std::string infile_A;

    // record cpu time
    time_t begin_t = clock(), end_t = clock();

    // measuring quantities
    std::vector<double> binEntropy;
    double meanEntropy = 0.0;
    double errEntropy = 0.0;

    std::vector<double> binChi2;
    double meanChi2 = 0.0;
    double errChi2 = 0.0;

    bool is_read_data = false;

    MonteCarloSAC() = default;

    ~MonteCarloSAC();

    /* set parameters */
    void set_SAC_params(int lt, double beta, int nOmega, double omegaMin, double omegaMax);

    void set_meas_params(int nbin, int nBetweenBins, int nstep, int nwarm);

    void set_sampling_params(const double &theta, const int &nCst);

    void set_input_file(const std::string &infile_Green, const std::string &infile_A = "");

    /* prepare for measuring: reading data from input file */
    void prepare();

    /* measuring process */
    void measure();

    /* analyse data and output results */
    void analyse_Stats();

    void print_Stats() const;

    void output_Stats(const std::string &outfilename) const;

    void output_Config(const std::string &outfilename) const;

    void clear_Stats();

    double calculate_Entropy();

private:


};

#endif //STOCHASTIC_ANALYTIC_CONTINUATION_MONTECARLOSAC_H
