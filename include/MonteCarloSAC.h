#ifndef STOCHASTIC_ANALYTIC_CONTINUATION_MONTECARLOSAC_H
#define STOCHASTIC_ANALYTIC_CONTINUATION_MONTECARLOSAC_H
#pragma once

/*
 *  This head file includes MonteCarloSAC class.
 *    This is the critical part for our SAC calculation,
 *    and always cooperates with a simulated annealing procedure.
 */

#include <ctime>
#include <chrono>
#include <random>

// random engine
static std::default_random_engine gen_MC_SAC(time(nullptr));

#include "SAC.h"
#include "Measure.h"


class MonteCarloSAC : public SAC {

public:

    // measuring parameters
    int nbin = 20;
    int nstep_1bin = 100;
    int step_between_bins = 10;

    int nwarm = 1e5;

    // Measure class
    Measure measure;

    // record cpu time
    std::chrono::steady_clock::time_point begin_t;
    std::chrono::steady_clock::time_point end_t;

public:

    MonteCarloSAC() = default;

    /* set up parameters for measurements */
    void set_measure_params(int nbin, int nstep_1bin, int step_between_bins, int nwarm);

    /* prepare for SAC and measurements */
    void prepare() override;

    /* Monte Carlo procedure */
    void run_Monte_Carlo();

    /* print out statistical data onto the terminal */
    void print_stats() const;

    /* file output of the goodness of fit */
    void output_stats(const std::string& filename) const;

    /* file output of averaged field configurations n(x) */
    void output_Configs(const std::string& filename) const;

    /* file output of the recovered spectrum */
    void output_averaged_spectrum(const std::string& filename) const;

};

#endif //STOCHASTIC_ANALYTIC_CONTINUATION_MONTECARLOSAC_H
