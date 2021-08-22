#ifndef SAC_H
#define SAC_H
#pragma once

/**
 *  This head file includes SAC class for the implementation of stochastic analytic continuation method,
 *  proposed by Anders.W. Sandvik.
 *  A Monte Carlo process and simulated annealing are performed
 *  to extract real-frequency information from imaginary-time correlations,
 *  which are obtained previously by QMC calculations.
 *
 */

#include <random>

#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#include <Eigen/Core>

#include "ReadInModule.h"
#include "FrequencyGrid.h"
#include "Kernel.h"

// random engine
static std::default_random_engine rand_engine_sac(time(nullptr));


class SAC {
public:

    /* model params */
    int nt{};                           // number of time slices
    double beta{};                      // inverse temperature
    double scale_factor{};              // scaling factor G(0)
    int nbootstrap{};                   // number of bootstrap samples

    /* griding params */
    FrequencyGrid grid;

    /* sampling params */
    int ndelta{};                       // number of delta functions
    int window_width{};                 // width of random move window
    double delta_amplitude{};           // amplitude of delta functions
    std::vector<int> delta_locations;   // locations of delta functions

    double theta{};                     // sampling temperature
    double accept_radio{};              // accepting radio of updates

    Kernel kernel;                      // kernel
    Eigen::VectorXd tau;                // tau points
    Eigen::VectorXd corr;               // time correlations from QMC
    Eigen::VectorXd sigma;              // standard deviation of transformed correlations

    Eigen::VectorXd corr_current;       // current time correlations from spectrum
    Eigen::VectorXd corr_update;        // updated time correlations from spectrum

    ReadInModule *readin;                // data from QMC input

public:

    SAC();

    /** subroutine for parameter settings */
    /* set up parameters for read in module */
    void set_read_in_params(int lt, double beta, int nbin, int rebin_pace, int num_bootstrap);

    /* set up file which contains data of tau points */
    void set_filename_tau(const std::string &infile_tau);

    /* set up file which contains data of time correlations */
    void set_filename_corr(const std::string &infile_corr);

    /* set up parameters for grids of frequency domain */
    void set_griding_params(double grid_interval, double spec_interval, double omega_min, double omega_max);

    /* set up parameters for sampling procedure */
    void set_sampling_params(double ndelta, double theta);

    /** initialization */
    void init();


private:

    /* read QMC data (transformed) from read in module */
    void init_from_module();

    /* initialize spectrum */
    void init_spectrum();

};

#endif //SAC_H
