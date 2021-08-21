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

#include "ReadInModule.h"

#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#include <Eigen/Core>

// random engine
static std::default_random_engine rand_engine_sac(time(nullptr));


class SAC {
public:

    /* model params */
    int lt{};                           // number of time slices
    double beta{};                      // inverse temperature

    /* sampling params */
    int num_delta{};                    // number of delta functions

    /* griding params */
    double omega_min{};                 // min and max of frequency
    double omega_max{};
    int int_omega_min{};                // min and max of frequency number
    int int_omega_max{};
    double fine_delta_omega{};          // delta frequency of sampling space
    double spec_delta_omega{};          // delta frequency of spectrum (frequency spacing in accumulated histogram)

    double theta{};                     // sampling temperature
    double accept_radio{};              // accepting radio of updates

    Eigen::MatrixXd kernel;             // kernel
    Eigen::VectorXd corr;               // current time correlations from spectrum
    Eigen::VectorXd corr_updated;       // updated time correlations from spectrum

    ReadInModule readin;                // data from QMC input

public:

    SAC() = default;

    void set_read_in_params(int lt, double beta, int nbin, int rebin_pace, int num_bootstrap);

    void set_filename_tau(const std::string &infile_tau);

    void set_filename_corr(const std::string &infile_corr);

    void set_griding_params(double spec_delta_omega, double fine_delta_omega, double omega_min, double omega_max);

    void set_sampling_params(double num_delta);

    void init();


private:

    void init_spectrum();

    void init_kernel();

};

#endif //SAC_H
