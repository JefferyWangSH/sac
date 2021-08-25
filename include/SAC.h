#ifndef STOCHASTIC_ANALYTIC_CONTINUATION_SAC_H
#define STOCHASTIC_ANALYTIC_CONTINUATION_SAC_H
#pragma once

/*
 *  This head file includes SAC class to extract fermion spectrum information from imaginary-time DQMC data.
 *  With method of Stochastic Analytic Continuation (SAC), we managed to
 *  obtain the real frequency fermion spectrum or LDOS from imaginary-time Matsubara Green's functions.
 *  Reference:
 *      K. S. D. Beach, arXiv:cond-mat/0403055
 */

#include <random>

// random engine
static std::default_random_engine generate_SAC(time(nullptr));


class SAC {

public:
    int lt = 80;
    double beta = 4.0;
    double dtau = 0.05;

    // slices in configuration space, eventually map to energy space
    // FIXME: mind here not include 0 and 1
    int nconfig = 50;
    double omega_min = -5.0;
    double omega_max = 5.0;
    double delta_n_config = 1.0 / (50 + 1);

    // vectors and matrices for SAC
    std::vector<double> x_list, n_list, tau_list;
    std::vector<double> omega_list, A_config;
    std::vector<double> g_tau, sigma_tau;
    std::vector<std::vector<double>> kernel;

    // Hamiltonian density h(\tau)
    std::vector<double> h_tau;

    // fictitious inverse temperature
    double alpha = exp(-10);

    // sampling parameter to accelerate convergence
    // conserve the first ( n_moment ) moments when update configurations
    int nMoment = 1;

    // name of file which contains QMC data
    std::string filename_greens = "default.txt";
    std::string filename_configs;

    /*
    // params for calculate accepting or swap rate of configurations
    int total_step = 0;
    int accept_step = 0;
    double accept_rate = 0.0;
    */

    SAC() = default;

    /* set up params for SAC */
    void set_SAC_params(int lt, double beta, int nconfig, double omega_min, double omega_max, int nMoment);

    /* set up inverse temperature alpha */
    void set_alpha(const double& alpha);

    /* set up filename for reading QMC data*/
    void set_QMC_filename(const std::string& filename);

    void set_Configs_filename(const std::string& filename);

    /* read DQMC data of dynamic measurements from file  */
    void read_QMC_data(const std::string& filename);

    void read_Configs_data(const std::string& filename = "");

    /* prepare for simulation, including reading data from input file */
    virtual void prepare();

    /* perform one step of updates */
    void Metropolis_update_1step();

    /* calculate Hamiltonian for a specific configuration */
    void cal_Config_Hamiltonian(double& H);

    friend class measure;

public:

    /* calculate Hamiltonian density h(\tau) for a specific configuration */
    void cal_Config_Hamiltonian_Density(std::vector<double> &h);

    /* randomly update weight configs */
    void update_Configs(std::vector<int> &index_selected, std::vector<double> &n_selected);

    /* MC update according to Metropolis algorithm */
    void Metropolis_update();

};

#endif //STOCHASTIC_ANALYTIC_CONTINUATION_SAC_H
