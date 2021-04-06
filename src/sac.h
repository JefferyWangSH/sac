#ifndef STOCHASTIC_ANALYTIC_CONTINUATION_SAC_H
#define STOCHASTIC_ANALYTIC_CONTINUATION_SAC_H
#pragma once

/*
 *  This head file includes SAC class to analyse imaginary-time DQMC data.
 *  With method of Stochastic Analytic Continuation (SAC), we managed to
 *  obtain the real frequency fermion spectrum functions from imaginary-time Matsubara Green's functions.
 */

#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2

#include <Eigen/Core>
#include <random>

typedef Eigen::MatrixXd matXd;
typedef Eigen::VectorXd vecXd;

// random engine
static std::default_random_engine gen(time(nullptr));


class SAC {

public:
    int lt = 80;
    double beta = 4.0;
    double dtau = 0.05;

    int nOmega = 50;
    double omegaMin = 0.0, omegaMax = 5.0;
    double deltaOmega = 0.1;

    vecXd omega_list, tau_list;
    vecXd A_omega, g_tau;
    vecXd g_tau_QMC, err_g_tau_QMC;
    matXd KernelMat;
    double chi_square = 0.0;

    double theta = exp(7.0);    // sampling temperature of SAC system, act as regularization parameter
    int n_constraint = 4;          // sampling parameter to accelerate convergence


    SAC() = default;

    /* set up params for SAC */
    void set_SAC_params(int lt, double beta, int nOmega, double omegaMin, double omegaMax);

    /* set up sampling temperature \theta */
    void set_sampling_params(const double &theta, const int &n_constraint);

    /* read DQMC data of dynamic measurements from file  */
    void read_QMC_data(const std::string &filename);

    /* prepare for simulation (enquire reading QMC data first) */
    void initialSAC();

    /* perform one step of updates */
    void Metropolis_update_1step();

private:
    /* calculate chi ^ 2, for a specific weight configuration A(\omega) */
    void cal_chi_square(const vecXd &A, vecXd &g_tau_fitted, double &chi_2);

    /* randomly update weight configs */
    void update_A_omega(vecXd &A_omega_old);

    /* MC update according to Metropolis algorithm */
    void Metropolis_update();

};

#endif //STOCHASTIC_ANALYTIC_CONTINUATION_SAC_H
