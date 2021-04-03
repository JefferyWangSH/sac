#ifndef STOCHASTIC_ANALYTIC_CONTINUATION_STOCHASTICAC_H
#define STOCHASTIC_ANALYTIC_CONTINUATION_STOCHASTICAC_H
#pragma once

/*
 *  This head file includes StochasticAC class to analyse imaginary-time DQMC data.
 *  With method of Stochastic Analytic Continuation (SAC), we managed to
 *  obtain the real frequency fermion spectrum functions from imaginary-time Matsubara Green's functions.
 */

#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2

#include <Eigen/Core>

typedef Eigen::MatrixXd matXd;
typedef Eigen::VectorXd vecXd;

class StochasticAC
{
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

    StochasticAC() = default;

    /* set up params for SAC */
    void set_SAC_params(int lt, double beta, int nOmega, double omegaMin, double omegaMax);

    /* read DQMC data of dynamic measurements from file  */
    void read_QMC_data(const std::string& filename);

    /* initialize information for simulation */
    void initialSAC();

    /* calculate chi^2: the deviation of G(\tau), for a specific weight config A(i) */
    double cal_chi_square();

    /* MC update according to Metropolis algorithm */
    void Metropolis_update();


};

#endif //STOCHASTIC_ANALYTIC_CONTINUATION_STOCHASTICAC_H
