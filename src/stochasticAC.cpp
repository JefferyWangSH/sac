#include "stochasticAC.h"
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <fstream>
#include <iostream>
#include <vector>

void StochasticAC::set_SAC_params(int lt, double beta, int nOmega, double omegaMin, double omegaMax) {
    /*
     *  Set SAC parameters.
     */

    // we restrict ourselves to the subspace where \omega > 0
    assert(lt > 0 && beta > 0 && nOmega > 0);
    assert(omegaMin >= 0);
    assert(omegaMin < omegaMax);

    this->lt = lt;
    this->beta = beta;
    this->dtau = beta/lt;
    this->nOmega = nOmega;
    this->omegaMin = omegaMin;
    this->omegaMax = omegaMax;
    this->deltaOmega = (omegaMax - omegaMin) / nOmega;

    omega_list.resize(nOmega);
    tau_list.resize(lt);
    A_omega.resize(nOmega);
    g_tau.resize(lt);
    g_tau_QMC.resize(lt);
    err_g_tau_QMC.resize(lt);

    KernelMat.resize(lt, nOmega);
}

void StochasticAC::set_sampling_params(const double &theta, const int &n_constraint) {
    this->theta = theta;
    this->n_constraint = n_constraint;
}

void StochasticAC::read_QMC_data(const std::string& filename) {
    /*
     *  Read imaginary-time QMC data from file.
     */

    std::ifstream infile;
    infile.open(filename, std::ios::in);

    if (!infile.is_open()) {
        std::cerr << "Fail to open file " + filename + " !"<<  std::endl;
    }

    std::string line;
    int i = 0;
    while(getline(infile, line)) {
        /* call boost library */       std::vector<std::string> data;
        boost::split(data, line, boost::is_any_of(" "), boost::token_compress_on);
        tau_list(i) = boost::lexical_cast<double>(data[0]);
        g_tau_QMC(i) = boost::lexical_cast<double>(data[1]);
        err_g_tau_QMC(i) = boost::lexical_cast<double>(data[2]);
        ++i;
    }
    assert(i == lt);
    infile.close();
}

void StochasticAC::initialSAC() {
    /*
     *  Initialize temp vector and matrices for simulation use.
     */

    // initialize list of omega and A(\omega)
    for (int i = 0; i < nOmega; ++i) {
        omega_list(i) = omegaMin + (i + 1) * deltaOmega;

        // initialize A(\omega) as a flat spectrum, thus normalized condition \sum A(\omega) * deltaOmega = 1.
        // suppose system is invariant under space transformation, we here restrict omega to be positive.
        A_omega(i) = 1 / (2 * deltaOmega);
    }

    // calculate kernel matrix
    for (int t = 0; t < lt; ++t) {
     const int tau = tau_list(t);
        for (int i = 0; i < nOmega; ++i) {
            const int omega = omega_list(i);
            KernelMat(t, i) = (exp(-tau*omega) + exp(-(beta-tau)*omega)) / M_PI;
        }
    }

    // initialize the fitting Green function g_tau and chi ^ 2
    cal_chi_square(A_omega, g_tau, chi_square);
}

void StochasticAC::cal_chi_square(const vecXd& A, vecXd& g_tau_fitted, double& chi_2) {
    /*
     *  Calculate chi ^ 2 and the fitting g(\tau) from a given weight config A(\omega)
     *  definition: chi ^ 2 = \sum ( g(\tau) - g_QMC(\tau) )^2 / \delta(\tau) ^ 2
     *          g (\tau) = \sum_{\omega} Kernel(\tau, \omega) * A(\omega)
     */

    // update the fitting g(\tau)
    g_tau_fitted = KernelMat * A;

    // calculate chi ^ 2
    chi_2 = 0.0;
    for (int t = 0; t < lt; ++t) {
        chi_2 += pow((g_tau_fitted(t) - g_tau_QMC(t)), 2) / pow(err_g_tau_QMC(t), 2);
    }
}

void StochasticAC::Metropolis_update() {
    /*
     *  Local update of weight configuration A(\omega) via Metropolis algorithm.
     *  Configurations are updated in place, according a probability distribution proportional to
     *          P (A) \propto exp(- chi^2 / \theta)
     *  Update are performed in such a way that the first few ( n_constraint ) frequency moments are conserved.
     */

}
