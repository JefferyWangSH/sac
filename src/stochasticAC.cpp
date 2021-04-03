#include "stochasticAC.h"
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <fstream>
#include <iostream>
#include <vector>

void StochasticAC::set_SAC_params(int lt, double beta, int nOmega, double omegaMin, double omegaMax) {
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

void StochasticAC::read_QMC_data(const std::string& filename) {
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
    // perform list of omega
    for (int i = 0; i < nOmega; ++i) {
        omega_list(i) = omegaMin + i * deltaOmega;
    }

    // calculate kernel matrix
    for (int t = 0; t < lt; ++t) {
        const int tau = tau_list(t);
        for (int i = 0; i < nOmega; ++i) {
            const int omega = omega_list(i);
            KernelMat(t, i) = (exp(-tau*omega) + exp(-(beta-tau)*omega)) / M_PI;
        }
    }

    // TODO: initialize G_tau and chi^2
}

void StochasticAC::Metropolis_update() {

}

