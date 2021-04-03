#include "stochasticAC.h"
#include "strings.h"
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
    spectrum_omega.resize(nOmega);
    green_tau.resize(lt);
    green_tau_QMC.resize(lt);
    err_tau_QMC.resize(lt);

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
        std::vector<std::string> data;
        split(line, data, " ");
        tau_list(i) = str2double(data[0]);
        green_tau_QMC(i) = str2double(data[1]);
        err_tau_QMC(i) = str2double(data[2]);
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

