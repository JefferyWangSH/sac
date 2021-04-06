#include <iostream>

#include "measure.h"

void Measure::set_SAC_params(int lt, double beta, int nOmega, double omegaMin, double omegaMax) {
    this->sac.set_SAC_params(lt, beta, nOmega, omegaMin, omegaMax);
}

void Measure::set_meas_params(int nbin, int nBetweenBins, int nstep, int nwarm) {
    this->nbin = nbin;
    this->nBetweenBins = nBetweenBins;
    this->nstep = nstep;
    this->nwarm = nwarm;
}

void Measure::set_sampling_params(const double &theta, const int &nCst) {
    this->theta = theta;
    this->nCst = nCst;
    this->sac.set_sampling_params(theta, nCst);
}

void Measure::set_file_name(const std::string &infilename, const std::string &outfilename) {
    this->infilename = infilename;
    this->outfilename = outfilename;
    is_data_read = false;
}

void Measure::prepare() {
    sac.read_QMC_data(infilename);
    is_data_read = true;
    sac.initialSAC();   // pre-read of data is needed for process of initialization
}

