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

    binChi2.resize(nbin);
    binEntropy.resize(nbin);
}

void Measure::set_sampling_params(const double &theta, const int &nCst) {
    this->theta = theta;
    this->nCst = nCst;
    this->sac.set_sampling_params(theta, nCst);
}

void Measure::set_input_file(const std::string &infilename) {
    this->infilename = infilename;
    is_data_read = false;
}

void Measure::prepare() {
    sac.read_QMC_data(infilename);
    is_data_read = true;
    sac.initialSAC();   // pre-read of data is needed for process of initialization
}

void Measure::measure(bool bool_display_process) {

}

void Measure::print_Stats() const {
    const double time = (double)(end_t - begin_t)/CLOCKS_PER_SEC;
    const int minute = floor(time / 60);
    const double sec = time - 60 * minute;

    std::cout << "===================================================================" << std::endl;

    std::cout << "Simulation Parameters: " << std::endl
              << "ln(1 / \\Theta):     " << log(1/theta) << std::endl
              << "nCst:                " << nCst << std::endl
              << std::endl;

    std::cout << "Measurements: " << std::endl
              << "Entropy S:           " << meanEntropy << "    err: " << errEntropy << std::endl
              << "Chi square chi^2:    " << meanChi2 << "    err: " << errChi2 << std::endl
              << std::endl;

    std::cout << "Time Cost:      " << minute << " min " << sec << " s" << std::endl;

    std::cout << "===================================================================" << std::endl
              << std::endl;

}

void Measure::output_Stats(const std::string &outfilename) {

}

void Measure::clear_Stats() {
    assert(binEntropy.size() == nbin);
    assert(binChi2.size() == nbin);

    for (int bin = 0; bin < nbin; ++bin) {
        binEntropy[bin] = 0.0;
        binChi2[bin] = 0.0;
    }

    meanEntropy = 0.0;
    errEntropy = 0.0;
    meanChi2 = 0.0;
    errChi2 = 0.0;
}

void Measure::calculate_Entropy() {

}
