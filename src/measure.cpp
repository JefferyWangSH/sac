#include <iostream>
#include <fstream>
#include <iomanip>

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

void Measure::measure() {
    assert(is_data_read);

    clear_Stats();

    // record cpu time
    begin_t = clock();

    // warm up process
    for (int nwm = 0; nwm < nwarm; ++nwm) {
        sac.Metropolis_update_1step();
    }

    // measuring process
    for (int bin = 0; bin < nbin; ++bin) {
        // loop for bin
        for (int nstp = 0; nstp < nstep; ++nstp) {
            sac.Metropolis_update_1step();
            binChi2[bin] += log(sac.chi_square);
            binEntropy[bin] += calculate_Entropy();
        }

        // ensure incoherence of samples between bins
        for (int nbb = 0; nbb < nBetweenBins; ++nbb) {
            sac.Metropolis_update_1step();
        }
    }
}

void Measure::analyse_Stats() {

    // alyse data: mean and error
    for (int bin = 0; bin < nbin; ++bin) {
        meanChi2 += binChi2[bin] / nstep;
        errChi2 += pow(binChi2[bin] / nstep, 2);
        meanEntropy += binEntropy[bin] / nstep;
        errEntropy += pow(binEntropy[bin] / nstep, 2);
    }

    meanChi2 /= nbin;
    errChi2 /= nbin;
    meanEntropy /= nbin;
    errEntropy /= nbin;

    errChi2 = pow(errChi2 - meanChi2 * meanChi2, 0.5) / pow(nbin - 1, 0.5);
    errEntropy = pow(errEntropy - meanEntropy * meanEntropy, 0.5) / pow(nbin - 1, 0.5);

    sac.accept_rate = (double)sac.accept_step / (double)sac.total_step;

    end_t = clock();
}

void Measure::print_Stats() const {

    const double time = (double)(end_t - begin_t)/CLOCKS_PER_SEC;
    const int minute = floor(time / 60);
    const double sec = time - 60 * minute;

    std::cout << "====================================================" << std::endl;

    std::cout << "  Simulation Parameters: " << std::endl
              << "    Temperature ln(1/theta):      " << log(1/theta) << std::endl
              << "    Number of constraints nCst:   " << nCst << std::endl
              << std::endl;

    std::cout.precision(8);
    std::cout << "  Measurements: " << std::endl
              << "    average rate of update being accepted: " << sac.accept_rate << std::endl
              << "    ln(chi^2):    " << meanChi2 << "    err: " << errChi2 << std::endl
              << "    Entropy S:    " << meanEntropy << "    err: " << errEntropy << std::endl
              << std::endl;
    std::cout.precision(-1);

    std::cout << "  Time Cost:      " << minute << " min " << sec << " s" << std::endl;

    std::cout << "====================================================" << std::endl
              << std::endl;

}

void Measure::output_Stats(const std::string &outfilename) const {
    std::ofstream outfile;
    outfile.open(outfilename, std::ios::out | std::ios::app);

    outfile.precision(8);
    outfile << std::setiosflags(std::ios::right)
            << std::setw(15) << log(1/theta)
            << std::setw(15) << nCst
            << std::setw(15) << meanChi2
            << std::setw(15) << meanEntropy
            << std::setw(15) << errChi2
            << std::setw(15) << errEntropy
            << std::endl;
    outfile.precision(-1);
    outfile.close();
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

    sac.accept_step = 0;
    sac.total_step = 0;
    sac.accept_rate = 0.0;
}

double Measure::calculate_Entropy() {
    /*
     *  Definition of entropy S:
     *
     *      S = - \sum_{i} A(i) * ln(A(i)) * K(0, i)
     *
     *  In SAC, entropy is used to single out the most reasonable spectrum.
     */
    double s = 0.0;
    for (int i = 0; i < sac.nOmega; ++i) {
        s += - sac.A_omega(i) * log(sac.A_omega(i)) * sac.KernelMat(0, i);
    }
    return s;
}

Measure::~Measure() {
    std::cout << "The simulation was done :)" << std::endl;
}

