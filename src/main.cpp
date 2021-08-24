#include <iostream>
#include <fstream>
#include <iomanip>

#include <boost/program_options.hpp>

#include "SAC.h"

#include <chrono>


/**
  *  TODO:
  *   1. Optimize random module (missing)
  *   2. ...
  */


/** The main program */
int main(int argc, char *argv[]) {

    int lt = 100;
    double beta = 4.0;
    int nbin = 1000;
    int rebin_pace = 1;
    int num_boostrap = 5000;

    double grid_interval = 1e-5;
    double spec_interval = 1e-2;
    double omega_min = -10.0;
    double omega_max = 10.0;

    int ndelta = 1000;
    double theta = 200.0;
    int max_annealing_steps = 5e3;
    int bin_size = 8e3;
    int bin_num = 30;

    std::chrono::steady_clock::time_point begin_t{}, end_t{};

    begin_t = std::chrono::steady_clock::now();

    Simulation::SAC *sac = new Simulation::SAC();

    sac->set_read_in_params(lt, beta, nbin, rebin_pace, num_boostrap);
    sac->set_filename_tau("../input/tau.dat");
    sac->set_filename_corr("../input/corr.dat");
    sac->set_griding_params(grid_interval, spec_interval, omega_min, omega_max);
    sac->set_sampling_params(ndelta, theta, max_annealing_steps, bin_num, bin_size);
    sac->set_mode_params("fermion", "single");

    sac->init();

    sac->annealing();

    end_t = std::chrono::steady_clock::now();
    std::cout << "total cost "
              << (double)std::chrono::duration_cast<std::chrono::milliseconds>(end_t-begin_t).count()/1000
              << "s." << std::endl;

    delete sac;
    return 0;
}