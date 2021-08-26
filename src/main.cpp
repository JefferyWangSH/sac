#include <chrono>
#include <fstream>
#include <iostream>
#include <boost/program_options.hpp>

#include "SAC.h"


/** The main program */
int main(int argc, char *argv[]) {

    int lt = 100;
    double beta = 4.0;
    int nbin = 1e3;
    int rebin_pace = 1;
    int nboostrap = 5e3;

    double grid_interval = 1e-5;
    double spec_interval = 1e-2;
    double omega_min = -10.0;
    double omega_max = 10.0;

    int ndelta = 1e3;
    double theta = 20;
    int max_annealing_steps = 5e3;
    int bin_size = 4e3;
    int bin_num = 5;
    int collecting_steps = 1e5;

    std::string infile_corr = "../input/cor.dat";
    std::string infile_tau = "../input/tau.dat";
    std::string outfile_spec = "../results/spec.dat";


    /** read params from command line */
    boost::program_options::options_description opts("Program options");
    boost::program_options::variables_map vm;

    opts.add_options()
            ("help,h", "display this information")
            ("lt", boost::program_options::value<int>(&lt)->default_value(100),
                    "number of time points of input QMC correlations, default: 100")
            ("beta", boost::program_options::value<double>(&beta)->default_value(4.0),
                    "inverse temperature of QMC, default: 4.0")
            ("nbin-qmc", boost::program_options::value<int>(&nbin)->default_value(1e3),
                    "total number of bins in QMC measurements, default: 1e3")
            ("rebin-pace", boost::program_options::value<int>(&rebin_pace)->default_value(1),
                    "pace of rebin, default: 1")
            ("nbootstrap", boost::program_options::value<int>(&nboostrap)->default_value(5e3),
                    "number of bootstrap samples for analysing QMC data, default: 5e3")
            ("grid-interval", boost::program_options::value<double>(&grid_interval)->default_value(1e-5),
                    "minimum interval of fine frequency grids in sampling space, default: 1e-5")
            ("spec-interval", boost::program_options::value<double>(&spec_interval)->default_value(1e-2),
                    "minimum interval of frequency grids in accumulated spectrum, default: 1e-2")
            ("freq-min", boost::program_options::value<double>(&omega_min)->default_value(-10.0),
                    "lower bound of frequency domain, default: -10.0")
            ("freq-max", boost::program_options::value<double>(&omega_max)->default_value(+10.0),
                    "upper bound of frequency domain, default: 10.0")
            ("ndelta", boost::program_options::value<int>(&ndelta)->default_value(1e3),
                    "number of delta functions, default: 1e3")
            ("theta", boost::program_options::value<double>(&theta)->default_value(20.0),
                    "initial sampling temperature, default: 20.0")
            ("max-anneal-steps", boost::program_options::value<int>(&max_annealing_steps)->default_value(5e3),
                    "maximum MC steps for simulated annealing precess, default: 5e3")
            ("sbin", boost::program_options::value<int>(&bin_size)->default_value(4e3),
                    "number of bootstrap samples within one bin, default: 4e3")
            ("nbin-sac", boost::program_options::value<int>(&bin_num)->default_value(5),
                    "total number of bins for SAC measurements, default: 5")
            ("collect-steps", boost::program_options::value<int>(&collecting_steps)->default_value(1e5),
                    "maximum MC steps for spectrum collecting precess, default: 1e5")
            ("itau", boost::program_options::value<std::string>(&infile_tau)->default_value("../input/tau.dat"),
                    "input filename of QMC time points, default: ../input/tau.dat")
            ("icor", boost::program_options::value<std::string>(&infile_corr)->default_value("../input/cor.dat"),
                    "input filename of QMC correlations data, default: ../input/cor.dat")
            ("ospec", boost::program_options::value<std::string>(&outfile_spec)->default_value("../results/spec.dat"),
                    "output filename of recovered spectrum data, default: ../results/spec.dat");

    try {
        boost::program_options::store(parse_command_line(argc, argv, opts), vm);
    }
    catch (...) {
        std::cerr << "Got undefined options from command line! "<< std::endl;
        exit(1);
    }
    boost::program_options::notify(vm);

    if (vm.count("help")) {
        std::cerr << argv[0] << std::endl;
        std::cerr << opts << std::endl;
        exit(1);
    }

    std::chrono::steady_clock::time_point begin_t{}, end_t{};
    begin_t = std::chrono::steady_clock::now();


    /** SAC simulations */
    Simulation::SAC *sac = new Simulation::SAC();

    sac->set_read_in_params(lt, beta, nbin, rebin_pace, nboostrap);
    sac->set_filename_tau(infile_tau);
    sac->set_filename_corr(infile_corr);
    sac->set_griding_params(grid_interval, spec_interval, omega_min, omega_max);
    sac->set_sampling_params(ndelta, theta, max_annealing_steps, bin_num, bin_size, collecting_steps);
    sac->set_mode_params("fermion", "single");

    sac->init();

    sac->perform_annealing();

    sac->decide_sampling_theta();

    sac->sample_and_collect();

    sac->output(outfile_spec);

    end_t = std::chrono::steady_clock::now();
    std::cout << "total cost "
              << (double)std::chrono::duration_cast<std::chrono::milliseconds>(end_t-begin_t).count()/1000
              << " s." << std::endl;

    delete sac;
    return 0;
}