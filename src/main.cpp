#include <chrono>
#include <fstream>
#include <iostream>
#include <boost/program_options.hpp>

#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/constants.hpp>

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

    // customized input folder
//    std::string folder_name = "L8b4U-4k0.500.50";
    std::string folder_name = "benchmark";
    std::string in_path = "../input/" + folder_name;
    std::string out_path = "../results/" + folder_name;
    if ( access(out_path.c_str(), 0) != 0 ) {
        std::string command = "mkdir " + out_path;
        if ( system(command.c_str()) != 0 ) {
            std::cerr << "fail to create " + out_path << std::endl;
        }
    }
    // delete previous log
    std::string out_log = out_path + "/log.log";
    if ( access(out_log.c_str(), 0) == 0 ) {
        std::string command = "rm " + out_log;
        int status = system(command.c_str());
    }

    std::cout << "Initialization starts ..." << std::endl;

    sac->set_read_in_params(lt, beta, nbin, rebin_pace, nboostrap);
    sac->set_filename_tau(in_path + "/tau.dat");
    sac->set_filename_corr(in_path + "/cor.dat");
    sac->set_filename_log(out_log);
    sac->set_griding_params(grid_interval, spec_interval, omega_min, omega_max);
    sac->set_sampling_params(ndelta, theta, max_annealing_steps, bin_num, bin_size, collecting_steps);
    sac->set_mode_params("fermion", "single");

    sac->init();

    std::cout << "Initialization finished." << std::endl
              << "Annealing starts with params :" << std::endl
              << "  nt     = " << sac->nt << std::endl
              << "  theta  = " << sac->data->theta << std::endl
              << "  nsweep = " << sac->measure->sbin << std::endl
              << "  nbin   = " << sac->measure->nbin << std::endl
              << "  omega  = " << sac->grid->GridIndex2Freq(0) << ", "
                               << sac->grid->GridIndex2Freq(sac->grid->GridsNum()-1) << std::endl
              << "......" << std::endl;

    sac->perform_annealing();

    sac->decide_sampling_theta();

    std::cout << "Annealing finished." << std::endl
              << "Start collecting spectrum ... " << std::endl;

    sac->sample_and_collect();

    sac->output(out_path + "/spec.dat");

    std::cout << "Accumulated spectrum has been writen into " + out_path + "/spec.dat ." << std::endl;

    end_t = std::chrono::steady_clock::now();
    std::cout << "Total cost: "
              << (double)std::chrono::duration_cast<std::chrono::milliseconds>(end_t-begin_t).count()/1000
              << " s." << std::endl;

    /** test **/
//    std::vector<double> freq_vec;
//    std::vector<double> spec_vec;
//
//    std::ifstream infile;
//    infile.open("../input/benchmark/spec_base.dat", std::ios::in);
//    if (!infile.is_open()) {
//        exit(1);
//    }
//    std::string line;
//    std::vector<std::string> data;
//    getline(infile, line);
//    boost::split(data, line, boost::is_any_of(" "), boost::token_compress_on);
//    data.erase(std::remove(std::begin(data), std::end(data), ""), std::end(data));
//    while(getline(infile, line)) {
//        boost::split(data, line, boost::is_any_of(" "), boost::token_compress_on);
//        data.erase(std::remove(std::begin(data), std::end(data), ""), std::end(data));
//        freq_vec.push_back(boost::lexical_cast<double>(data[0]));
//        spec_vec.push_back(boost::lexical_cast<double>(data[1]));
//    }
//    infile.close();
//
//
//    Eigen::VectorXd freq(freq_vec.size()), spec(freq_vec.size());
//    for (int f = 0; f < freq_vec.size(); ++f) {
//        freq(f) = freq_vec[f];
//        spec(f) = spec_vec[f];
//    }
//
//    Eigen::MatrixXd kernel(sac->nt, freq.size());
//    for (int t = 0; t < sac->nt; ++t) {
//        for (int f = 0; f < freq.size(); ++f) {
//            kernel(t, f) = exp(-freq[f]*sac->tau[t]) / (2*M_PI*(1+exp(-sac->beta*freq[f])));
//        }
//    }
//    kernel = sac->readin->rotate_mat * kernel;
//    spec = spec/sac->scale_factor;
//
//    std::cout << ((spec_interval*kernel*spec-sac->corr).array() * sac->sigma.array()).square().sum() << std::endl;
//
//    Eigen::VectorXd err = spec_interval*kernel*spec - sac->corr;
//    for (int i = 0; i < sac->readin->cov_mat_dim; ++i) {
//        std::cout << sac->readin->corr_mean[i] << "     "
//                  << sac->readin->corr_err[i] << "     "
//                  << sac->readin->cov_eig[i] << "     "
//                  << err[i] << "     "
//                  << sac->sigma[i] << std::endl;
//    }


    delete sac;
    return 0;
}