#include <chrono>
#include <fstream>
#include <iostream>
#include <iomanip>

#include <boost/program_options.hpp>
#include <boost/format.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/constants.hpp>
#include <boost/algorithm/string/classification.hpp>

#include "sac.h"
#include "qmc_data_reader.h"
#include "freq_grids.h"
#include "annealing_chain.h"
#include "measure.h"
#include "random.h"


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
    double theta = 1e6;
    int max_annealing_steps = 5e3;
    int bin_size = 5e3;
    int bin_num = 5;
    int collecting_steps = 1e5;

    std::string kernel_type = "fermion";
    std::string update_type = "single";
    std::string tau_file_path = "../input/benchmark/tau.dat";
    std::string corr_file_path = "../input/benchmark/cor.dat";
    std::string log_file_path = "../output/benchmark/log.log";
    std::string spec_file_path = "../output/benchmark/spec.dat";
    std::string report_file_path = "../output/benchmark/report.dat";

    /** read params from command line */
    boost::program_options::options_description opts("Program options");
    boost::program_options::variables_map vm;

    opts.add_options()
            ("help,h", "display this information")
            ("lt", boost::program_options::value<int>(&lt)->default_value(100),
                    "number of imgainary-time grids of input QMC correlations, default: 100")
            ("beta", boost::program_options::value<double>(&beta)->default_value(4.0),
                    "inverse temperature of QMC, default: 4.0")
            ("nbin-qmc", boost::program_options::value<int>(&nbin)->default_value(1e3),
                    "number of measuring bins in QMC measurements, default: 1e3")
            ("rebin-pace", boost::program_options::value<int>(&rebin_pace)->default_value(1),
                    "pace of rebin, default: 1")
            ("nboostrap", boost::program_options::value<int>(&nboostrap)->default_value(5e3),
                    "number of bootstrap samples for the preprocessing of input QMC data, default: 5e3")
            ("grid-interval", boost::program_options::value<double>(&grid_interval)->default_value(1e-5),
                    "minimum interval of fine frequency grids in the sampling space of delta functions, default: 1e-5")
            ("spec-interval", boost::program_options::value<double>(&spec_interval)->default_value(1e-2),
                    "minimum interval of frequency grids in accumulated spectral function \
                    (should be larger than the intervl of grids), default: 1e-2")
            ("freq-min", boost::program_options::value<double>(&omega_min)->default_value(-10.0),
                    "lower bound of frequency domain, default: -10.0")
            ("freq-max", boost::program_options::value<double>(&omega_max)->default_value(+10.0),
                    "upper bound of frequency domain, default: 10.0")
            ("ndelta", boost::program_options::value<int>(&ndelta)->default_value(1e3),
                    "number of delta functions, default: 1e3")
            ("theta", boost::program_options::value<double>(&theta)->default_value(1e6),
                    "initial sampling temperature, default: 1e6")
            ("max-anneal-steps", boost::program_options::value<int>(&max_annealing_steps)->default_value(5e3),
                    "maximum MC steps for simulated annealing precess, default: 5e3")
            ("max-collecting-steps", boost::program_options::value<int>(&collecting_steps)->default_value(1e5),
                    "maximum MC steps for spectrum collecting precess, default: 1e5")
            ("sbin", boost::program_options::value<int>(&bin_size)->default_value(5e3),
                    "number of bootstrap samples within one bin, default: 5e3")
            ("nbin-sac", boost::program_options::value<int>(&bin_num)->default_value(5),
                    "number of measuring bins for SAC measurements, default: 5")
            ("kernel-type", boost::program_options::value<std::string>(&kernel_type)->default_value("fermion"),
                    "kernel type which relates QMC correlations with spectral function, default: fermion")
            ("update-type", boost::program_options::value<std::string>(&update_type)->default_value("single"),
                    "updating type of MC updates in SAC simulation, default: single")
            ("tau-file-path", boost::program_options::value<std::string>(&tau_file_path)->default_value("../input/benchmark/tau.dat"),
                    "path of input file containing imaginary-time grids from QMC, default: ../input/benchmark/tau.dat")
            ("corr-file-path", boost::program_options::value<std::string>(&corr_file_path)->default_value("../input/benchmark/cor.dat"),
                    "path of input file containing correlation functions measured by QMC, default: ../input/benchmark/cor.dat")
            ("log-file-path", boost::program_options::value<std::string>(&log_file_path)->default_value("../output/benchmark/log.dat"),
                    "output path of logging file during simualtion of SAC, default: ../output/benchmark/log.log")
            ("spec-file-path", boost::program_options::value<std::string>(&spec_file_path)->default_value("../output/benchmark/spec.dat"),
                    "output path of recovered spectral functions, default: ../output/benchmark/spec.dat")
            ("report-file-path", boost::program_options::value<std::string>(&report_file_path)->default_value("../output/benchmark/report.dat"),
                    "output path of the report of recovery quality, default: ../output/benchmark/report.dat");      

    try {
        boost::program_options::store(parse_command_line(argc, argv, opts), vm);
    }
    catch (...) {
        std::cerr << " Got undefined options from command line! \n"<< std::endl;
        exit(1);
    }
    boost::program_options::notify(vm);

    if (vm.count("help")) {
        std::cerr << argv[0] << std::endl;
        std::cerr << opts << std::endl;
        return 0;
    }

    // time record
    std::chrono::steady_clock::time_point begin_t{}, end_t{};
    double duration;

    /** SAC simulations */
    Simulation::SAC *sac = new Simulation::SAC();
    
    // set up random seeds for simulation
    // fixed seed for debug
    // Random::set_seed_fix(12345);
    // Random::set_seed_fix(time(nullptr));

    // print current date and time
    auto current_time = boost::posix_time::second_clock::local_time();
    std::cout << boost::format(" Current time : %s \n") % current_time << std::endl;

    begin_t = std::chrono::steady_clock::now();
    std::cout << " Initialization starts. \n" << std::endl;

    // set up simulating params
    sac->set_qmc_reader_params(lt, beta, nbin, rebin_pace, nboostrap);
    sac->set_file_path_tau(tau_file_path);
    sac->set_file_path_corr(corr_file_path);
    sac->set_outfile_path(log_file_path, spec_file_path, report_file_path);
    sac->set_griding_params(grid_interval, spec_interval, omega_min, omega_max);
    sac->set_sampling_params(ndelta, theta, max_annealing_steps, bin_num, bin_size, collecting_steps);
    sac->set_kernel_type(kernel_type);
    sac->set_update_type(update_type);

    // initialization
    sac->init();

    end_t = std::chrono::steady_clock::now();
    duration = (double)std::chrono::duration_cast<std::chrono::milliseconds>(end_t-begin_t).count()/1000;
    std::cout << boost::format(" Initialization finished in %.2f s. \n") % duration << std::endl;

    begin_t = std::chrono::steady_clock::now();
    boost::format fmt_param_int("%| 35s|%| 5s|%| 18d|");
    boost::format fmt_param_double("%| 35s|%| 5s|%| 18.3e|");
    boost::format fmt_param_range("%| 35s|%| 5s|%| 10.3f|,%| 7.3f|");
    const std::string& joiner = "->";
    std::cout << " Annealing starts with following parameters :\n" << std::endl;
    std::cout << fmt_param_int % "Number of tau points `nt`" % joiner % sac->nt << std::endl;
    std::cout << fmt_param_double % "Sampling temperature `theta`" % joiner % sac->annealing_data->theta << std::endl;
    std::cout << fmt_param_int % "Number of MC sweep `nsweep`" % joiner % sac->measure->size_of_bin << std::endl;
    std::cout << fmt_param_int % "Number of bins `nbin`" % joiner % sac->measure->nbin << std::endl;
    std::cout << fmt_param_range % "Range of spectrum `omega`" % joiner % sac->grids->FreqIndex2Freq(0) % sac->grids->FreqIndex2Freq(sac->grids->FreqNum()-1) << std::endl;
    std::cout << "\n Annealing process ... \n " << std::endl;
    std::cout << boost::format(" Log information of annealing written into %s ... \n") % sac->log_file_path << std::endl;

    // annealing process
    sac->perform_annealing();

    // deciding sampling temperature
    sac->decide_sampling_theta();

    end_t = std::chrono::steady_clock::now();
    duration = (double)std::chrono::duration_cast<std::chrono::milliseconds>(end_t-begin_t).count()/1000;
    auto minute = std::floor(duration/60);
    auto sec = duration - 60 * minute;
    std::cout << boost::format(" Annealing finished in %d min %.2f s. \n") % minute % sec << std::endl;

    begin_t = std::chrono::steady_clock::now();
    std::cout << " Start collecting spectrum ... \n" << std::endl;

    // sampling and collecting
    sac->sample_and_collect();

    // output collected spectrum
    sac->output_recovered_spectrum();

    end_t = std::chrono::steady_clock::now();
    duration = (double)std::chrono::duration_cast<std::chrono::milliseconds>(end_t-begin_t).count()/1000;
    std::cout << boost::format(" Accumulated spectrum stored in %s , costing %.2f s. \n") % sac->spec_file_path % duration << std::endl;

    // output report of recovery quality
    sac->report_recovery_quality();
    std::cout << boost::format(" Quality report of recovered spectrum stored in %s . \n") % sac->report_file_path << std::endl;

    delete sac;
    return 0;
}