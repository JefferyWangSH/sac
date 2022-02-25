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


namespace Debug {

    /**
      * Subroutine to calculate fitting errors of correlations function after a SAC simulation was finished.
      * The differences between recovered correlations and QMC measured correlations are calculated,
      * which to some extent indicates the correctness of our accumulated spectrum.
      *
      * @param sac - SAC object which performed the SAC simulations
      * @param spec_file - name of file which includes the recovered spectrum
      * @retval None
      */
    void calculate_fitting_error(const Simulation::SAC &sac, const std::string &spec_file) {
        // calculate fitting error of correlations
        // first read spectrum information from file
        std::vector<double> freq_vec;
        std::vector<double> spec_vec;

        std::ifstream infile;
        infile.open(spec_file, std::ios::in);
        if (!infile.is_open()) {
            std::cerr << boost::format(" Fail to open file %s ! \n") % spec_file << std::endl;
            exit(1);
        }
        std::string line;
        std::vector<std::string> data;
        getline(infile, line);
        boost::split(data, line, boost::is_any_of(" "), boost::token_compress_on);
        data.erase(std::remove(std::begin(data), std::end(data), ""), std::end(data));
        while(getline(infile, line)) {
            boost::split(data, line, boost::is_any_of(" "), boost::token_compress_on);
            data.erase(std::remove(std::begin(data), std::end(data), ""), std::end(data));
            freq_vec.push_back(boost::lexical_cast<double>(data[0]));
            spec_vec.push_back(boost::lexical_cast<double>(data[1]));
        }
        infile.close();

        // convert to eigen type
        Eigen::VectorXd freq = Eigen::Map<Eigen::VectorXd>(freq_vec.data(), freq_vec.size());
        Eigen::VectorXd spec = Eigen::Map<Eigen::VectorXd>(spec_vec.data(), spec_vec.size());

        // generate kernel
        Eigen::MatrixXd kernel(sac.nt, freq.size());
        for (int t = 0; t < sac.nt; ++t) {
            for (int f = 0; f < freq.size(); ++f) {
                kernel(t, f) = exp(-freq[f]*sac.tau_from_qmc[t]) / (2*M_PI*(1+exp(-sac.beta*freq[f])));
            }
        }
        // rotate and refactor
        kernel = sac.qmc_data_reader->rotate_mat*kernel;
        spec = spec/sac.scale_factor;

        Eigen::VectorXd err = sac.grids->SpecInterval()*kernel*spec - sac.corr_from_qmc;
        for (int i = 0; i < sac.qmc_data_reader->cov_mat_dim; ++i) {
            std::cout << std::setiosflags(std::ios::right)
                      << std::setw(15) << sac.qmc_data_reader->corr_mean_qmc[i]
                      << std::setw(15) << sac.qmc_data_reader->corr_err_qmc[i]
                      << std::setw(15) << sac.qmc_data_reader->cov_eig[i]
                      << std::setw(15) << sac.corr_from_qmc[i]
                      << std::setw(15) << err[i]
                      << std::setw(15) << sac.sigma_from_qmc[i] << std::endl;
        }
    }

}

namespace Customize {
    
    /**
      * Subroutine for user customized input.
      *
      * @param  input_folder_name: folder name which provides input data
      * @param  input_file_tau_path: path of input tau data
      * @param  input_file_cor_path: path of input cor data
      * @param  output_file_log_path: path of output log
      * @param  output_file_spec_path: path of output specturm
      * @retval None
      */
    void customize_input_folder(const std::string& input_folder_name, 
                                std::string& input_file_tau_path, std::string& input_file_cor_path,
                                std::string& output_file_log_path, std::string& output_file_spec_path) 
    {
        // folder names
        const std::string& input_folder_path = "../input/" + input_folder_name;
        const std::string& output_folder_path = "../output/" + input_folder_name;
        input_file_tau_path = input_folder_path + "/tau.dat";
        input_file_cor_path = input_folder_path + "/cor.dat";
        output_file_spec_path = output_folder_path + "/spec.dat";
        output_file_log_path = output_folder_path + "/log.log";

        // confirm existence of input data
        if ( access(input_folder_path.c_str(), 0) != 0 ) {
            std::cerr << boost::format(" Input folder %s not exists. \n") % input_folder_path << std::endl;
            exit(1);
        }
        else {
            if ( access(input_file_tau_path.c_str(), 0) != 0 ) {
                std::cerr << boost::format(" Cannot open input file %s . \n") % input_file_tau_path << std::endl;
                exit(1);
            }
            else if ( access(input_file_cor_path.c_str(), 0) != 0 ) {
                std::cerr << boost::format(" Cannot open input file %s . \n") % input_file_cor_path << std::endl;
                exit(1);
            }
            else {
                std::cout << boost::format(" Input files detected in folder %s . \n") % input_folder_path << std::endl;
            }
        }
        
        // try to create corresponding output folder
        if ( access(output_folder_path.c_str(), 0) != 0 ) {
            const std::string& command = "mkdir " + output_folder_path;
            if ( system(command.c_str()) != 0 ) {
                std::cerr << boost::format(" Fail to create folder %s \n") % output_folder_path << std::endl;
                exit(1);
            }
        }

        // delete previous log file if exist.
        if ( access(output_file_log_path.c_str(), 0) == 0 ) {
            const std::string& command = "rm " + output_file_log_path;
            if ( system(command.c_str()) != 0 ) {
                std::cerr << boost::format(" Fail to delete folder %s \n") % output_file_log_path << std::endl;
                exit(1);
            }
        }
    }
}


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

    std::string input_folder_name = "benchmark";

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
            ("nboostrap", boost::program_options::value<int>(&nboostrap)->default_value(5e3),
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
            ("input-folder-name", boost::program_options::value<std::string>(&input_folder_name)->default_value("benchmark"),
                    "name of input folder containing QMC data, default: benchmark");      

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

    // customized input folder
    std::string infile_tau_path, infile_cor_path;
    std::string outfile_log_path, outfile_spec_path;
    Customize::customize_input_folder(input_folder_name, infile_tau_path, infile_cor_path, outfile_log_path, outfile_spec_path);

    begin_t = std::chrono::steady_clock::now();
    std::cout << " Initialization starts . \n" << std::endl;

    // set up simulating params
    sac->set_read_in_params(lt, beta, nbin, rebin_pace, nboostrap);
    sac->set_filename_tau(infile_tau_path);
    sac->set_filename_corr(infile_cor_path);
    sac->set_filename_log(outfile_log_path);
    sac->set_griding_params(grid_interval, spec_interval, omega_min, omega_max);
    sac->set_sampling_params(ndelta, theta, max_annealing_steps, bin_num, bin_size, collecting_steps);
    sac->set_mode_params("fermion", "single");

    // initialization
    sac->init();

    // // output initialization results
    // for (int i = 0; i < sac->qmc_data_reader->cov_mat_dim; ++i) {
    //     std::cout << std::setiosflags(std::ios::right)
    //               << std::setw(15) << sac->qmc_data_reader->corr_mean[i]
    //               << std::setw(15) << sac->qmc_data_reader->corr_err[i]
    //               << std::setw(15) << sac->qmc_data_reader->cov_eig[i]
    //               << std::setw(15) << sac->corr[i]
    //               << std::setw(15) << sac->sigma[i] << std::endl;
    // }

    end_t = std::chrono::steady_clock::now();
    duration = (double)std::chrono::duration_cast<std::chrono::milliseconds>(end_t-begin_t).count()/1000;
    std::cout << boost::format(" Initialization finished in %.2f s. \n") % duration << std::endl;

    begin_t = std::chrono::steady_clock::now();
    boost::format fmt_param_int("%| 35s|%| 5s|%| 10d|");
    boost::format fmt_param_double("%| 35s|%| 5s|%| 10.2e|");
    boost::format fmt_param_range("%| 35s|%| 5s|%| 10.2f|, %.2f");
    const std::string& joiner = "->";
    std::cout << " Annealing starts with following parameters :" << std::endl;
    std::cout << fmt_param_int % "Number of tau points `nt`" % joiner % sac->nt << std::endl;
    std::cout << fmt_param_double % "Sampling temperature `theta`" % joiner % sac->annealing_data->theta << std::endl;
    std::cout << fmt_param_int % "Number of MC sweep `nsweep`" % joiner % sac->measure->size_of_bin << std::endl;
    std::cout << fmt_param_int % "Number of bins `nbin`" % joiner % sac->measure->nbin << std::endl;
    std::cout << fmt_param_range % "Range of spectrum `omega`" % joiner % sac->grids->FreqIndex2Freq(0) % sac->grids->FreqIndex2Freq(sac->grids->FreqNum()-1) << std::endl;
    std::cout << " Annealing process ... \n " << std::endl;

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
    sac->output(outfile_spec_path);

    end_t = std::chrono::steady_clock::now();
    duration = (double)std::chrono::duration_cast<std::chrono::milliseconds>(end_t-begin_t).count()/1000;
    std::cout << boost::format(" Accumulated spectrum has been stored into %s , costing %.2f s. \n") % outfile_spec_path % duration << std::endl;

    // output the fitting error
    Debug::calculate_fitting_error(*sac, outfile_spec_path);

    delete sac;
    return 0;
}