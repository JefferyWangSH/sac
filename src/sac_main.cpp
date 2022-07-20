#include "sac_core.h"
#include "sac_kernel.h"
#include "sac_annealing.h"
#include "sac_measure.h"
#include "sac_writer.h"
#include "freq_grids.h"
#include "qmc_reader.h"
#include "random.h"
#include "utils/toml.hpp"

#include <chrono>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>

#include <boost/program_options.hpp>
#include <boost/format.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>


int main( int argc, char *argv[] ) {


    // -------------------------------------------------------------------------------------------------------------
    //                              Read Input and Output Files from Prigram Options
    // -------------------------------------------------------------------------------------------------------------

    // input and output files, read from program options
    std::string config_file = "../benchmark/config.toml";
    std::string tgrids_file = "../benchmark/tgrids.in";
    std::string corr_file = "../benchmark/corr.in";
    std::string log_file = "../benchmark/log.out";
    std::string spec_file = "../benchmark/spec.out";
    std::string report_file = "../benchmark/report.out";

    // read params from command line
    boost::program_options::options_description opts( "Program options" );
    boost::program_options::variables_map vm;

    opts.add_options()
        ( "help,h", "display this information" )

        ( "config", boost::program_options::value<std::string>( &config_file )->default_value( "../benchmark/config.toml" ),
                "toml configuration file for SAC, default: ../benchmark/config.toml" )
        ( "tgrids", boost::program_options::value<std::string>( &tgrids_file )->default_value( "../benchmark/tgrids.in" ),
                "input file which contains the QMC imaginary-time grids, default: ../benchmark/tgrids.in" )
        ( "corr", boost::program_options::value<std::string>( &corr_file )->default_value( "../benchmark/corr.in" ),
                "input file which contains the QMC correlation functions, default: ../benchmark/corr.in" )
        ( "log", boost::program_options::value<std::string>( &log_file )->default_value( "../benchmark/log.out" ),
                "log file which records the history of the simulated annealing in SAC, default: ../benchmark/log.out" )
        ( "spec", boost::program_options::value<std::string>( &spec_file )->default_value( "../benchmark/spec.out" ),
                "output file which contains the recovered spectral functions, default: ../benchmark/spec.out" )
        ( "report", boost::program_options::value<std::string>( &report_file )->default_value( "../benchmark/report.out" ),
                "output file which contains the quality report of SAC, default: ../benchmark/report.out" );

    try {
        boost::program_options::store( parse_command_line(argc, argv, opts), vm );
    }
    catch ( ... ) {
        std::cerr << "main(): undefined program options got from the command line." << std::endl;
        exit(1);
    }
    boost::program_options::notify( vm );

    // show the helping message
    if ( vm.count( "help" )) {
        std::cerr << argv[0] << std::endl;
        std::cerr << opts << std::endl;
        return 0;
    }
    

    // ------------------------------------------------------------------------------------------------------------
    //                           Helping Variables and Functions to Organize the Program
    // ------------------------------------------------------------------------------------------------------------

    // parse the configuration file
    auto config = toml::parse_file( config_file );

    // time recorders
    std::chrono::steady_clock::time_point begin_time{}, end_time{};
    double seconds{};

    // lambda function for the standard time output
    auto StandardTimeOutput = [](double seconds) {
        auto day = std::floor(seconds/86400);
        auto hour = std::floor((seconds - day * 86400) / 3600);
        auto minute = std::floor((seconds - day * 86400 - hour * 3600) / 60);
        auto sec = seconds - 86400 * day - 3600 * hour - 60 * minute;
        if ( day ) { return boost::format("%d d %d h %d m %.2f s") % day % hour % minute % sec; }
        else if ( hour ) { return boost::format("%d h %d m %.2f s") % hour % minute % sec; }
        else if ( minute ) { return boost::format("%d m %.2f s") % minute % sec; }
        else { return boost::format("%.2f s") % sec; }
    };

    // // set up random seeds for the Monte Carlo simulation
    // // fixed seed for debug
    // Utils::Random::set_seed( 12345 );
    // Utils::Random::set_seed( time(nullptr) );


    // ------------------------------------------------------------------------------------------------------------
    //                        Initialize SAC Modules and Preprocess the QMC Correlations
    // ------------------------------------------------------------------------------------------------------------

    // print current date and time
    auto current_time = boost::posix_time::second_clock::local_time();
    std::cout << boost::format(" Current time : %s \n") % current_time << std::endl;

    begin_time = std::chrono::steady_clock::now();
    std::cout << " Initialization starts. ( Preprocessing of input QMC correlations ) \n" << std::endl;


    // ---------------------------------  Initialize the QmcReader module  ----------------------------------------

    const double beta = config["QmcReader"]["beta"].value_or(8.0);
    const int num_tgrids_qmc = config["QmcReader"]["number_of_tgrids"].value_or(160);
    const int num_bin_qmc = config["QmcReader"]["number_of_bins"].value_or(1e3);
    const int num_bootstrap = config["QmcReader"]["number_of_bootstraps"].value_or(5e3);
    const int rebin_pace = config["QmcReader"]["pace_of_rebin"].value_or(1);
    
    SAC::Initializer::QmcReader* qmc_reader = new SAC::Initializer::QmcReader();
    qmc_reader->set_params( num_tgrids_qmc, beta, num_bin_qmc, rebin_pace, num_bootstrap );
    qmc_reader->read_tgrids_from_file( tgrids_file );
    qmc_reader->read_corr_from_file( corr_file );
    
    // transform the data into the diagonal space of the QMC orrelation functions
    qmc_reader->analyse_corr();
    qmc_reader->filter_and_rotate();


    // ----------------------------------  Initialize the FreqGrids module  ---------------------------------------

    const double freq_interval = config["SAC"]["FreqGrids"]["freq_interval"].value_or(1e-5);
    const double spec_interval = config["SAC"]["FreqGrids"]["spec_interval"].value_or(1e-2);
    const double freq_min = config["SAC"]["FreqGrids"]["freq_min"].value_or(-10.0);
    const double freq_max = config["SAC"]["FreqGrids"]["freq_max"].value_or(+10.0);
    Grids::FreqGrids* grids = new Grids::FreqGrids();
    grids->set_grids_params( freq_interval, spec_interval, freq_min, freq_max );
    grids->initial();


    // -----------------------------------  Initialize the Kernel module  -----------------------------------------

    const std::string kernel_type = config["SAC"]["Types"]["kernel_type"].value_or("fermion");
    SAC::Kernel* kernel = new SAC::Kernel();
    kernel->set_kernel_params( qmc_reader->cov_mat_dim(), grids->FreqNum(), kernel_type );
    kernel->initial( *qmc_reader, *grids );


    // -----------------------------------  Initialize the Measure module  ----------------------------------------

    const int num_bin_sac = config["SAC"]["Measure"]["number_of_bins"].value_or(5);
    const int size_bin_sac = config["SAC"]["Measure"]["size_of_bins"].value_or(1e3);
    SAC::Measure* measure = new SAC::Measure();
    measure->resize( num_bin_sac, size_bin_sac );


    // -------------------------------  Initialize the Annealing::Chain module  -----------------------------------

    const int max_annealing_steps = config["SAC"]["Annealing"]["max_annealing_steps"].value_or(5e3);
    SAC::Annealing::Chain* chain = new SAC::Annealing::Chain();
    chain->set_max_length( max_annealing_steps );


    // -----------------------------------  Initialize the SacCore module  ----------------------------------------

    const int num_deltas = config["SAC"]["Annealing"]["number_of_deltas"].value_or(1e3);
    const int collecting_steps = config["SAC"]["Measure"]["collecting_steps"].value_or(1e5);
    const int stabilization_pace = config["SAC"]["Annealing"]["pace_of_stabilization"].value_or(10);
    const double theta = config["SAC"]["Annealing"]["theta"].value_or(1e8);
    const double annealing_rate = config["SAC"]["Annealing"]["annealing_rate"].value_or(0.9);
    const std::string update_type = config["SAC"]["Types"]["update_type"].value_or("single");
    
    SAC::SacCore* core = new SAC::SacCore();
    core->set_sampling_params( num_deltas, collecting_steps, stabilization_pace, update_type );
    core->set_annealing_params( theta, annealing_rate, log_file );
    core->initial( *kernel, *qmc_reader, *grids );


    // ----------------------------------  Output the Simulation Parameters  --------------------------------------

    end_time = std::chrono::steady_clock::now();
    seconds = (double)std::chrono::duration_cast<std::chrono::milliseconds>( end_time - begin_time ).count() / 1000;
    std::cout << boost::format(" Initialization finished in %s. \n") % StandardTimeOutput(seconds) << std::endl;

    // output formats
    boost::format fmt_str("%| 42s|%| 5s|%| 18s|");
    boost::format fmt_int("%| 42s|%| 5s|%| 18d|");
    boost::format fmt_double("%| 42s|%| 5s|%| 18.2f|");
    boost::format fmt_science("%| 42s|%| 5s|%| 18.2e|");
    boost::format fmt_range("%| 42s|%| 5s|%| 10.2f|,%| 7.2f|");
    const std::string& joiner = "->";
    std::cout << " Annealing starts with the following parameters :\n" << std::endl;

    std::cout << fmt_int % "Number of tau points in QMC 'lt'" % joiner % qmc_reader->time_num() << std::endl;
    std::cout << fmt_double % "Inverse temperature in QMC 'beta'" % joiner % qmc_reader->beta() << std::endl;
    std::cout << fmt_int % "Number of bins in QMC 'nbin-qmc'" % joiner % qmc_reader->bin_num_total() << std::endl;
    std::cout << fmt_int % "Rebin pace of QMC bin data 'rebin_pace'" % joiner % qmc_reader->rebin_pace() << std::endl;
    std::cout << fmt_int % "Number of bootstrap samples 'nbootstrap'" % joiner % qmc_reader->bootstrap_num() << std::endl << std::endl;
    
    std::cout << fmt_str % "Type of kernel 'kernel_type'" % joiner % kernel_type << std::endl;
    std::cout << fmt_str % "Type of MC updates 'update_type'" % joiner % update_type << std::endl << std::endl;
    
    std::cout << fmt_int % "Number of tau points (truncated) 'nt'" % joiner % core->TimeSize() << std::endl;
    std::cout << fmt_science % "Initial sampling temperature 'theta'" % joiner % core->Theta() << std::endl;
    std::cout << fmt_double % "Annealing rate 'annealing_pace'" % joiner % core->AnnealingRate() << std::endl;
    std::cout << fmt_double % "Stablization rate 'stablization_pace'" % joiner % core->StabilizationPace() << std::endl;
    std::cout << fmt_int % "Number of delta functions 'ndelta'" % joiner % core->NumDeltas() << std::endl;
    std::cout << fmt_int % "Number of meausring bins 'nbin-sac'" % joiner % measure->number_of_bin() << std::endl;
    std::cout << fmt_int % "Capacity of a measuring bin 'sbin-sac'" % joiner % measure->size_of_bin() << std::endl;
    std::cout << fmt_int % "Number of spec samples 'collecting_steps'" % joiner % core->CollectingSteps() << std::endl << std::endl;
    
    std::cout << fmt_science % "Interval of fine grids 'freq_interval'" % joiner % grids->FreqInterval() << std::endl;
    std::cout << fmt_science % "Interval of spectrum 'spec_interval'" % joiner % grids->SpecInterval() << std::endl;
    std::cout << fmt_range % "Frequency range of spectrum 'freq'" % joiner 
                           % grids->FreqIndex2Freq(0) % grids->FreqIndex2Freq( grids->FreqNum()-1 ) << std::endl;


    // ------------------------------------------------------------------------------------------------------------
    //                                         Perform SAC Simulations
    // ------------------------------------------------------------------------------------------------------------
    
    std::cout << "\n Annealing process ... \n" << std::endl;
    std::cout << boost::format(" Log information of annealing written into %s ... \n") % log_file << std::endl;
    begin_time = std::chrono::steady_clock::now();

    // annealing process
    core->perform_annealing( *kernel, *grids, *measure, *chain );

    // decide the sampling temperature
    core->decide_sampling_theta( *kernel, *chain );

    end_time = std::chrono::steady_clock::now();
    seconds = (double)std::chrono::duration_cast<std::chrono::milliseconds>( end_time - begin_time ).count()/1000;
    std::cout << boost::format(" Annealing finished in %s. \n") % StandardTimeOutput(seconds) << std::endl;

    begin_time = std::chrono::steady_clock::now();
    std::cout << " Start collecting spectrum ... \n" << std::endl;

    // sampling and collecting
    core->sample_and_collect( *kernel, *grids, *measure, *chain );

    // output collected spectrum
    SAC::Writer::write_spectrum( spec_file, *core );

    end_time = std::chrono::steady_clock::now();
    seconds = (double)std::chrono::duration_cast<std::chrono::milliseconds>( end_time - begin_time ).count()/1000;
    std::cout << boost::format(" Accumulated spectrum stored in %s , costing %s. \n") % spec_file % StandardTimeOutput(seconds) << std::endl;

    // output report of recovery quality
    SAC::Writer::write_quality_report( report_file, *core, *kernel, *grids, *qmc_reader );
    std::cout << boost::format(" Quality report of recovered spectrum stored in %s . \n") % report_file << std::endl;


    // memory release
    delete qmc_reader;
    delete core;
    delete kernel;
    delete measure;
    delete chain;
    delete grids;

}
