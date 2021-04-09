#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>
#include <iostream>
#include <fstream>

#include "sac.h"
#include "measure.h"


/**
 *  TODO: some problems are list here
 *   1. A(\omega) = A(-\omega) ? if so, modify kernel matrix.
 *   2. Displaced Green's functions obtained from QMC seem not to be anti-periodic. ???
 *   3. Update scheme doesn't work if \omega = 0.
 *   4. Regulate theta and nCst to obtain the most reasonable spectrum;
 *      FIXED: during the annealing process, one should keep theta being of order chi^2,
 *             thus configurations can be updated with moderate accepting rate.
 *   5. ...
 */


/** The main program */
int main(int argc, char *argv[]) {

    /* default SAC and Measure parameters */
    int lt = 80;
    double beta = 4.0;
    int nOmega = 50;
    double omegaMin = 0.0, omegaMax = 5.0;

    double theta = exp(20);
    int nCst = 4;

    int nbin = 20;
    int nBetweenBins = (int)pow(10, 3);
    int nstep = 50;
    int nwarm = (int)(pow(10, 5));

    std::string infile_Green = "../results/benchmark_g.txt";
    std::string infile_A = "../results/configs/config_A_" + boost::lexical_cast<std::string>(log(theta) + 0.1) + ".txt";
    std::string outfile_Stats = "../results/stats/stats_A_" + boost::lexical_cast<std::string>(log(theta)) + ".txt";
    std::string outfile_Config =  "../results/configs/config_A_" + boost::lexical_cast<std::string>(log(theta)) + ".txt";

    /* read params from command line */
    boost::program_options::options_description opts("Program options");
    boost::program_options::variables_map vm;

    opts.add_options()
            ("help,h", "display this information")
            ("lt,t", boost::program_options::value<int>(&lt), "imaginary-time lattice size of QMC simulation, default: 80")
            ("beta,b", boost::program_options::value<double>(&beta), "inverse temperature of quantum system, default: 4.0")
            ("nOmega,n", boost::program_options::value<int>(&nOmega), "number of slices in frequency space, default: 50")
            ("omegaMin", boost::program_options::value<double>(&omegaMin), "lower bound of frequency space, default: 0.0")
            ("omegaMax", boost::program_options::value<double>(&omegaMax), "upper bound of frequency space, default: 5.0")
            ("infile,i", boost::program_options::value<std::string>(&infile_Green), "input filename, default: ../results/benchmark_g.txt");
            // TODO

    try {
        boost::program_options::store(parse_command_line(argc, argv, opts), vm);
    }
    catch (...) {
        std::cerr << "Undefined options got from command line."<< std::endl;
        exit(1);
    }

    boost::program_options::notify(vm);

    if (vm.count("help")) {
        std::cerr << argv[0] << std::endl;
        std::cerr << opts << std::endl;
        exit(1);
    }


    /* start SAC and measuring process */
    Measure sacMeasure;

    sacMeasure.set_SAC_params(lt, beta, nOmega, omegaMin, omegaMax);
    sacMeasure.set_meas_params(nbin, nBetweenBins, nstep, nwarm);

    // sacMeasure.set_input_file(infile_Green, infile_A);
    sacMeasure.set_input_file(infile_Green);

    sacMeasure.prepare();

    while (theta > exp(-5)) {
        // annealing process
        sacMeasure.set_sampling_params(theta, nCst);

        sacMeasure.measure();

        sacMeasure.analyse_Stats();

        sacMeasure.print_Stats();

        // sacMeasure.output_Stats(outfile_Stats);

        // sacMeasure.output_Config(outfile_Config);

        theta *= exp(-0.1);
    }

    return 0;
}