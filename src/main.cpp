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
 *      FIXED: during the annealing process, one should keep theta being of order hi^2,
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

    double theta = exp(25);
    int nCst = 4;

    int nbin = 20;
    int nBetweenBins = (int)pow(10, 3);
    int nstep = 50;
    int nwarm = (int)pow(10, 5);

    std::string infilename = "../results/benchmark_g.txt";
    std::string outfilename = "../results/output.txt";

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
            ("infile,i", boost::program_options::value<std::string>(&infilename), "input filename, default: ../results/benchmark_g.txt")
            ("outfile,o", boost::program_options::value<std::string>(&outfilename), "output filename, default: ../results/output.txt");

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
    sacMeasure.set_input_file(infilename);

    sacMeasure.prepare();

    sacMeasure.set_sampling_params(theta, nCst);

    sacMeasure.measure();

    sacMeasure.analyse_Stats();

    sacMeasure.print_Stats();

    sacMeasure.output_Stats(outfilename);

    std::cout << sacMeasure.sac.A_omega << std::endl;


    // annealing process
    /*
    std::ofstream outfile;
    outfile.open(outfilename, std::ios::out | std::ios::trunc);

    while (theta >= exp(0)) {
        sacMeasure.set_sampling_params(theta, nCst);

        sacMeasure.measure();

        sacMeasure.analyse_Stats();

        sacMeasure.print_Stats();

        sacMeasure.output_Stats(outfilename);

        theta /= exp(0.5);
    }

    std::cout << sacMeasure.sac.A_omega << std::endl;
    std::cout << sacMeasure.sac.A_omega.sum() * sacMeasure.sac.deltaOmega << std::endl;
    */

    return 0;
}