#include <iostream>
#include <fstream>
#include <iomanip>

#include <boost/program_options.hpp>

#include "SAC.h"
#include "MonteCarloSAC.h"


/**
 *  TODO:
 *   1. Support reading configurations from file
 *   2. update in parallel
 *   3. ...
 */


/** The main program */
int main(int argc, char *argv[]) {

    int lt = 80;
    double beta = 4;

    double omega_min = -8.0;
    double omega_max =  8.0;

    int nMoment = 1;
    int nconfig = 100;

    int nbin = 2000;
    int nstep_1bin = 100;
    int step_between_bins = 10;
    int nwarm = 3e5;

    /** Recover averaged spectrum */
    /*
    MonteCarloSAC sac;

    sac.set_SAC_params(lt, beta, nconfig, omega_min, omega_max, nMoment);
    sac.set_QMC_filename("../results/data/gt_l4_lt80_u-4.0_b4.0_k_pi2pi2.txt");

    sac.set_Configs_filename("../results/configs/alpha_2.50.txt");

    double alpha = exp(-5);
    sac.set_alpha(alpha);

    sac.prepare();

    std::vector<double> gt(lt);
    for (int t = 0; t < lt; ++t) {
        gt[t] = 0.0;
        for (int i = 0; i < nconfig; ++i) {
            gt[t] += sac.kernel[t][i] * sac.n_list[i];
        }
    }
    for (auto g : gt) {
        std::cout << g << std::endl;
    }
    double H = 0.0;
    sac.cal_Config_Hamiltonian(H);
    std::cout << H << std::endl;
    */

    /** Monte Carlo procedure */

    MonteCarloSAC sac;

    sac.set_SAC_params(lt, beta, nconfig, omega_min, omega_max, nMoment);

    sac.set_measure_params(nbin, nstep_1bin, step_between_bins, nwarm);

    sac.set_QMC_filename("../results/data/gt_l4_lt80_u-4.0_b4.0_k_pi2pi2.txt");

    double alpha = exp(-29.9);
    double dalpha = exp(0.1);

    for (int n = 0; n < 400; ++n) {

        sac.set_alpha(alpha);

        std::stringstream ss;
        ss << std::setiosflags(std::ios::fixed) << std::setprecision(2) << log(alpha) - log(dalpha);
        std::string lnAlpha_old = ss.str();
        ss.str("");
        std::string infile_configs = "../results/configs/alpha_" + lnAlpha_old + ".txt";
        sac.set_Configs_filename(infile_configs);

        sac.prepare();

        sac.run_Monte_Carlo();

        sac.print_stats();

        ss << std::setiosflags(std::ios::fixed) << std::setprecision(2) << log(alpha);
        std::string lnAlpha = ss.str();
        ss.str("");
        std::string outfile_configs = "../results/configs/alpha_" + lnAlpha + ".txt";
        sac.output_Configs(outfile_configs);

        std::string outfile_hamilton = "../results/h-alpha.txt";
        sac.output_stats(outfile_hamilton);

        alpha *= dalpha;
    }

    return 0;
}