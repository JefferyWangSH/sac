#include <boost/program_options.hpp>
#include <iostream>
#include <fstream>
#include <iomanip>

#include "SAC.h"
#include "MonteCarloSAC.h"


/**
 *  TODO:
 *   1. Support reading configurations from file
 *   2. ...
 */


/** The main program */
int main(int argc, char *argv[]) {

    /*
    SAC sac;

    sac.set_SAC_params(80, 4, 49, -5, 5, 1);

    sac.set_alpha(exp(-10));

    sac.set_QMC_filename("../results/data/gt_l4_lt80_u-4.0_b4.0_k_pi2pi2.txt");

    sac.prepare();

    for (int n = 0; n < 50; ++n) {
        sac.Metropolis_update_1step();
    }

    double sum = 0.0;
    for (int i = 0; i < sac.n_config; ++i) {
        std::cout << sac.x_list[i] << "   " << sac.n_list[i] << std::endl;
        sum += sac.n_list[i];
    }
    std::cout << std::endl;
    std::cout << sum << std::endl;
    */

    MonteCarloSAC sac;

    sac.set_SAC_params(80, 4, 49, -5, 5, 1);

    sac.set_measure_params(100, 50, 10, 5e5, 100);

    sac.set_input_file("../results/data/gt_l4_lt80_u-4.0_b4.0_k_pi2pi2.txt");

    std::vector<double> alpha_list(60);
    for (int i = 0; i < 60; ++i) {
        alpha_list[i] = ( i == 0 )? exp(-6) : alpha_list[i-1] * exp(0.1);
    }
    sac.set_tempering_profile(60, alpha_list);

    sac.prepare();

    std::chrono::steady_clock::time_point begin_t = std::chrono::steady_clock::now();
    std::chrono::steady_clock::time_point end_t;

    sac.run_Monte_Carlo();


    std::cout << std::endl;
    for (auto p : sac.p_list[10]) {
        std::cout << p << std::endl;
    }
    std::cout << std::endl;

    end_t = std::chrono::steady_clock::now();
    auto time = std::chrono::duration_cast<std::chrono::milliseconds>(end_t - begin_t).count();
    const int minute = std::floor((double)time / 1000 / 60);
    const double sec = (double)time / 1000 - 60 * minute;

    std::cout << minute << " min " << sec << " s" << std::endl;


    std::ofstream outfile;
    outfile.open("../results/test.txt", std::ios::out | std::ios::trunc);

    std::vector<double> H_list(sac.nalpha);
    for (int i = 0; i < sac.nalpha; ++i) {
        sac.sac_list[i].cal_Config_Hamiltonian(H_list[i]);
        outfile << std::setiosflags(std::ios::right)
                << std::setw(15) << i
                << std::setw(15) << log(H_list[i])
                << std::endl;
    }
    outfile.close();

    return 0;
}