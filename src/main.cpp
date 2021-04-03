#include <iostream>
#include "stochasticAC.h"


int main() {

    StochasticAC sac;

    sac.set_SAC_params(80, 4.0, 50, 0.0, 5.0);

    sac.initialSAC();

    sac.read_QMC_data("../results/benchmark_g.txt");

    std::cout.precision(16);
    std::cout << sac.tau_list(75) << std::endl;
    std::cout << sac.green_tau_QMC(75) << std::endl;
    std::cout << sac.err_tau_QMC(75) << std::endl;

    return 0;
}