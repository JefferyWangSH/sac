#include <iostream>
#include <fstream>
#include <iomanip>

#include <boost/program_options.hpp>

#include "SAC.h"
#include "ReadInModule.h"

#include <chrono>


/**
 *  TODO:
 *   1. Optimize random module (missing)
 *   2. ...
 */



/** The main program */
int main(int argc, char *argv[]) {

    int lt = 101;
    double beta = 4.0;
    int nbin = 1000;
    int rebin_pace = 1;
    int num_boostrap = 5000;

//    std::chrono::steady_clock::time_point begin_t{}, end_t{};
//    begin_t = std::chrono::steady_clock::now();

    ReadInModule read_in;

    read_in.set_params(lt, beta, nbin, rebin_pace, num_boostrap);
    read_in.read_tau_from_file("../input/benchmark/tgrid.dat");
    read_in.read_corr_from_file("../input/benchmark/cor.dat");
    read_in.compute_corr_means();
    read_in.compute_cov_matrix();
//    read_in.deallocate_memory();

//    std::cout << read_in.lt << std::endl;
//    std::cout << read_in.tau_seq << std::endl;
//
//    std::cout << read_in.nbin << std::endl;
//    std::cout << read_in.nbin_total << std::endl;
//
//    for (int l = 0; l < read_in.lt; ++l) {
//        std::cout << read_in.corr_tau_bin[0][l] << std::endl;
//    }
//
//    std::cout << read_in.tau_seq << std::endl << std::endl;
//    std::cout << read_in.corr_mean_seq << std::endl << std::endl;
//    std::cout << read_in.corr_err_seq << std::endl;

    for (auto i = 0; i < read_in.cov_eig.size(); ++i) {
        read_in.cov_eig[i] = pow(read_in.cov_eig[i]/num_boostrap, 0.5);
    }

    std::cout << read_in.cov_eig << std::endl;
    std::cout << read_in.cov_eig.size() << std::endl;
//    std::cout << read_in.rotate_mat << std::endl;

//    std::cout << (read_in.cov_mat - read_in.cov_mat.transpose()).maxCoeff() << std::endl;
//    std::cout << ( read_in.rotate_mat.transpose() * read_in.cov_eig.asDiagonal() * read_in.rotate_mat - read_in.cov_mat ).maxCoeff() << std::endl;

//    end_t = std::chrono::steady_clock::now();
//    std::cout << (double)std::chrono::duration_cast<std::chrono::milliseconds>(end_t - begin_t).count()/1000 << std::endl;



    return 0;
}