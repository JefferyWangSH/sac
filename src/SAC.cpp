#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

#include "SAC.h"

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>


void SAC::set_SAC_params(int lt, double beta, int nOmega, double omegaMin, double omegaMax) {
    /*
     *  Set SAC parameters.
     */

    // we restrict ourselves to the subspace where \omega > 0
    assert(lt > 0 && beta > 0 && nOmega > 0);
    assert(omegaMin >= 0);
    assert(omegaMin < omegaMax);

    this->lt = lt;
    this->beta = beta;
    this->dtau = beta/lt;
    this->nOmega = nOmega;
    this->omegaMin = omegaMin;
    this->omegaMax = omegaMax;
    this->deltaOmega = (omegaMax - omegaMin) / nOmega;

    omega_list.resize(nOmega);
    tau_list.resize(lt);
    A_omega.resize(nOmega);
    g_tau.resize(lt);
    g_tau_QMC.resize(lt);
    err_g_tau_QMC.resize(lt);

    KernelMat.resize(lt, nOmega);
}

void SAC::set_sampling_params(const double &theta, const int &n_constraint) {
    this->theta = theta;
    this->n_constraint = n_constraint;
}

void SAC::read_QMC_data(const std::string &filename) {
    /*
     *  Read imaginary-time QMC data (time-displaced Green's function) from file
     */

    std::ifstream infile;
    infile.open(filename, std::ios::in);

    if (!infile.is_open()) {
        std::cerr << "fail to open file " + filename + " !" << std::endl;
        exit(1);
    }

    std::string line;
    int t = 0;
    while(getline(infile, line)) {
        /* call boost library */
        std::vector<std::string> data;
        boost::split(data, line, boost::is_any_of(" "), boost::token_compress_on);
        tau_list(t) = boost::lexical_cast<double>(data[0]);
        g_tau_QMC(t) = boost::lexical_cast<double>(data[1]);
        err_g_tau_QMC(t) = boost::lexical_cast<double>(data[2]);
        ++t;
    }
    assert(t == lt);
    infile.close();
}

void SAC::read_Config_data(const std::string &filename) {
    /*
     *  Read weight configurations from file.
     *  if no file input, initialize configs A(\omega) according to a uniform distribution.
     */

    if (filename.empty()) {
        // no configs input, initialize A(\omega) as a flat spectrum
        std::cerr << " no configs input, initialize wight configs by mean distribution. " << std::endl;
        for (int i = 0; i < nOmega; ++i) {
            omega_list(i) = omegaMin + (i + 1) * deltaOmega;

            // normalized condition \sum A(\omega) * deltaOmega = 1.
            A_omega(i) = 1 / (deltaOmega * nOmega);
        }
    }

    else {
        std::ifstream infile;
        infile.open(filename, std::ios::in);

        if (!infile.is_open()) {
            std::cerr << "fail to open file " + filename + " !" << std::endl;
            exit(1);
        }
        else {
            std::string line;
            int i = 0;
            while(getline(infile, line)) {
                // specific infile format
                if (i < 6) { i++; continue;}

                std::vector<std::string> data;
                boost::split(data, line, boost::is_any_of(" "), boost::token_compress_on);
                omega_list(i - 6) = boost::lexical_cast<double>(data[1]);   // todo: check what is data[0] ???
                A_omega(i - 6) = boost::lexical_cast<double>(data[2]);
                ++i;
            }
            assert(i - 6 == nOmega);
            infile.close();
            std::cerr << "succeed to read configs from file " + filename + " !" << std::endl;
        }
    }
}

void SAC::initialSAC() {
    /*
     *  initialize temp vector and matrices for simulation use,
     *  after the imaginary-time QMC data or weight configurations have been read from file.
     */

    /** Important pre-work: read infile data first */
    // calculate kernel matrix
    for (int t = 0; t < lt; ++t) {
        const double tau = tau_list(t);
        for (int i = 0; i < nOmega; ++i) {
            const double omega = omega_list(i);
            KernelMat(t, i) = exp(-tau*omega) / (1 + exp(beta*omega));
        }
    }

    // initialize the fitting Green function g_tau and chi ^ 2
    cal_chi_square(A_omega, g_tau, chi_square);
}

void SAC::cal_chi_square(const vecXd &A, vecXd &g_tau_fitted, double &chi_2) {
    /*
     *  Calculate chi ^ 2 and the fitting g(\tau) from a given weight config A(\omega)
     *  Definition:
     *      chi ^ 2 = \sum ( g(\tau) - g_QMC(\tau) )^2 / \delta(\tau) ^ 2
     *      g (\tau) = \sum_{\omega} Kernel(\tau, \omega) * A(\omega)
     */
    assert(A.size() == nOmega);
    assert(g_tau_fitted.size() == lt);

    // update the fitting g(\tau)
    g_tau_fitted = KernelMat * A;

    // calculate chi ^ 2
    chi_2 = 0.0;
    for (int t = 0; t < lt; ++t) {
        chi_2 += pow((g_tau_fitted(t) - g_tau_QMC(t)), 2) / pow(err_g_tau_QMC(t), 2);
    }
}

void SAC::update_A_omega(vecXd &A_omega_new) {
    /*
     *  Randomly update weight configurations A(\omega) in such a way that
     *  the first few ( n_constraint ) frequency moments are conserved with fixed \rho(m):
     *
     *      \rho (m) = \sum_{i} \omega(i)^m * A(i) * [ 1 ± exp(-\beta * \omega(i))],    0 <= m < n_constraint
     *
     *  And (n_constraint + 1) wights are updated at the same time.
     */
    assert(A_omega_new.size() == nOmega);

    // randomly select ( n_constraint + 1 ) weights.
    std::vector<int> vector_aux(nOmega);
    std::iota(vector_aux.begin(), vector_aux.end(), 0);

    std::vector<int> select;
    std::sample(vector_aux.begin(), vector_aux.end(), std::back_inserter(select), n_constraint + 1, gen);

    vecXd A_select(n_constraint + 1);
    for (int i = 0; i < A_select.size(); ++i) {
        A_select(i) = A_omega_new(select[i]);
    }

    // construct augmented matrix
    matXd augMat(n_constraint, n_constraint + 2);
    for (int m = 0; m < n_constraint; ++m) {
        /*
         *  Constraint \rho (m), 0 <= m < n_constraint :
         *  C_m (i) = \omega(i)^m * ± exp(-\beta * \omega(i))]
         *  \rho (m) = \sum_{i} C_m(i) * A(i)
         */
        for (int i = 0; i < n_constraint + 1; ++i) {
            augMat(m, i) = pow(omega_list(select[i]), m) * ( 1 + pow(-1, m) * exp(- beta * omega_list(select[i])));
        }
    }

    const vecXd rho_select = augMat.block(0, 0 ,n_constraint, n_constraint + 1) * A_select;
    augMat.col(n_constraint + 1) = rho_select;

    // randomly select a point at a finite line segment in an (n_constraint + 1)-dimensional hypercube
    // Gaussian elimination method
    for (int m = 0; m < n_constraint; ++m) {
        for (int n = 0; n < m; ++n) {
            augMat.row(m) = augMat.row(m) - augMat.row(n) / augMat(n, n) * augMat(m, n);
        }
    }

    // randomly select  A_min(i) <= A(i) <= A_max(i)
    double A_min = 0.0, A_max = 0.0;
    for (int m = n_constraint - 1; m >= 0; --m) {
        for (int n = n_constraint - 1; n > m; --n) {
            augMat.row(m) = augMat.row(m) - augMat.row(n) * augMat(m, n) / augMat(n, n);
        }
        augMat.row(m) /= augMat(m, n_constraint);

        if (augMat(m, m) > 0)
            A_max = (A_max == 0.0)? augMat(m, n_constraint + 1) : fmin(A_max, augMat(m, n_constraint + 1));
        else
            A_min = fmax(A_min, augMat(m, n_constraint + 1));
    }

    // generate new wight configurations A(\omega)
    static std::uniform_real_distribution<double> u(0, 1);
    const double A_new = A_min + (A_max - A_min) * u(gen);
    for (int i = n_constraint; i >= 0; --i) {
        A_omega_new(select[i]) = (i==n_constraint)? A_new : (augMat(i, n_constraint + 1) - A_new) / augMat(i, i);
    }
}

void SAC::Metropolis_update() {
    /*
     *  Local update of weight configurations A(\omega) via Metropolis algorithm,
     *  according to a probability distribution proportional to
     *
     *      P (A) \propto exp(- chi^2 / \theta)
     *
     *  Once accepted, the configurations are updated in place.
     */
    total_step++;

    // copy from current configuration
    vecXd A_omega_new = A_omega;
    vecXd g_tau_new = g_tau;
    double chi_square_new = chi_square;

    // randomly update A(\omega)
    update_A_omega(A_omega_new);
    cal_chi_square(A_omega_new, g_tau_new, chi_square_new);

    // accept or reject new configuration according to standard Metropolis algorithm
    const double p = exp(-(chi_square_new - chi_square) / theta);
    if (std::bernoulli_distribution(std::min(1.0, p))(gen)) {
        // accepted
        A_omega = A_omega_new;
        g_tau = g_tau_new;
        chi_square = chi_square_new;
        accept_step++;
    }
}

void SAC::Metropolis_update_1step() {
    for (int n = 0; n < nOmega / (n_constraint + 1); ++n) {
        Metropolis_update();
    }
}


