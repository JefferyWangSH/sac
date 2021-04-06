#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

#include "sac.h"

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#include <Eigen/QR>


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

void SAC::read_QMC_data(const std::string& filename) {
    /*
     *  Read imaginary-time QMC data from file.
     */

    std::ifstream infile;
    infile.open(filename, std::ios::in);

    if (!infile.is_open()) {
        std::cerr << "Fail to open file " + filename + " !"<<  std::endl;
    }

    std::string line;
    int i = 0;
    while(getline(infile, line)) {
        /* call boost library */
        std::vector<std::string> data;
        boost::split(data, line, boost::is_any_of(" "), boost::token_compress_on);
        tau_list(i) = boost::lexical_cast<double>(data[0]);
        g_tau_QMC(i) = boost::lexical_cast<double>(data[1]);
        err_g_tau_QMC(i) = boost::lexical_cast<double>(data[2]);
        ++i;
    }
    assert(i == lt);
    infile.close();
}

void SAC::initialSAC() {
    /*
     *  Initialize temp vector and matrices for simulation use.
     */

    // initialize list of omega and A(\omega)
    for (int i = 0; i < nOmega; ++i) {
        omega_list(i) = omegaMin + (i + 1) * deltaOmega;

        // initialize A(\omega) as a flat spectrum, thus normalized condition \sum A(\omega) * deltaOmega = 1.
        // suppose system is invariant under space transformation, we here restrict omega to be positive.
        A_omega(i) = 1 / (2 * deltaOmega);
    }

    // calculate kernel matrix
    for (int t = 0; t < lt; ++t) {
     const int tau = tau_list(t);
        for (int i = 0; i < nOmega; ++i) {
            const int omega = omega_list(i);
            KernelMat(t, i) = (exp(-tau*omega) + exp(-(beta-tau)*omega)) / M_PI;
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
    // FIXME: check whether it is ergodic
    std::vector<int> vector_aux(nOmega);
    std::iota(vector_aux.begin(), vector_aux.end(), 0);

    std::vector<int> select;
    std::sample(vector_aux.begin(), vector_aux.end(), std::back_inserter(select), n_constraint + 1, gen);

    // perform augmented matrix
    matXd augment_Mat(n_constraint, n_constraint + 2);
    for (int m = 0; m < n_constraint; ++m) {
        /*
         *  Constraint \rho (m), 0 <= m < n_constraint :
         *  C_m (i) = \omega(i)^m * ± exp(-\beta * \omega(i))]
         *  \rho (m) = \sum_{i} C_m(i) * A(i)
         */
        for (int i = 0; i < n_constraint + 1; ++i) {
            augment_Mat(m, i) = pow(omega_list(select[i]), m) * ( 1 + pow(-1, m) * exp(- beta * omega_list(select[i])));
        }
    }

    vecXd A_select(n_constraint + 1);
    for (int i = 0; i < A_select.size(); ++i) {
        A_select(i) = A_omega_new(select[i]);
    }

    const vecXd rho_select = augment_Mat.block(0, 0 ,n_constraint, n_constraint + 1) * A_select;
    augment_Mat.col(n_constraint + 1) = rho_select;

    // FIXME: use FullPivHouseholderQR: more accurate but slower
    // randomly select a point at a finite line segment in an (n_constraint + 1)-dimensional hypercube
    Eigen::HouseholderQR<matXd> qr;
    qr.compute(augment_Mat);
    matXd R = qr.matrixQR().triangularView<Eigen::Upper>();

    // randomly select  A_min(i) <= A(i) <= A_max(i)
    double A_min = 0.0, A_max = 0.0;
    for (int m = n_constraint - 1; m >= 0; --m) {
        for (int n = n_constraint - 1; n > m; --n) {
            R.row(m) = R.row(m) - R.row(n) * R(m, n) / R(n, n);
        }
        R.row(m) /= R(m, n_constraint);

        if (R(m, m) > 0)
            A_max = (A_max == 0.0)? R(m, n_constraint + 1) : fmin(A_max, R(m, n_constraint + 1));
        else
            A_min = fmax(A_min, R(m, n_constraint + 1));
    }

    // generate new wight configurations A(\omega)
    static std::uniform_real_distribution<double> u(0, 1);
    const double A_new = A_min + (A_max - A_min) * u(gen);
    for (int i = n_constraint; i >= 0; --i) {
        A_omega_new(select[i]) = (i==n_constraint)? A_new : (R(i, n_constraint + 1) - A_new) / R(i, i);
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

    // randomly update A(\omega)
    vecXd A_omega_new = A_omega;
    vecXd g_tau_new = g_tau;
    double chi_square_new = chi_square;

    update_A_omega(A_omega_new);
    cal_chi_square(A_omega_new, g_tau_new, chi_square_new);

    // accept or reject new configuration according to standard Metropolis algorithm
    const double p = exp(-(chi_square_new - chi_square) / theta);
    if (std::bernoulli_distribution(std::min(1.0, p))(gen)) {
        // accepted
        A_omega = A_omega_new;
        g_tau = g_tau_new;
        chi_square = chi_square_new;
    }
}

void SAC::Metropolis_update_1step() {
    for (int n = 0; n < nOmega / (n_constraint + 1); ++n) {
        Metropolis_update();
    }
}
