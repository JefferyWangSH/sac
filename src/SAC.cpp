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


void SAC::set_SAC_params(int _lt, double _beta, int _n_config, double _omega_min, double _omega_max, int _n_moment) {
    /*
     *  Set SAC parameters.
     */
    assert( _lt > 0 && _beta > 0 );
    assert( _n_config > 0 );
    assert( _omega_min < _omega_max );

    this->lt = _lt;
    this->beta = _beta;
    this->dtau = _beta / _lt;
    this->n_config = _n_config;
    this->omega_min = _omega_min;
    this->omega_max = _omega_max;
    this->delta_n_config = 1.0 / (_n_config + 1);
    this->n_moment = _n_moment;
}

void SAC::set_alpha(const double &_alpha) {
    this->alpha = _alpha;
}

void SAC::set_QMC_filename(const std::string &_filename) {
    this->filename = _filename;
}

void SAC::read_QMC_data(const std::string &_filename) {
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
        std::vector<std::string> data;
        boost::split(data, line, boost::is_any_of(" "), boost::token_compress_on);
        data.erase(std::remove(std::begin(data), std::end(data), ""), std::end(data));

        g_tau[t] = boost::lexical_cast<double>(data[1]);
        sigma_tau[t] = boost::lexical_cast<double>(data[2]);
        ++t;
    }
    // check the consistence between data and model settings
    assert(t == lt);
    infile.close();
}

/*
void SAC::read_Config_data(const std::string &filename) {
    // Read weight configurations from file.
    // if no file input, initialize configs A(\omega) according to a uniform distribution.


    if (filename.empty()) {
        // no configs input, initialize n(x) as a flat spectrum
        std::cerr << " no configs input, initialize wight configs by mean distribution. " << std::endl;
        for (int i = 0; i < n_config; ++i) {
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
*/

void SAC::prepare() {
    /*
     *  allocate memory for SAC and
     *  initialize temp vector and matrices for simulation use,
     *  TODO: if no configs file input,
     *        initialize configs n(x) according to a uniform distribution.
     */

    // clear previous data
    x_list.clear();
    n_list.clear();
    omega_list.clear();
    A_config.clear();
    tau_list.clear();
    g_tau.clear();
    sigma_tau.clear();
    h_tau.clear();
    kernel.clear();

    // reserve memory for vectors
    x_list.reserve(n_config);
    n_list.reserve(n_config);
    omega_list.reserve(n_config);
    A_config.reserve(n_config);

    for (int i = 0; i < n_config; ++i){
        x_list.emplace_back(0.0);
        n_list.emplace_back(0.0);
        omega_list.emplace_back(0.0);
        A_config.emplace_back(0.0);
    }

    tau_list.reserve(lt);
    g_tau.reserve(lt);
    sigma_tau.reserve(lt);
    h_tau.reserve(lt);
    kernel.reserve(lt);

    for (int t = 0; t < lt; ++t) {
        tau_list.emplace_back(0.0);
        g_tau.emplace_back(0.0);
        sigma_tau.emplace_back(0.0);
        h_tau.emplace_back(0.0);
        kernel.emplace_back(n_config, 0.0);
    }

    // shrink to fit
    x_list.shrink_to_fit();
    n_list.shrink_to_fit();
    omega_list.shrink_to_fit();
    A_config.shrink_to_fit();
    tau_list.shrink_to_fit();
    g_tau.shrink_to_fit();
    sigma_tau.shrink_to_fit();
    h_tau.shrink_to_fit();
    kernel.shrink_to_fit();

    // initialize tau
    for (int t = 0; t < lt; ++t) {
        tau_list[t] = (double)(t + 1) / lt * beta;
    }

    // initialize n(x) with a flat spectrum
    for (int i = 0; i < n_config; ++i) {
        x_list[i] = delta_n_config * (i + 1);
        omega_list[i] = omega_min + ( omega_max - omega_min ) * x_list[i];

        // normalized condition \sum_i n(i) = 1.
        n_list[i] = 1.0 / n_config;
    }

    /** QMC data should be read from file before initialize kernel. */
    read_QMC_data(filename);

    // calculate kernel matrix
    for (int t = 0; t < lt; ++t) {
        const double tau = tau_list[t];
        for (int i = 0; i < n_config; ++i) {
            const double omega = omega_list[i];
            kernel[t][i] = exp(- tau * omega) / (1 + exp(- beta * omega));
        }
    }

    // calculate hamiltonian density h(\tau) for default configurations
    cal_Config_Hamiltonian_Density(h_tau);
}

void SAC::cal_Config_Hamiltonian(double &H) {
    /*
     *  Effective Hamiltonian for configuration field:
     *     H[n(x)] = \int d\tau h(\tau)^2
     *  where h(\tau) = { ( \sum_{x} K(\tau, x) * n(x) ) - g(\tau) } / \sigma(\tau)
     */

    // make sure hamiltonian density h(\tau) previously updated
    H = 0;
    for (int t = 0; t < lt; ++t) {
        // FIXME: Mind integrate over tau, here could missing factor \delta \tau
        H += h_tau[t] * h_tau[t] * dtau;
    }
}


void SAC::cal_Config_Hamiltonian_Density(std::vector<double> &h) {
    /*
     *  Hamiltonian density for a specific field configuration C
     *      h(\tau) = { ( \sum_{x} K(\tau, x) * n(x) ) - g(\tau) } / \sigma(\tau)
     *  Critical quantity during Metropolis updates.
     */
    assert(h.size() == lt);

    for (int t = 0; t < lt; ++t) {
        h[t] = 0.0;
        for (int i = 0; i < n_config; ++i) {
            h[t] += kernel[t][i] * n_list[i];
        }
        h[t] = ( h[t] - g_tau[t] ) / sigma_tau[t];
    }
}

void SAC::update_Configs(std::vector<int> &index_selected, std::vector<double> &n_selected) {
    /*
     *  Randomly update weigh configs by conserving moments to accelerate convergence.
     *    1. Firstly, we randomly choice ( n_moment + 1 ) points in space of configs.
     *    2. Then we randomly select a point on the 1-dimensional line of constraint
     *       through the ( n_moment + 1 )-dimensional space of configs by conserving the first (n_moment) moments
     *          M^(n) = \int_{0}^{1}  dx n(x) x^{n} ,   0 <= n < n_moment
     *  Mind the density h(\tau) is changed in place.
     */
    assert( index_selected.size() == n_moment + 1 );
    assert( n_selected.size() == n_moment + 1 );

    // randomly select ( n_moment + 1 ) points in space of configs.
    std::vector<int> vector_aux(n_config);
    std::iota(vector_aux.begin(), vector_aux.end(), 0);

    std::vector<int> select;
    std::sample(vector_aux.begin(), vector_aux.end(), std::back_inserter(select), n_moment + 1, generate_SAC);

    // selected configs, we will then update them
    std::vector<double> n_selected_new(n_moment + 1);
    for (int i = 0; i < n_moment + 1; ++i) {
        n_selected_new[i] = n_list[select[i]];
    }

    // selected configs are updated in such a way that
    //   n’(x) = n(x) - s * Q(x)
    // where s parameterizes the 1-dimensional line in space of configs
    std::vector<double> Q(n_moment + 1);
    std::uniform_int_distribution<int> uniform_int(0, n_moment);
    const int basic = uniform_int(generate_SAC);

    for (int i = 0; i < n_moment + 1; ++i) {
        if ( i == basic ) {
            Q[i] = -1;
            continue;
        }
        Q[i] = 1.0;
        for (int j = 0; j < n_moment + 1; ++j) {
            if ( j != i && j != basic ) {
                Q[i] *= ( x_list[select[j]] - x_list[select[basic]] ) / ( x_list[select[j]] - x_list[select[i]] );
            }
        }
    }

    // parameterized variable s
    double s, s_max = 0.0, s_min = 0.0;
    for (int i = 0; i < n_moment + 1; ++i) {
        if ( Q[i] < 0.0 ) {
            s_min = ( s_min == 0 )? n_selected_new[i] / Q[i] : std::max(s_min, n_selected_new[i] / Q[i]);
        }
        if ( Q[i] > 0.0 ) {
            s_max = ( s_max == 0 )? n_selected_new[i] / Q[i] : std::min(s_max, n_selected_new[i] / Q[i]);
        }
    }
    std::uniform_real_distribution<double> uniform_double(0, 1);
    s = s_min + (s_max - s_min) * uniform_double(generate_SAC);

    // update selected configs
    for (int i = 0; i < n_moment + 1; ++i) {
        n_selected_new[i] -= s * Q[i];
    }

    for (int i = 0; i < n_moment + 1; ++i) {
        n_selected[i] = n_selected_new[i];
        index_selected[i] = select[i];
    }
}

void SAC::Metropolis_update() {
    /*
     *  Metropolis update of configuration fields
     *  if new configs are accepted, the corresponding fields n(x) and h(\tau) are updated in place.
     */

    // we try to randomly update configs in a moment-conserved manner
    std::vector<int> index(n_moment + 1);
    std::vector<double> n_new(n_moment + 1);
    update_Configs(index, n_new);

    // calculate difference of h(\tau): delta_h(\tau)
    std::vector<double> delta_h(lt, 0.0);
    for (int t = 0; t < lt; ++t) {
        for (int i = 0; i < n_moment + 1; ++i) {
            delta_h[t] += kernel[t][i] * ( n_new[index[i]] - n_list[index[i]] );
        }
        delta_h[t] /= sigma_tau[t];
    }

    // calculate differences of hamiltonian: delta_H
    double delta_H = 0.0;
    for (int t = 0; t < lt; ++t) {
        // FIXME: Mind integrate over tau, here could missing factor \delta \tau
        delta_H += ( delta_h[t] * delta_h[t] + 2 * delta_h[t] * h_tau[t]) * dtau;
    }

    // accept or reject new configs according to standard Metropolis algorithm
    // if accepted, update n(x) and h(\tau) in place.
    const double p = exp( - alpha * delta_H );
    if (std::bernoulli_distribution(std::min(1.0, p))(generate_SAC)) {
        // updates accepted
        for (int i = 0; i < n_moment + 1; ++i) {
            n_list[index[i]] = n_new[i];
        }
        for (int t = 0; t < lt; ++t) {
            h_tau[t] +=  delta_h[t];
        }
    }
}

void SAC::Metropolis_update_1step() {
    /*
     *  We define one Monte Carlo step as in which we roughly update all sites in the space of configs
     */
    for (int n = 0; n < ceil((double) n_config / (n_moment + 1)); ++n) {
        Metropolis_update();
    }
}

