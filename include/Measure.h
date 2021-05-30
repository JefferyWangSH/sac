#ifndef STOCHASTIC_ANALYTIC_CONTINUATION_MEASURE_H
#define STOCHASTIC_ANALYTIC_CONTINUATION_MEASURE_H
#pragma once

#include <vector>

class MonteCarloSAC;

class Measure {

public:

    int nbin{};
    int nalpha{};
    int nconfig{};

    int nn{0};

    std::vector<double> bin_Hamilton;
    double mean_Hamilton{}, err_Hamilton{};
    double tmp_Hamilton{};

    std::vector<std::vector<double>> bin_n_x;
    std::vector<double> mean_n_x, err_n_x;
    std::vector<double> tmp_n_x;

public:

    Measure() = default;

    explicit Measure(int nbin, int nconfig);

    void resize(int nbin, int nconfig);

    void prepare();

    void clear_tmp_stats();

    void measure(MonteCarloSAC& sac);

    void normalize_stats();

    void write_data_to_bin(const int& bin);

    void analyse_stats();

};

#endif //STOCHASTIC_ANALYTIC_CONTINUATION_MEASURE_H
