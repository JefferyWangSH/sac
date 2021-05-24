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

    std::vector<std::vector<double>> bin_H_alpha;
    std::vector<double> H_alpha, err_H_alpha;
    std::vector<double> tmp_H_alpha;

    std::vector<std::vector<std::vector<double>>> bin_n_x_alpha;
    std::vector<std::vector<double>> n_x_alpha, err_n_x_alpha;
    std::vector<std::vector<double>> tmp_n_x_alpha;

    // Todo: support keep track of swap rate p


public:

    Measure() = default;

    explicit Measure(int nbin, int nalpha, int nconfig);

    void resize(int nbin, int nalpha, int nconfig);

    void prepare();

    void clear_tmp_stats();

    void measure(MonteCarloSAC &sac);

    void normalize_stats();

    void write_data_to_bin(const int &bin);

    void analyse_stats();

private:

};

#endif //STOCHASTIC_ANALYTIC_CONTINUATION_MEASURE_H
