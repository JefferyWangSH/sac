#ifndef SAC_MEASURE_H
#define SAC_MEASURE_H
#pragma once

/**
  *  This header file includes `Measure` class for monitoring and measuring `chi2` and average accepting rate
  *  during the MC update of delta functions at a given fictitious temperature `theta`.
  *  Bin calculation of means and errors are performed in this file.
  */

#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#include <Eigen/Core>


namespace Measure{
    class Measure {

    public:

        // number of bins
        int nbin{};

        // bin size, number of samples within one bin
        int sbin{};

        // measurements of chi2 and accepting radio
        double chi2_mean{}, chi2_err{};
        double accept_radio_mean{}, accept_radio_err{};

        // data samples collected from Monte Carlo, index labels collection sequence
        Eigen::VectorXd sample_chi2;
        Eigen::VectorXd sample_accept_radio;

        // means of data sorted by number of bin
        Eigen::VectorXd bin_chi2;
        Eigen::VectorXd bin_accept_radio;

    public:

        Measure() = default;

        Measure(int nbin, int sbin);

        /* resize dimensions */
        void resize(int nbin, int sbin);

        /* fill data */
        void fill(int s, double chi2, double accept_radio);

        /* analyse measurements for one certain bin */
        void bin_analyse(int n);

        /* compute means and errors using bin analysis */
        void analyse();

        /* clear previous statistics data */
        void clear();

        /* return mean of chi2 */
        const double chi2() const;

        /* return mean of accepting radio */
        const double accept_radio() const;
    };
}

#endif //SAC_MEASURE_H
