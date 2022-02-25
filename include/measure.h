#ifndef SAC_MEASURE_H
#define SAC_MEASURE_H
#pragma once

/**
  *  This header file includes `Measure` class for monitoring and measuring `chi2` and average accepting rate
  *  during the MC update of delta functions at a given fictitious temperature `theta`.
  *  Bin calculation of mean values and error esitimation are implemented.
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
        int size_of_bin{};

        // measurements of chi2 and accepting radio
        double chi2_mean{}, chi2_err{};
        double accept_radio_mean{}, accept_radio_err{};

        // data samples collected during Monte Carlo process, index labels the time sequence of collection
        Eigen::VectorXd chi2_sample{};
        Eigen::VectorXd accept_radio_sample{};

        // means of data sorted by bin index
        Eigen::VectorXd chi2_bin{};
        Eigen::VectorXd accept_radio_bin{};

    public:
        Measure() = default;

        Measure(int nbin, int size_of_bin);

        /* resize dimensions */
        void resize(int nbin, int size_of_bin);

        /* collect data */
        void collect(int s, double chi2, double accept_radio);

        /* analyse measurements for one certain bin */
        void bin_analyse(int n);

        /* compute means and errors using bin analysis */
        void analyse();

        /* clear previous statistics data */
        void clear();

        /* return mean value of chi2 */
        double chi2() const;

        /* return mean value of accepting radio */
        double accept_radio() const;
    };

} // namespace Measure

#endif //SAC_MEASURE_H
