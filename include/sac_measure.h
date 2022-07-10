#ifndef SAC_MEASURE_H
#define SAC_MEASURE_H
#pragma once

/**
  *  This header file defines `SAC::Measure` class for measuring 
  *  `chi2` and the averaged accepting rate of new configurations
  *  during the MC update of delta functions at specific fictitious temperature `theta`.
  *  Bin calculations of the mean value and error-bar are implemented as member functions.
  */

#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#include <Eigen/Core>


namespace SAC {
    
    // --------------------------------------  SAC::Measure class  ----------------------------------------
    class Measure {
        
        private:
            // number of bins
            int m_number_of_bin{};

            // bin size, number of samples within one bin
            int m_size_of_bin{};

            // measurements of chi2 and accepting ratio
            double m_chi2_mean{}, m_chi2_err{};
            double m_accept_ratio_mean{}, m_accept_ratio_err{};

            // samples of metadata which is collected during the Monte Carlo process
            Eigen::VectorXd m_chi2_sample{};
            Eigen::VectorXd m_accept_ratio_sample{};

            // mean values sorted by bin indices
            Eigen::VectorXd m_chi2_bin{};
            Eigen::VectorXd m_accept_ratio_bin{};

        public:

            Measure() = default;
            Measure( int number_of_bin, int size_of_bin );

            // resize dimensions
            void resize( int number_of_bin, int size_of_bin );

            // collect data
            // s labels the index of sample
            void collect( int s, double chi2, double accept_ratio );

            // analyse the collected samples for one certain bin
            // n labels the index of bin
            void bin_analyse( int n );

            // compute the mean value and error-bar using bin_analysis()
            void analyse();

            // clear previous statistical data
            void clear();

            // interface member functions
            int number_of_bin() const;
            int size_of_bin()   const;

            double chi2()                const;
            double chi2_error()          const;
            double chi2( int n )         const;
            double accept_ratio()        const;
            double accept_ratio_error()  const;
            double accept_ratio( int n ) const;

    };

} // namespace SAC

#endif //SAC_MEASURE_H
