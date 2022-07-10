#ifndef SAC_SACCORE_H
#define SAC_SACCORE_H
#pragma once

/**
  *  This header file defines `SAC::SacCore` class as a central module of the SAC method.
  *  The stochastic analytic continuation method (SAC) 
  *  was first proposed by Anders.W. Sandvik in 1998. 
  *  It performs a Monte Carlo simulation combining with a simulated annealing process
  *  to extract the real-frequency spectral information 
  *  from imaginary-time QMC correlation functions.
  */

#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#include <Eigen/Core>


namespace SAC {

    // -------------------------------------------  SAC::SacCore class  ----------------------------------------------
    // Todo: rename Core to SacCore
    class Core {

        private:

            // information about the QMC data
            int m_time_size{};                    // number of imaginary-time points
            double m_beta{};                      // inverse temperature
            double m_scaling_factor{};            // scaling factor G(t=0)

            double m_annealing_pace{};            // pace of annealing (annealing rate)
            int m_stablization_pace{};            // pace of stablization (recompute chi2)

            Eigen::VectorXd m_tgrids_qmc{};       // imaginary-time points from QMC (processed)
            Eigen::VectorXd m_corr_qmc{};         // correlation functions from QMC (processed)
            Eigen::VectorXd m_sigma_qmc{};        // standard deviation of the transformed correlations

            Eigen::VectorXd m_corr_now{};         // current correlation functions
            Eigen::VectorXd m_corr_next{};        // updated correlation functions

            int m_collecting_steps{};             // number of MC steps for collecting spectrum
            Eigen::VectorXd m_freq{};             // frequency list of the recovered spectral function
            Eigen::VectorXd m_spec{};             // recovered spectral function

            double m_accept_ratio{};              // averaged accepting ratio of MC moves
            double m_chi2{};                      // chi2 (goodness of fitting) of the current spectral configs
            double m_chi2_min{};                  // minimum of chi2

        public:

            Core() = default;




    };

} // namespace SAC

#endif // SAC_SACCORE_H
