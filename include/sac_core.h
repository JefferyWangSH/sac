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

#include <string>
#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#include <Eigen/Core>
#include "sac_annealing.h"

// forward declaration
namespace Grids { class FreqGrids; }


namespace SAC {

    // forward declarations
    class Kernel;
    class Measure;
    namespace Initializer { class QmcReader; }

    // -------------------------------------------  SAC::SacCore class  ----------------------------------------------
    class SacCore {

        private:

            // information about the QMC data
            int m_time_size{};                    // number of imaginary-time points
            double m_beta{};                      // inverse temperature
            double m_scaling_factor{};            // scaling factor G(t=0)

            // SAC parameters
            int m_delta_num{};                    // number of delta functions
            double m_delta_amplitude{};           // amplitude of the delta peak
            double m_annealing_rate{};            // pace of annealing (annealing rate)
            int m_stabilization_pace{};           // pace of stabilization (recompute chi2)

            Eigen::VectorXd m_tgrids_qmc{};       // imaginary-time points from QMC (processed)
            Eigen::VectorXd m_corr_qmc{};         // correlation functions from QMC (processed)
            Eigen::VectorXd m_sigma_qmc{};        // standard deviation of the transformed correlations

            Eigen::VectorXd m_corr_now{};         // current correlation functions
            Eigen::VectorXd m_corr_next{};        // updated correlation functions

            int m_collecting_steps{};             // number of MC steps for collecting the spectrum
            Eigen::VectorXd m_freq{};             // frequency list of the recovered spectral function
            Eigen::VectorXd m_spec{};             // recovered spectral function

            double m_accept_ratio{};              // averaged accepting ratio of MC moves
            double m_chi2{};                      // chi2 (goodness of fitting) of the current spectral configs
            double m_chi2_min{};                  // minimum of chi2

            std::string m_update_type{};          // modes of MC update (single or pair)
            std::string m_kernel_type{};          // type of the SAC kernel

            // `Annealing::Metedata` class containing the basic data used for SAC updates and sampling
            // including:
            //   1. sampling temperature `theta`
            //   2. window width of random delta function moves `window_width`
            //   3. goodness of fittings `chi2` 
            //   4. locations of delta functions
            // the changes of these data are packed up and stored in memory,
            // serving as a reference to the history of the annealing process.
            Annealing::MetaData m_metadata{};


        public:

            SacCore() = default;

            // --------------------------------------  Set up SAC parameters  ---------------------------------------

            // set up parameters for the Monte Carlo process
            void set_sampling_params( int delta_num, int collecting_steps, int stabilization_pace, 
                                      const std::string& update_type = "single" );

            // set up paramters for the control of the annealing process
            void set_annealing_params( double theta, double annealing_rate );
            
            // initialize SacCore from QmcReader module
            void initial( const Kernel& kernel, 
                          const Initializer::QmcReader& qmc_reader, 
                          const Grids::FreqGrids& grids );

            
            // -------------------------------  Crucial member functions for SAC  -----------------------------------

            // annealing process
            void perform_annealing( const Kernel& kernel, const Grids::FreqGrids& grids , 
                                    Measure& measure, Annealing::Chain& chain );

            // determine the sampling temperature after the annealing process
            // by slightly raising the artificial temperature `theta` to avoid overfitting
            void decide_sampling_theta( const Kernel& kernel, const Annealing::Chain& chain );

            // sampling and collecting the spectrum
            void sample_and_collect( const Kernel& kernel, const Grids::FreqGrids& grids , Measure& measure );


            // ----------------------------------  Interface member functions  --------------------------------------
            
            const Eigen::VectorXd& FrequencyGrids()    const;
            const Eigen::VectorXd& RecoveredSpectrum() const;
            double FrequencyGrids( int i )    const;
            double RecoveredSpectrum( int i ) const;


        private:
            
            // -----------------------------------  Private member functions  ---------------------------------------

            // compute correlation functions from current configurations of the delta functions
            void compute_corr_from_deltas( const Kernel& kernel );

            // compute chi2, the goodness of fitting, for any input correlations
            double compute_goodness( const Eigen::VectorXd& computed_corr ) const;

            // one step of Monte Carlo update of the delta functions
            void update_deltas_1step( const Kernel& kernel, const Grids::FreqGrids& grids );

            // move a single delta function in an update attempt
            void update_deltas_1step_single( const Kernel& kernel, const Grids::FreqGrids& grids );

            // move a pair of delta functions in an update attempt
            void update_deltas_1step_pair( const Kernel& kernel, const Grids::FreqGrids& grids );

            // perform MC updates at a fixed artificial temperature `theta`
            // till reaching the equilibrium state
            void update_at_fixed_theta( const Kernel& kernel, const Grids::FreqGrids& grids , Measure& measure );

    };

} // namespace SAC

#endif // SAC_SACCORE_H
