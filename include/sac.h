#ifndef SAC_H
#define SAC_H
#pragma once

/**
  *  This head file includes SAC class for the implementation of stochastic analytic continuation method,
  *  proposed by Anders.W. Sandvik.
  *  A Monte Carlo process and simulated annealing are performed
  *  to extract real-frequency spectral information from imaginary-time correlation functions,
  *  which can be obtained previously by QMC calculations.
  */

#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#include <Eigen/Core>

namespace Grids { class FreqGrids; }
namespace Kernel { class Kernel; }
namespace Annealing { class AnnealingData; class AnnealingChain; }
namespace DataReader { class QMCDataReader; }
namespace Measure { class Measure; }


namespace Simulation {
    
    class SAC {
    public:
        /* model params */
        int nt{};                           // number of time slices
        double beta{};                      // inverse temperature
        double scale_factor{};              // scaling factor G(0)

        std::string update_type{};          // modes of MC update (single or pair)
        std::string kernel_type{};          // kernel type

        Eigen::VectorXd tau_from_qmc{};     // tau points from QMC
        Eigen::VectorXd corr_from_qmc{};    // correlation functions from QMC
        Eigen::VectorXd sigma_from_qmc{};   // standard deviation of transformed correlations

        Eigen::VectorXd corr_now{};         // current correlation function
        Eigen::VectorXd corr_next{};        // updated correlation function

        int collecting_steps{};             // number of MC steps for collecting spectrum
        Eigen::VectorXd freq{};             // recovered frequency
        Eigen::VectorXd spec{};             // recovered spectrum

        double accept_radio{};              // average accepting radio of MC move
        double chi2{};                      // chi2 (goodness of fitting) of current spectrum
        double chi2_min{};                  // minimum of chi2

        std::string log_name{};             // log file name

        /* griding params */
        Grids::FreqGrids *grids{};

        /* sampling params */
        Annealing::AnnealingData *annealing_data{};
        Annealing::AnnealingChain *annealing_chain{};

        /* kernel */
        Kernel::Kernel *kernel{};

        /* data from QMC input */
        DataReader::QMCDataReader *qmc_data_reader{};

        /* measuring module */
        Measure::Measure *measure{};

    public:
        /* construction and destruction */
        SAC();
        ~SAC();

        /** subroutine for parameter settings */
        /* set up parameters for read in module */
        void set_read_in_params(int lt, double beta, int nbin, int rebin_pace, int bootstrap_num);

        /* set up file which contains data of tau points */
        void set_filename_tau(const std::string &infile_tau);

        /* set up file which contains data of time correlations */
        void set_filename_corr(const std::string &infile_corr);

        /* set up file name of log output */
        void set_filename_log(const std::string &outfile_log);

        /* set up parameters for grids of frequency domain */
        void set_griding_params(double freq_interval, double spec_interval, double freq_min, double freq_max);

        /* set up parameters for sampling procedure */
        void set_sampling_params(int ndelta, double theta, int max_annealing_steps, int bin_num, int bin_size, int collecting_steps);

        /* set up parameters controlling simulation modes */
        void set_mode_params(const std::string &kernel_type, const std::string &update_type);

        /* initialization */
        void init();

        /* annealing process */
        void perform_annealing();

        /* determine sampling temperature after annealing */
        void decide_sampling_theta();

        /* sampling and collect spectrum */
        void sample_and_collect();

        /* file output of recovered spectrum */
        void output(const std::string &filename);

    private:
        /* read and preprocessing QMC data */
        void init_from_module();

        /* initialize spectrum */
        void init_spectrum();

        /* computing current correlations from spectrum */
        void compute_corr_from_spec();

        /* compute chi2, the goodness of fitting, for any input correlations */
        double compute_goodness(const Eigen::VectorXd &corr_from_spectrum) const;

        /* one step of Monte Carlo update of delta functions */
        void update_deltas_1step();

        /* move one single delta function in an attempted update */
        void update_deltas_1step_single();

        /* move one pair of delta functions in an attempted update */
        void update_deltas_1step_pair();

        /* equilibrium of system at a fixed theta */
        void update_fixed_theta();

        /* log output: n labels index of bin */
        void write_log(int n);
    };

} // namespace Simulation

#endif //SAC_H
