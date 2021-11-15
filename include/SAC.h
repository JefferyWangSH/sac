#ifndef SAC_H
#define SAC_H
#pragma once

/**
  *  This head file includes SAC class for the implementation of stochastic analytic continuation method,
  *  proposed by Anders.W. Sandvik.
  *  A Monte Carlo process and simulated annealing are performed
  *  to extract real-frequency information from imaginary-time correlations,
  *  which are obtained previously by QMC calculations.
  */

#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#include <Eigen/Core>

namespace Grid { class FrequencyGrid; }
namespace Kernel { class Kernel; }
namespace Annealing { class AnnealData; class AnnealChain; }
namespace QMCData { class ReadInModule; }
namespace Measure { class Measure; }


namespace Simulation {
    
    class SAC {
    public:

        /* model params */
        int nt{};                           // number of time slices
        double beta{};                      // inverse temperature
        double scale_factor{};              // scaling factor G(0)

        std::string update_mode{};          // modes of MC update (single or pair)
        std::string kernel_mode{};          // kernel types

        Eigen::VectorXd tau;                // tau points
        Eigen::VectorXd corr;               // time correlations from QMC
        Eigen::VectorXd sigma;              // standard deviation of transformed correlations

        Eigen::VectorXd corr_current;       // current time correlations from spectrum
        Eigen::VectorXd corr_update;        // updated time correlations from spectrum

        int collecting_steps{};             // number of MC steps for collecting spectrum
        Eigen::VectorXd freq;               // recovered frequency
        Eigen::VectorXd spectrum;           // recovered spectrum

        double accept_radio{};              // average accepting radio of MC move
        double chi2{};                      // chi2 (goodness of fitting) of current spectrum
        double chi2_minimum{};              // minimum of chi2

        std::string log_name{};             // log file name

        /* griding params */
        Grid::FrequencyGrid *grid;

        /* sampling params */
        Annealing::AnnealData *data;
        Annealing::AnnealChain *anneal;

        /* kernel */
        Kernel::Kernel *kernel;

        /* data from QMC input */
        QMCData::ReadInModule *readin;

        /* measuring module */
        Measure::Measure *measure;

    public:

        /* construction and destruction */
        SAC();
        ~SAC();

        /** subroutine for parameter settings */
        /* set up parameters for read in module */
        void set_read_in_params(int lt, double beta, int nbin, int rebin_pace, int num_bootstrap);

        /* set up file which contains data of tau points */
        void set_filename_tau(const std::string &infile_tau);

        /* set up file which contains data of time correlations */
        void set_filename_corr(const std::string &infile_corr);

        /* set up file name of log output */
        void set_filename_log(const std::string &outfile_log);

        /* set up parameters for grids of frequency domain */
        void set_griding_params(double grid_interval, double spec_interval, double omega_min, double omega_max);

        /* set up parameters for sampling procedure */
        void set_sampling_params(int ndelta, double theta, int max_annealing_steps, int bin_num, int bin_size, int collecting_steps);

        /* set up parameters controlling simulation modes */
        void set_mode_params(const std::string &kernel_mode, const std::string &update_mode);

        /** initialization */
        void init();

        /** annealing process */
        void perform_annealing();

        /* determine sampling temperature after annealing */
        void decide_sampling_theta();

        /** sampling and collect spectrum */
        void sample_and_collect();

        /** file output of recovered spectrum */
        void output(const std::string &filename);

    private:

        /* read QMC data (transformed) from read in module */
        void init_from_module();

        /* initialize spectrum */
        void init_spectrum();

        /* computing current correlations from spectrum */
        void compute_corr_from_spec();

        /* compute chi2, the goodness of fitting, for any input correlations */
        const double compute_goodness(const Eigen::VectorXd &corr_from_spectrum) const;

        /** one step of Monte Carlo update of delta functions */
        void update_deltas_1step();

        /* move a single delta function in one moving attempt */
        void update_deltas_1step_single();

        /* move a pair of delta functions in one moving attempt */
        void update_deltas_1step_pair();

        /** equilibrium of system at a fixed theta */
        void update_fixed_theta();

        /* log output: n labels index of bin number */
        void write_log(int n);

    };
}


#endif //SAC_H
