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

#include <memory>
#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#include <Eigen/Core>

#include "qmc_data_reader.h"
#include "freq_grids.h"
#include "kernel.h"
#include "annealing_chain.h"
#include "measure.h"


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

        std::string log_file_path{};        // path of logging out file
        std::string spec_file_path{};       // path of spectrum out file
        std::string report_file_path{};     // path of quality report file

        /* griding params */
        std::unique_ptr<Grids::FreqGrids> grids{};

        /* sampling params */
        std::unique_ptr<Annealing::AnnealingData> annealing_data{};
        std::unique_ptr<Annealing::AnnealingChain> annealing_chain{};

        /* kernel */
        std::unique_ptr<Kernel::Kernel> kernel{};

        /* data from QMC input */
        std::unique_ptr<DataReader::QMCDataReader> qmc_data_reader{};

        /* measuring module */
        std::unique_ptr<Measure::Measure> measure{};

    public:
        /* construction function */
        SAC() = default;

        /** subroutine for parameter settings */
        /* set up parameters for read in module */
        void set_qmc_reader_params(int lt, double beta, int nbin, int rebin_pace, int bootstrap_num);

        /* set up file which contains data of tau points */
        void set_file_path_tau(const std::string &tau_file_path);

        /* set up file which contains data of time correlations */
        void set_file_path_corr(const std::string &corr_file_path);

        /* set up path of output file, including output of logs, recovered spectrum and quality report */
        void set_outfile_path(const std::string &log_file_path, const std::string &spec_file_path, const std::string &report_file_path);

        /* set up parameters for grids of frequency domain */
        void set_griding_params(double freq_interval, double spec_interval, double freq_min, double freq_max);

        /* set up parameters for sampling procedure */
        void set_sampling_params(int ndelta, double theta, int max_annealing_steps, int bin_num, int bin_size, int collecting_steps);

        /* set up kernel type */
        void set_kernel_type(const std::string &kernel_type);

        /* set up updating type */
        void set_update_type(const std::string &update_type);

        /* initialization */
        void init();

        /* annealing process */
        void perform_annealing();

        /* determine sampling temperature after annealing */
        void decide_sampling_theta();

        /* sampling and collect spectrum */
        void sample_and_collect();

        /* file output of recovered spectrum */
        void output_recovered_spectrum();

        /* output quality information about the recovered spectrum */
        void report_recovery_quality();

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
