#include "sac.h"
#include "qmc_data_reader.h"
#include "freq_grids.h"
#include "kernel.h"
#include "annealing_chain.h"
#include "measure.h"
#include "random.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <boost/format.hpp>

namespace Simulation {

    SAC::SAC() {
        this->qmc_data_reader = new DataReader::QMCDataReader();
        this->grids = nullptr;
        this->kernel = nullptr;
        this->annealing_data = nullptr;
        this->annealing_chain = nullptr;
        this->measure = nullptr;
    }

    SAC::~SAC() {
        delete this->grids;
        delete this->kernel;
        delete this->annealing_chain;
        delete this->measure;
        // avoid double free here, no need to delete data
    //    delete this->annealing_data;
    }

    void SAC::set_read_in_params(int lt, double beta, int nbin, int rebin_pace, int bootstrap_num) {
        this->qmc_data_reader->set_params(lt, beta, nbin, rebin_pace, bootstrap_num);
    }

    void SAC::set_filename_tau(const std::string &infile_tau) {
        this->qmc_data_reader->read_tau_from_file(infile_tau);
    }

    void SAC::set_filename_corr(const std::string &infile_corr) {
        this->qmc_data_reader->read_corr_from_file(infile_corr);
    }

    void SAC::set_filename_log(const std::string &outfile_log) {
        this->log_name = outfile_log;
    }

    void SAC::set_griding_params(double freq_interval, double spec_interval, double freq_min, double freq_max) {
        this->grids = new Grids::FreqGrids(freq_interval, spec_interval, freq_min, freq_max);
    }

    void SAC::set_sampling_params(int ndelta, double theta, int max_annealing_steps, int bin_num, int bin_size, int collecting_steps) {
        this->annealing_data = new Annealing::AnnealingData();
        this->annealing_data->ndelta = ndelta;
        this->annealing_data->theta = theta;
        this->collecting_steps = collecting_steps;

        this->annealing_chain = new Annealing::AnnealingChain(max_annealing_steps);
        this->measure = new Measure::Measure(bin_num, bin_size);
    }

    void SAC::set_mode_params(const std::string &kernel_type, const std::string &update_type) {
        assert( update_type == "single" || update_type == "pair" );
        this->update_type = update_type;
        this->kernel_type = kernel_type;
    }


    void SAC::init() {
        // initialize read in module
        this->init_from_module();

        // initialize grids
        this->grids->init();

        // initialize spectrum
        this->init_spectrum();

        // initialize kernel
        this->kernel = new Kernel::Kernel(this->nt, this->grids->FreqNum());
        this->kernel->init(*this, *this->grids, this->kernel_type);
        this->kernel->rotate(qmc_data_reader->rotate_mat);

        // initialize correlations
        this->corr_now.resize(this->nt);
        this->corr_next.resize(this->nt);
        this->compute_corr_from_spec();

        // initialize goodness of fitting (chi2)
        this->chi2 = this->compute_goodness(this->corr_now);
        this->chi2_min = this->chi2;

        this->accept_radio = 0;

        // initialize containers for recovering spectrum
        this->freq.resize(this->grids->SpecNum());
        this->spec.resize(this->grids->SpecNum());

        // free memory
    //    // FIXME
    //    this->qmc_data_reader->deallocate_memory();
    //    delete this->qmc_data_reader;
    //    this->qmc_data_reader = nullptr;
    }

    void SAC::init_from_module() {
        this->qmc_data_reader->analyse_corr();
        this->qmc_data_reader->discard_and_rotate();

        this->nt = this->qmc_data_reader->cov_mat_dim;
        this->beta = this->qmc_data_reader->beta;
        this->scale_factor = this->qmc_data_reader->g0;

        this->tau_from_qmc = this->qmc_data_reader->tau_qmc;
        this->corr_from_qmc = this->qmc_data_reader->rotate_mat * this->qmc_data_reader->corr_mean_qmc;
        this->sigma_from_qmc = (sqrt(this->qmc_data_reader->bootstrap_num) / this->qmc_data_reader->cov_eig.array().sqrt()).matrix();
    }

    void SAC::init_spectrum() {
        // initialize locations of delta functions
        // customized initializing strategies according to different priori information of spectrum

    //    // delta-like distribution
    //    const double delta_freq = -2.0;
    //    this->annealing_data->locations = Eigen::VectorXi::Constant(this->annealing_data->ndelta, this->grids->Freq2GridIndex(delta_freq));

        // rectangle-like distribution
        const double left_edge = -2.0;
        const double right_edge = 2.0;
        const int left_edge_index = this->grids->Freq2FreqIndex(left_edge);
        const int right_edge_index = this->grids->Freq2FreqIndex(right_edge);
        const int rectangle_len = right_edge_index - left_edge_index + 1;
        this->annealing_data->locations.resize(this->annealing_data->ndelta);
        for (int i = 0; i < this->annealing_data->locations.size(); ++i) {
            this->annealing_data->locations(i) = left_edge_index + i % rectangle_len;
        }

    //    // random distribution
    //    std::uniform_int_distribution<> rand_delta(0, this->grids->FreqNum()-1);
    //    this->annealing_data->locations.resize(this->annealing_data->ndelta);
    //    for (int i = 0; i < this->annealing_data->locations.size(); ++i) {
    //        this->annealing_data->locations(i) = rand_delta(Random::Engine);
    //     //    assert( this->annealing_data->locations(i) >=0 );
    //     //    assert( this->annealing_data->locations(i) < this->grids->FreqNum() );
    //    }

    //    // TODO: gaussian-like distribution
    //    const double gaussian_peak = 1.0;
    //    const double gaussian_sigma = 1.0;
    //    // filling from peak to edges

        // initialize amplitudes of delta functions
        // equal amplitudes (scaled)
        // Caution: 1.0 comes from normalization of spectrum
        // TODO: extend to arbitrary normalization factor
        this->annealing_data->amplitude = 1.0 / (this->scale_factor * this->annealing_data->ndelta);

        // initialize width of random move window
        // FIXME: 1/10 of average frequency ?
        double average_freq = std::abs(log(1.0/this->qmc_data_reader->corr_mean_qmc[this->nt-1]) / this->tau_from_qmc[this->nt-1]);
        this->annealing_data->window_width = std::ceil( 0.1 * average_freq / this->grids->FreqInterval() );
        // this->annealing_data->window_width = std::std::ceil( 0.5 * this->grids->FreqNum() );
    }


    void SAC::compute_corr_from_spec() {
        // Prerequisite: Eigen version > 3.3.9
        const Eigen::MatrixXd& tmp_kernel = this->kernel->kernel(Eigen::all, this->annealing_data->locations);
        this->corr_now = tmp_kernel * Eigen::VectorXd::Constant(this->annealing_data->ndelta, this->annealing_data->amplitude);
    }

    double SAC::compute_goodness(const Eigen::VectorXd &corr) const {
        assert( corr.size() == this->nt );
        return ((corr - this->corr_from_qmc).array() * this->sigma_from_qmc.array()).square().sum();
    }

    void SAC::update_deltas_1step() {
        if ( this->update_type == "single" ) {
            this->update_deltas_1step_single();
        }
        else if ( this->update_type == "pair" ) {
            this->update_deltas_1step_pair();
        }
    }

    /**
     *  One Monte Carlo step of updates of delta functions (ndelta number of moving attempt)
     *  randomly move of one delta function for one single attempt
     */
    void SAC::update_deltas_1step_single() {
        // helping params
        std::uniform_int_distribution<> rand_delta(0, this->annealing_data->ndelta-1);
        std::uniform_int_distribution<> rand_width(1, this->annealing_data->window_width);
        std::uniform_int_distribution<> rand_location(0, this->grids->FreqNum()-1);
        int select_delta;
        int move_width;
        int location_now;
        int location_next;
        double chi2_next;
        double p;
        int accept_count = 0;

        // attempt to move for ndelta times
        for (int i = 0; i < this->annealing_data->ndelta; ++i) {
            // randomly select one delta function
            select_delta = rand_delta(Random::Engine);
            location_now = this->annealing_data->locations[select_delta];

            if (this->annealing_data->window_width > 0 && this->annealing_data->window_width < this->grids->FreqNum()) {
                // randomly move within window
                move_width = rand_width(Random::Engine);
                if (std::bernoulli_distribution(0.5)(Random::Engine)) {
                    location_next = location_now + move_width;
                }
                else {
                    location_next = location_now - move_width;
                }

                // abort updates if out of frequency domain
                if (location_next < 0 || location_next >= this->grids->FreqNum()) {
                    i -= 1;
                    continue;
                }
            }
            else if (this->annealing_data->window_width == this->grids->FreqNum()) {
                // randomly move over frequency domain
                location_next = rand_location(Random::Engine);
            }
            else { std::cerr << " Wrong occurs, check the width of updating window ! " << std::endl; exit(1); }

            // compute updated correlation
            this->corr_next = this->corr_now + this->annealing_data->amplitude *
                    ( this->kernel->kernel.col(location_next) - this->kernel->kernel.col(location_now) );

            // compute updated chi2 and accepting radio
            chi2_next = this->compute_goodness(this->corr_next);
            p = exp( (this->chi2 - chi2_next) / (2.0 * this->annealing_data->theta) );

            if ( std::bernoulli_distribution(std::min(p, 1.0))(Random::Engine) ) {
                // accepted
                this->annealing_data->locations[select_delta] = location_next;
                this->corr_now = this->corr_next;
                this->chi2 = chi2_next;
                if ( this->chi2 < this->chi2_min ) {
                    this->chi2_min = this->chi2;
                }
                accept_count++;
            }
        }
        // compute accepting radio
        this->accept_radio = (double)accept_count / this->annealing_data->ndelta;
    }

    /**
     *  One Monte Carlo step of updates of delta functions (ndelta/2 number of moving attempt)
     *  randomly move of two delta functions for one single attempt
     */
    void SAC::update_deltas_1step_pair() {
        assert( this->annealing_data->ndelta >= 2 );

        // helping params
        std::uniform_int_distribution<> rand_delta(0, this->annealing_data->ndelta-1);
        std::uniform_int_distribution<> rand_width(1, this->annealing_data->window_width);
        std::uniform_int_distribution<> rand_location(0, this->grids->FreqNum()-1);
        int select_delta1, select_delta2;
        int move_width1, move_width2;
        int location_now1, location_now2;
        int location_next1, location_next2;
        bool out_of_domain1, out_of_domain2;
        double chi2_next;
        double p;
        int accept_count = 0;

        // attempt to move for ndelta/2 times
        // moving pair of deltas in an attempt
        for (int i = 0; i < std::ceil(this->annealing_data->ndelta/2); i++) {
            // randomly select two different delta functions
            select_delta1 = rand_delta(Random::Engine);
            select_delta2 = select_delta1;
            while ( select_delta1 == select_delta2 ) {
                select_delta2 = rand_delta(Random::Engine);
            }
            location_now1 = this->annealing_data->locations[select_delta1];
            location_now2 = this->annealing_data->locations[select_delta2];

            if (this->annealing_data->window_width >=0 && this->annealing_data->window_width < this->grids->FreqNum()) {
                // randomly select width of moving within window
                move_width1 = rand_width(Random::Engine);
                move_width2 = rand_width(Random::Engine);
                if (std::bernoulli_distribution(0.5)(Random::Engine)) {
                    location_next1 = location_now1 + move_width1;
                    location_next2 = location_now2 - move_width2;
                }
                else {
                    location_next1 = location_now1 - move_width1;
                    location_next2 = location_now2 + move_width2;
                }

                // abort updates if out of frequency domain
                out_of_domain1 = (location_next1 < 0 || location_next1 >= this->grids->FreqNum());
                out_of_domain2 = (location_next2 < 0 || location_next2 >= this->grids->FreqNum());
                if ( out_of_domain1 || out_of_domain2 ) {
                    i -= 1;
                    continue;
                }
            }
            else if (this->annealing_data->window_width == this->grids->FreqNum()) {
                // randomly move over frequency domain
                location_next1 = rand_location(Random::Engine);
                location_next2 = rand_location(Random::Engine);
            }
            else { std::cerr << " Wrong occurs, check the width of updating window ! " << std::endl; exit(1); }

            // compute updated correlation
            this->corr_next = this->corr_now + this->annealing_data->amplitude *
                    ( this->kernel->kernel.col(location_next1) - this->kernel->kernel.col(location_now1)
                    + this->kernel->kernel.col(location_next2) - this->kernel->kernel.col(location_now2) );

            // compute updated chi2 and accepting radio
            chi2_next = this->compute_goodness(this->corr_next);
            p = exp( (this->chi2 - chi2_next) / (2.0 * this->annealing_data->theta) );

            if ( std::bernoulli_distribution(std::min(p, 1.0))(Random::Engine) ) {
                // accepted
                this->annealing_data->locations[select_delta1] = location_next1;
                this->annealing_data->locations[select_delta2] = location_next2;
                this->corr_now = this->corr_next;
                this->chi2 = chi2_next;
                if ( this->chi2 < this->chi2_min ) {
                    this->chi2_min = this->chi2;
                }
                accept_count++;
            }
        }
        // compute accepting radio
        this->accept_radio = (double)accept_count / std::ceil(this->annealing_data->ndelta/2);
    }


    void SAC::update_fixed_theta() {
        // total steps in a fixe theta: nbin * sbin
        for ( int n = 0; n < this->measure->nbin; ++n) {
            // n corresponds to index of bins
            for (int s = 0; s < this->measure->size_of_bin; ++s) {
                // s corresponds to index of samples in one bin
                // recalculate goodness chi2 every 10 steps
                if ( s % 10 == 1 ) {
                    this->chi2 = this->compute_goodness(this->corr_now);
                }
                this->update_deltas_1step();
                this->measure->collect(s, this->chi2, this->accept_radio);
            }

            // compute means for bin of index `n`
            this->measure->bin_analyse(n);

            // writing log
            this->write_log(n);

            // adjust width of window
            // make sure the accepting radio of random move is around 0.5
            if (this->measure->accept_radio_bin(n) > 0.5) {
                this->annealing_data->window_width = std::min((int)std::ceil(this->annealing_data->window_width * 1.5), this->grids->FreqNum());
            }
            if (this->measure->accept_radio_bin(n) < 0.4) {
                this->annealing_data->window_width = std::ceil(this->annealing_data->window_width / 1.5);
            }
        }
    }

    void SAC::write_log(int n) {
        // n labels index of bin number
        // check if log file is empty
        std::ifstream check_empty(this->log_name, std::ios::in|std::ios::app);
        if (!check_empty.is_open()) {
            std::cerr << boost::format(" Fail to open file %s .\n") % this->log_name << std::endl;
            exit(1);
        }
        bool is_empty = (check_empty.peek() == EOF);
        check_empty.close();

        std::ofstream log_out(this->log_name, std::ios::out|std::ios::app);
        if (!log_out.is_open()) {
            std::cerr << boost::format(" Fail to open file %s .\n") % this->log_name << std::endl;
            exit(1);
        }
        if (!is_empty) {
            boost::format log_format("%| 13d|%| 13d|%| 15.3e|%| 15.3e|%| 18.3e|%| 15.3e|%| 15.5f|%| 15.5e|");
            log_out << log_format % (this->annealing_chain->len()+1)
                                  % (n+1) 
                                  % this->annealing_data->theta
                                  % (this->chi2_min/this->nt)
                                  % (this->measure->chi2_bin(n)/this->nt)
                                  % (this->measure->chi2_bin(n)-this->chi2_min)
                                  % this->measure->accept_radio_bin(n)
                                  % (this->annealing_data->window_width * this->grids->FreqInterval())
                    << std::endl;
        }
        else {
            // if empty, print header infomation
            boost::format header_format("%| 13s|%| 13s|%| 15s|%| 15s|%| 18s|%| 15s|%| 15s|%| 15s|");
            log_out << header_format % "AnnealStep" % "BinIndex" % "Theta" % "MinChi2/nt" % "AverageChi2/nt" 
                                     % "DeltaChi2" % "AcceptRadio" % "WindowWidth" << std::endl;
        }
        log_out.close();
    }

    void SAC::perform_annealing() {
        // annealing process, no more than `max_length` steps
        for (int i = 0; i < this->annealing_chain->max_length; ++i) {
            // updating
            this->update_fixed_theta();

            // record simulating information for current theta
            this->measure->analyse();
            this->annealing_data->chi2 = this->measure->chi2();
            this->annealing_chain->push(*this->annealing_data);

            // exit condition
            if (this->measure->chi2() - this->chi2_min < 1e-3) {
                break;
            }

            // lower down sampling temperature
            this->annealing_data->theta /= 1.1;
        }
    }

    void SAC::decide_sampling_theta() {
        // decide sampling temperature by slightly increasing theta
        for (int i = this->annealing_chain->len()-1; i >= 0; --i) {
            // raise chi2 by a standard deviation with respect to the minimum
            if ( this->annealing_chain->chain[i].chi2 > this->chi2_min + 2.0 * sqrt(this->chi2_min) ) {
                this->annealing_data = &this->annealing_chain->chain[i];
                break;
            }
        }
        this->compute_corr_from_spec();
        this->chi2 = this->compute_goodness(this->corr_now);

        // clear up useless information
        this->annealing_chain->clear();
    }


    void SAC::sample_and_collect() {
        // equilibrate at current sampling temperature
        this->update_fixed_theta();

        // generate freq sequence and clear spectrum
        for (int n = 0; n < this->grids->SpecNum(); ++n) {
            this->freq(n) = this->grids->SpecIndex2Freq(n);
        }
        this->spec.setZero();

        // sampling and collecting spectrum
        for (int i = 0; i < this->collecting_steps; ++i) {
            // recompute chi2 every 10 steps
            if ( i % 10 == 1 ) {
                this->chi2 = this->compute_goodness(this->corr_now);
            }
            // updating
            this->update_deltas_1step();

            // collecting
            for (int j = 0; j < this->annealing_data->ndelta; ++j) {
                this->spec( this->grids->FreqIndex2SpecIndex(this->annealing_data->locations(j)) ) += this->annealing_data->amplitude;
            }
        }

        // scaling and recovering spectrum
        this->spec *=  this->scale_factor / (this->collecting_steps * this->grids->SpecInterval());
    }

    void SAC::output(const std::string &filename) {
        std::ofstream outfile(filename, std::ios::out|std::ios::trunc);
        if (!outfile.is_open()) {
            std::cerr << boost::format(" Fail to open file %s .\n") % filename << std::endl;
            exit(1);
        }
        outfile << std::setiosflags(std::ios::right);
            for (int i = 0; i < this->grids->SpecNum(); ++i) {
                outfile << std::setw(10) << i
                        << std::setw(15) << this->freq(i)
                        << std::setw(15) << this->spec(i)
                        << std::endl;
            }
        outfile.close();
    }


} // namespace Simualtion