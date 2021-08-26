#include "SAC.h"

#include <iostream>
#include <fstream>
#include <iomanip>


Simulation::SAC::SAC() {
    this->readin = new QMCData::ReadInModule();
    this->grid = nullptr;
    this->kernel = nullptr;
    this->data = nullptr;
    this->anneal = nullptr;
    this->measure = nullptr;
}


Simulation::SAC::~SAC() {
    delete this->grid;
    delete this->kernel;
    delete this->anneal;
    delete this->measure;
    // avoid double free here, no need to delete data
    //    delete this->data;
}


void Simulation::SAC::set_read_in_params(int _lt, double _beta, int _nbin, int _rebin_pace, int _num_bootstrap) {
    this->readin->set_params(_lt, _beta, _nbin, _rebin_pace, _num_bootstrap);
}


void Simulation::SAC::set_filename_tau(const std::string &infile_tau) {
    this->readin->read_tau_from_file(infile_tau);
}


void Simulation::SAC::set_filename_corr(const std::string &infile_corr) {
    this->readin->read_corr_from_file(infile_corr);
}


void Simulation::SAC::set_griding_params(double grid_interval, double spec_interval, double omega_min, double omega_max) {
    this->grid = new Grid::FrequencyGrid(grid_interval, spec_interval, omega_min, omega_max);
}


void Simulation::SAC::set_sampling_params(double ndelta, double theta, int max_annealing_steps, int _bin_num, int _bin_size, int _collecting_steps) {
    this->data = new Annealing::AnnealData();
    this->data->ndelta = ndelta;
    this->data->theta = theta;
    this->collecting_steps = _collecting_steps;

    this->anneal = new Annealing::AnnealChain(max_annealing_steps);
    this->measure = new Measure::Measure(_bin_num, _bin_size);
}


void Simulation::SAC::set_mode_params(const std::string &_kernel_mode, const std::string &_update_mode) {
    assert( _update_mode == "single" || _update_mode == "pair" );
    this->update_mode = _update_mode;
    this->kernel_mode = _kernel_mode;
}


void Simulation::SAC::init() {
    // initialize read in module
    this->init_from_module();

    // initialize grids
    this->grid->init();

    // initialize spectrum
    this->init_spectrum();

    // initialize kernel
    this->kernel = new Kernel::Kernel(this->nt, this->grid->GridsNum());
    this->kernel->init(*this, *this->grid, this->kernel_mode);
    this->kernel->rotate(readin->rotate_mat);

    // initialize correlations
    this->corr_current.resize(this->nt);
    this->corr_update.resize(this->nt);
    this->compute_corr_from_spec();

    // initialize goodness of fitting (chi2)
    this->chi2 = this->compute_goodness(this->corr_current);
    this->chi2_minimum = this->chi2;

    this->accept_radio = 0;

    // initialize containers for recovering spectrum
    this->freq.resize(this->grid->SpecNum());
    this->spectrum.resize(this->grid->SpecNum());

    // free memory
    this->readin->deallocate_memory();
    delete this->readin;
    this->readin = nullptr;
}


void Simulation::SAC::init_from_module() {
    this->readin->analyse_corr();
    this->readin->discard_and_rotate();

    this->nt = this->readin->cov_mat_dim;
    this->beta = this->readin->beta;
    this->scale_factor = this->readin->g0;

    this->tau = this->readin->tau;
    this->corr = this->readin->rotate_mat * this->readin->corr_mean;
    this->sigma = (sqrt(this->readin->num_bootstrap) / this->readin->cov_eig.array().sqrt()).matrix();
}


void Simulation::SAC::init_spectrum() {
    // initialize locations of delta functions
    // for symmetric spectrum, initialize locations near middle of frequency domain
    this->data->locations = Eigen::VectorXi::Constant(this->data->ndelta, ceil(0.5 * (0 + this->grid->GridsNum())));

    // equal amplitudes
    this->data->amplitude = 1.0 / this->data->ndelta;

    // width of random move window
    // FIXME: 1/10 of average frequency ?
    const double average_freq = log(1.0/this->readin->corr_mean[nt-1]) / this->tau[nt-1];
    this->data->window_width = ceil( 0.1 * average_freq / this->grid->GridInterval() );
}


void Simulation::SAC::compute_corr_from_spec() {
    // Prerequisite: Eigen version > 3.3.9
    const Eigen::MatrixXd tmp_kernel = this->kernel->kernel(Eigen::VectorXi::LinSpaced(this->nt, 0, this->nt), this->data->locations);
    this->corr_current = tmp_kernel * Eigen::VectorXd::Constant(this->data->ndelta, this->data->amplitude);
}


const double Simulation::SAC::compute_goodness(const Eigen::VectorXd &corr_from_spectrum) const{
    assert( corr_from_spectrum.size() == this->nt );
    return ((corr_from_spectrum - this->corr).array() * this->sigma.array()).square().sum();
}


void Simulation::SAC::update_deltas_1step() {
    if ( this->update_mode == "single" ) {
        this->update_deltas_1step_single();
    }
    else if ( this->update_mode == "pair" ) {
        this->update_deltas_1step_pair();
    }
}


/**
  *  One Monte Carlo step of updates of delta functions (ndelta number of moving attempt)
  *  randomly move of one delta function for one single attempt
  */
void Simulation::SAC::update_deltas_1step_single() {
    // helping params
    std::uniform_int_distribution<> rand_delta(0, this->data->ndelta-1);
    std::uniform_int_distribution<> rand_width(1, this->data->window_width);
    int select_delta;
    int move_width;
    int location_current;
    int location_updated;
    double chi2_updated;
    double p;
    int accept_count = 0;

    // attempt to move for ndelta times
    for (int i = 0; i < this->data->ndelta; ++i) {
        // randomly select one delta function
        select_delta = rand_delta(rand_engine_sac);
        move_width = rand_width(rand_engine_sac);
        location_current = this->data->locations[select_delta];

        if (std::bernoulli_distribution(0.5)(rand_engine_sac)) {
            location_updated = location_current + move_width;
        }
        else {
            location_updated = location_current - move_width;
        }

        // abort updates if out of frequency domain
        if (location_updated < 0 || location_updated >= this->grid->GridsNum()) {
            i -= 1;
            continue;
        }

        // compute updated correlation
        this->corr_update = this->corr_current + this->data->amplitude *
                ( this->kernel->kernel.col(location_updated) - this->kernel->kernel.col(location_current) );

        // compute updated chi2 and accepting radio
        chi2_updated = this->compute_goodness(this->corr_update);
        p = exp( (this->chi2 - chi2_updated) / (2.0 * this->data->theta) );

        if ( std::bernoulli_distribution(std::min(p, 1.0))(rand_engine_sac) ) {
            // accepted
            this->data->locations[select_delta] = location_updated;
            this->corr_current = this->corr_update;
            this->chi2 = chi2_updated;
            if ( this->chi2 < this->chi2_minimum ) {
                this->chi2_minimum = this->chi2;
            }
            accept_count++;
        }
    }
    // compute accepting radio
    this->accept_radio = (double)accept_count / this->data->ndelta;
}


/**
  *  One Monte Carlo step of updates of delta functions (ndelta/2 number of moving attempt)
  *  randomly move of two delta functions for one single attempt
  */
void Simulation::SAC::update_deltas_1step_pair() {
    assert( this->data->ndelta >= 2 );

    // helping params
    std::uniform_int_distribution<> rand_delta(0, this->data->ndelta-1);
    std::uniform_int_distribution<> rand_width(1, this->data->window_width);
    int select_delta1, select_delta2;
    int move_width1, move_width2;
    int location_current1, location_current2;
    int location_updated1, location_updated2;
    bool out_of_domain1, out_of_domain2;
    double chi2_updated;
    double p;
    int accept_count = 0;

    // attempt to move for ndelta/2 times
    // moving pair of deltas in an attempt
    for (int i = 0; i < ceil(this->data->ndelta/2); i++) {
        // randomly select two different delta functions
        select_delta1 = rand_delta(rand_engine_sac);
        select_delta2 = select_delta1;
        while ( select_delta1 == select_delta2 ) {
            select_delta2 = rand_delta(rand_engine_sac);
        }

        // randomly select width of moving within window
        move_width1 = rand_width(rand_engine_sac);
        move_width2 = rand_width(rand_engine_sac);
        location_current1 = this->data->locations[select_delta1];
        location_current2 = this->data->locations[select_delta2];

        if (std::bernoulli_distribution(0.5)(rand_engine_sac)) {
            location_updated1 = location_current1 + move_width1;
            location_updated2 = location_current2 - move_width2;
        }
        else {
            location_updated1 = location_current1 - move_width1;
            location_updated2 = location_current2 + move_width2;
        }

        // abort updates if out of frequency domain
        out_of_domain1 = (location_updated1 < 0 || location_updated1 >= this->grid->GridsNum());
        out_of_domain2 = (location_updated2 < 0 || location_updated2 >= this->grid->GridsNum());
        if ( out_of_domain1 || out_of_domain2 ) {
            i -= 1;
            continue;
        }

        // compute updated correlation
        this->corr_update = this->corr_current + this->data->amplitude *
                ( this->kernel->kernel.col(location_updated1) - this->kernel->kernel.col(location_current1)
                + this->kernel->kernel.col(location_updated2) - this->kernel->kernel.col(location_current2) );

        // compute updated chi2 and accepting radio
        chi2_updated = this->compute_goodness(this->corr_update);
        p = exp( (this->chi2 - chi2_updated) / (2.0 * this->data->theta) );

        if ( std::bernoulli_distribution(std::min(p, 1.0))(rand_engine_sac) ) {
            // accepted
            this->data->locations[select_delta1] = location_updated1;
            this->data->locations[select_delta2] = location_updated2;
            this->corr_current = this->corr_update;
            this->chi2 = chi2_updated;
            if ( this->chi2 < this->chi2_minimum ) {
                this->chi2_minimum = this->chi2;
            }
            accept_count++;
        }
    }
    // compute accepting radio
    this->accept_radio = (double)accept_count / ceil(this->data->ndelta/2);
}


void Simulation::SAC::update_fixed_theta() {
    // total steps in a fixe theta: nbin * sbin
    for ( int n = 0; n < this->measure->nbin; ++n) {
        // n corresponds to index of bins
        for (int s = 0; s < this->measure->sbin; ++s) {
            // s corresponds to index of samples in one bin
            // recalculate goodness chi2 every 10 steps
            if ( s % 10 == 1 ) {
                this->chi2 = this->compute_goodness(this->corr_current);
            }
            this->update_deltas_1step();
            this->measure->fill(s, this->chi2, this->accept_radio);
        }

        // compute means for bin of index `n`
        this->measure->bin_analyse(n);

        // writing log
        this->write_log(n);

        // adjust width of window
        // make sure the accepting radio of random move is around 0.5
        if (this->measure->bin_accept_radio(n) > 0.5) {
            this->data->window_width = (ceil(this->data->window_width * 1.5) < this->grid->GridsNum())?
                                        ceil(this->data->window_width * 1.5) : this->grid->GridsNum();
        }
        if (this->measure->bin_accept_radio(n) < 0.4) {
            this->data->window_width = ceil(this->data->window_width / 1.5);
        }
    }
}


void Simulation::SAC::write_log(int n) {
    // n labels index of bin number
    std::ofstream log;
    log.open("../results/log.log", std::ios::out|std::ios::app);
    log << std::setiosflags(std::ios::right)
        << std::setw(10) << this->anneal->len() + 1
        << std::setw(10) << n + 1
        << std::setw(15) << this->data->theta
        << std::setw(15) << this->chi2_minimum / this->nt
        << std::setw(15) << this->measure->bin_chi2(n) / this->nt
        << std::setw(15) << this->measure->bin_chi2(n) - this->chi2_minimum
        << std::setw(15) << this->measure->bin_accept_radio(n)
        << std::setw(15) << this->data->window_width * this->grid->GridInterval()
        << std::endl;
    log.close();
}


void Simulation::SAC::perform_annealing() {
    // annealing process, no more than `max_length` steps
    for (int i = 0; i < this->anneal->max_length; ++i) {
        // updating
        this->update_fixed_theta();

        // record simulating information for current theta
        this->measure->analyse();
        this->data->chi2 = this->measure->chi2();
        this->anneal->push(*this->data);

        // exit condition
        if (this->measure->chi2() - this->chi2_minimum < 1e-3) {
            break;
        }

        // lower down sampling temperature
        this->data->theta /= 1.1;
    }
}


void Simulation::SAC::decide_sampling_theta() {
    // decide sampling temperature by slightly increasing theta
    for (int i = this->anneal->len()-1; i >= 0; --i) {
        // raise chi2 by a standard deviation with respect to the minimum
        if ( this->anneal->chain[i].chi2 > this->chi2_minimum + 2.0 * sqrt(this->chi2_minimum) ) {
            this->data = &this->anneal->chain[i];
            break;
        }
    }

    // clear up useless information
    this->anneal->clear();
}


void Simulation::SAC::sample_and_collect() {
    // equilibrate at current sampling temperature
    this->update_fixed_theta();

    // generate freq sequence and clear spectrum
    for (int n = 0; n < this->grid->SpecNum(); ++n) {
        this->freq(n) = this->grid->SpecIndex2Freq(n);
    }
    this->spectrum.setZero();

    // sampling and collecting spectrum
    for (int i = 0; i < this->collecting_steps; ++i) {
        // recompute chi2 every 10 steps
        if ( i % 10 == 1 ) {
            this->chi2 = this->compute_goodness(this->corr_current);
        }
        // updating
        this->update_deltas_1step();

        // collecting
        for (int j = 0; j < this->data->ndelta; ++j) {
            this->spectrum( this->grid->Grid2Spec(this->data->locations(j)) ) += this->data->amplitude;
        }
    }

    // scaling and recovering spectrum
    this->spectrum *=  this->scale_factor / (this->collecting_steps * this->grid->SpecInterval());
}


void Simulation::SAC::output(const std::string &filename) {
    std::ofstream outfile;
    outfile.open(filename, std::ios::out|std::ios::trunc);
    outfile << std::setiosflags(std::ios::right);

    for (int i = 0; i < this->grid->SpecNum(); ++i) {
        outfile << std::setw(10) << i
                << std::setw(15) << this->freq(i)
                << std::setw(15) << this->spectrum(i)
                << std::endl;
    }
    outfile.close();
}
