#include "SAC.h"

#include <iostream>


Simulation::SAC::SAC() {
    this->readin = new QMCData::ReadInModule();
}

Simulation::SAC::~SAC() {
    delete this->grid;
    delete this->kernel;
    delete this->delta;
    delete this->anneal;
    std::cout << "simulation was done！ :P" << std::endl;
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

void Simulation::SAC::set_sampling_params(double ndelta, double theta, int max_annealing_steps) {
    this->delta = new Annealing::DeltaData();
    this->delta->ndelta = ndelta;
    this->delta->theta = theta;

    this->anneal = new Annealing::AnnealChain(max_annealing_steps);
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
    this->kernel = new Kernel::Kernel(nt, grid->upper());
    this->kernel->init(*this, *grid, kernel_mode);
    this->kernel->rotate(readin->rotate_mat);

    // initialize correlations
    this->corr_current.resize(nt);
    this->corr_update.resize(nt);
    this->compute_corr_from_spec();

    // initialize goodness of fitting (chi2)
    this->chi2 = compute_goodness(this->corr_current);
    this->chi2_minimum = this->chi2;

    this->accept_count = 0;
    this->accept_radio = 0.0;

    // free memory
    this->readin->deallocate_memory();
    delete this->readin;
}

void Simulation::SAC::init_from_module() {
    this->readin->analyse_corr();
    this->readin->discard_and_rotate();

    this->nt = readin->cov_mat_dim;
    this->beta = readin->beta;
    this->scale_factor = readin->g0;
    this->nbootstrap = readin->num_bootstrap;

    this->tau = readin->tau;
    this->corr = readin->rotate_mat * readin->corr_mean;
    this->sigma = (sqrt(readin->num_bootstrap) / readin->cov_eig.array().sqrt()).matrix();
}

void Simulation::SAC::init_spectrum() {
    // initialize locations of delta functions
    // for symmetric spectrum, initialize locations near middle of frequency domain
    this->delta->locations = Eigen::VectorXi::Constant(delta->ndelta, ceil(0.5 * (grid->lower() + grid->upper())));

    // equal amplitudes
    this->delta->amplitude = 1.0 / delta->ndelta;

    // width of random move window
    // FIXME: 1/10 of average frequency ?
    const double average_freq = log(1.0/readin->corr_mean[nt-1]) / tau[nt-1];
    this->delta->window_width = ceil( 0.1 * average_freq / grid->interval() );
}

void Simulation::SAC::compute_corr_from_spec() {
    // Prerequisite: Eigen version > 3.3.9
    const Eigen::MatrixXd tmp_kernel = kernel->kernel(Eigen::VectorXi::LinSpaced(nt, 0, nt), delta->locations);
    corr_current = tmp_kernel * Eigen::VectorXd::Constant(delta->ndelta, delta->amplitude);
}

double Simulation::SAC::compute_goodness(const Eigen::VectorXd &corr_from_spectrum) {
    assert( corr_from_spectrum.size() == nt );
    return ((corr_from_spectrum - corr).array() * sigma.array()).square().sum();
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
    std::uniform_int_distribution<> rand_delta(0, delta->ndelta-1);
    std::uniform_int_distribution<> rand_width(1, delta->window_width);
    int select_delta;
    int select_width;
    int location_current;
    int location_updated;
    double chi2_updated;
    double p;

    // attempt to move for ndelta times
    this->accept_count = 0;
    for (int i = 0; i < delta->ndelta; ++i) {
        // randomly select one delta function
        select_delta = rand_delta(rand_engine_sac);
        select_width = rand_width(rand_engine_sac);
        location_current = delta->locations[select_delta];

        if (std::bernoulli_distribution(0.5)(rand_engine_sac)) {
            location_updated = location_current + select_width;
        }
        else {
            location_updated = location_current - select_width;
        }

        // abort updates if out of frequency domain
        if (location_updated < grid->lower() || location_updated >= grid->upper()) {
//             FIXME : check potential differences here
            i -= 1;
            continue;
        }

        // compute updated correlation
        this->corr_update = this->corr_current + delta->amplitude *
                ( kernel->kernel.col(location_updated) - kernel->kernel.col(location_current) );

        // compute updated chi2 and accepting radio
        chi2_updated = compute_goodness(this->corr_update);
        p = exp( (this->chi2 - chi2_updated) / (2.0 * delta->theta) );

        if ( std::bernoulli_distribution(std::min(p, 1.0))(rand_engine_sac) ) {
            // accepted
            this->delta->locations[select_delta] = location_updated;
            this->corr_current = this->corr_update;
            this->chi2 = chi2_updated;
            if ( this->chi2 < this->chi2_minimum ) {
                this->chi2_minimum = this->chi2;
            }
            this->accept_count++;
        }
    }

    // compute accepting radio
    this->accept_radio = accept_count / delta->ndelta;
    this->accept_count = 0;
}

/**
  *  One Monte Carlo step of updates of delta functions (ndelta/2 number of moving attempt)
  *  randomly move of two delta functions for one single attempt
  */
void Simulation::SAC::update_deltas_1step_pair() {
    assert( delta->ndelta >= 2 );

    // helping params
    std::uniform_int_distribution<> rand_delta(0, delta->ndelta-1);
    std::uniform_int_distribution<> rand_width(1, delta->window_width);



    // attempt to move for ndelta/2 times
    // moving pair of deltas in an attempt
    this->accept_count = 0;
    for (int i = 0; i < delta->ndelta/2; i++) {
        // randomly select two delta functions

    }


//    // helping params
//    std::uniform_int_distribution<> rand_delta(0, delta->ndelta-1);
//    std::uniform_int_distribution<> rand_width(1, delta->window_width);
//    int select_delta;
//    int select_width;
//    int location_current;
//    int location_updated;
//    double chi2_updated;
//    double p;
//
//    // attempt to move for ndelta times
//    this->accept_count = 0;
//    for (int i = 0; i < delta->ndelta; ++i) {
//        select_delta = rand_delta(rand_engine_sac);
//        select_width = rand_width(rand_engine_sac);
//        location_current = delta->locations[select_delta];
//
//        if (std::bernoulli_distribution(0.5)(rand_engine_sac)) {
//            location_updated = location_current + select_width;
//        }
//        else {
//            location_updated = location_current - select_width;
//        }
//
//        // abort updates if out of frequency domain
//        if (location_updated < grid->lower() || location_updated >= grid->upper()) {
//            //             FIXME : check potential differences here
//            i -= 1;
//            continue;
//        }
//
//        // compute updated correlation
//        this->corr_update = this->corr_current + delta->amplitude *
//                ( kernel->kernel.col(location_updated) - kernel->kernel.col(location_current) );
//
//        // compute updated chi2 and accepting radio
//        chi2_updated = compute_goodness(this->corr_update);
//        p = exp( (this->chi2 - chi2_updated) / (2.0 * delta->theta) );
//
//        if ( std::bernoulli_distribution(std::min(p, 1.0))(rand_engine_sac) ) {
//            // accepted
//            this->delta->locations[select_delta] = location_updated;
//            this->corr_current = this->corr_update;
//            this->chi2 = chi2_updated;
//            if ( this->chi2 < this->chi2_minimum ) {
//                this->chi2_minimum = this->chi2;
//            }
//            this->accept_count++;
//        }
//    }
//
//    // compute accepting radio
//    this->accept_radio = accept_count / delta->ndelta;
//    this->accept_count = 0;
}



