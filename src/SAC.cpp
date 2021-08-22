#include "SAC.h"


SAC::SAC() {
    this->readin = new ReadInModule();
}

void SAC::set_read_in_params(int _lt, double _beta, int _nbin, int _rebin_pace, int _num_bootstrap) {
    this->readin->set_params(_lt, _beta, _nbin, _rebin_pace, _num_bootstrap);
}

void SAC::set_filename_tau(const std::string &infile_tau) {
    this->readin->read_tau_from_file(infile_tau);
}

void SAC::set_filename_corr(const std::string &infile_corr) {
    this->readin->read_corr_from_file(infile_corr);
}

void SAC::set_griding_params(double _grid_interval, double _spec_interval, double _omega_min, double _omega_max) {
    this->grid.set_params(_grid_interval, _spec_interval, _omega_min, _omega_max);
}

void SAC::set_sampling_params(double _ndelta, double _theta) {
    this->ndelta = _ndelta;
    this->theta = _theta;
}

void SAC::init() {
    // initialize read in module
    this->init_from_module();

    // initialize grids
    this->grid.init();

    // initialize spectrum
    this->init_spectrum();

    // initialize kernel
    this->kernel.set_params(nt, grid.upper());
    this->kernel.init(*this, grid, "fermion");
    this->kernel.rotate(readin->rotate_mat);

    // free memory
    this->readin->deallocate_memory();
    delete this->readin;
}

void SAC::init_from_module() {
    this->readin->analyse_corr();
    this->readin->discard_and_rotate();

    this->nt = readin->cov_mat_dim;
    this->beta = readin->beta;
    this->scale_factor = readin->g0;
    this->nbootstrap = readin->num_bootstrap;

    this->tau = readin->tau_seq;
    this->corr = readin->rotate_mat * readin->corr_mean_seq;
    this->sigma = (sqrt(readin->num_bootstrap) / readin->cov_eig.array().sqrt()).matrix();
}

void SAC::init_spectrum() {

    // initialize locations of delta functions
    // for symmetric spectrum, initialize locations near middle of frequency domain
    this->delta_locations.reserve(ndelta);
    for (int i = 0; i < ndelta; ++i) {
        this->delta_locations.emplace_back( ceil(0.5 * (grid.lower() + grid.upper())) );
    }

    // equal amplitudes
    this->delta_amplitude = 1.0 / ndelta;

    // width of random move window
    // FIXME: 1/10 of average frequency ?
    const double average_freq = log(1.0/readin->corr_mean_seq[nt-1]) / tau[nt-1];
    this->window_width = ceil( 0.1 * average_freq / grid.interval() );
}
