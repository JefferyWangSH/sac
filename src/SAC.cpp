#include "SAC.h"


void SAC::set_read_in_params(int _lt, double _beta, int _nbin, int _rebin_pace, int _num_bootstrap) {
    this->readin.set_params(_lt, _beta, _nbin, _rebin_pace, _num_bootstrap);
}

void SAC::set_filename_tau(const std::string &infile_tau) {
    this->readin.read_tau_from_file(infile_tau);
}

void SAC::set_filename_corr(const std::string &infile_corr) {
    this->readin.read_corr_from_file(infile_corr);
}

void SAC::set_griding_params(double _spec_delta_omega, double _fine_delta_omega, double _omega_min, double _omega_max) {
    this->spec_delta_omega = _spec_delta_omega;
    this->fine_delta_omega = _fine_delta_omega;
    this->omega_min = _omega_min;
    this->omega_max = _omega_max;
}

void SAC::set_sampling_params(double _num_delta) {
    this->num_delta = _num_delta;
}

void SAC::init() {
    // initialize read in module
    readin.analyse_corr();
    readin.discard_and_rotate();
    readin.deallocate_memory();
}
