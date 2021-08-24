#include "Measure.h"

Measure::Measure::Measure(int nbin, int sbin) {
    this->nbin = nbin;
    this->sbin = sbin;

    this->sample_accept_radio.resize(sbin);
    this->sample_chi2.resize(sbin);
    this->bin_chi2.resize(nbin);
    this->bin_accept_radio.resize(nbin);
}

void Measure::Measure::clear() {
    // clear previous data
    this->sample_accept_radio.setZero();
    this->sample_chi2.setZero();
    this->bin_chi2.setZero();
    this->bin_accept_radio.setZero();
}

void Measure::Measure::fill(int s, double chi2, double accept_radio) {
    assert( s >= 0 && s < sbin );
    // s corresponds to index of samples in one bin
    this->sample_accept_radio(s) = accept_radio;
    this->sample_chi2(s) = chi2;
}

void Measure::Measure::bin_analyse(int n) {
    // compute means for one bin
    // n corresponds to index of certain one bin
    this->bin_chi2(n) = this->sample_chi2.sum() / this->sbin;
    this->bin_accept_radio(n) = this->sample_accept_radio.sum() / this->sbin;
}

void Measure::Measure::analyse() {
    // computing means and errors of all bins
    this->chi2_mean = this->bin_chi2.sum() / this->nbin;
    this->accept_radio_mean = this->bin_accept_radio.sum() / this->nbin;
    this->chi2_err = sqrt( (this->bin_chi2.array() - this->chi2_mean).square().sum() / (this->nbin - 1) );
    this->accept_radio_err = sqrt( (this->bin_accept_radio.array() - this->accept_radio_mean).square().sum() / (this->nbin - 1) );
}

const double Measure::Measure::chi2() const{
    return this->chi2_mean;
}

const double Measure::Measure::accept_radio() const {
    return this->accept_radio_mean;
}


