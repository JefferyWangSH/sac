#include "measure.h"

namespace Measure {
    
    Measure::Measure(int nbin, int size_of_bin) {
        this->nbin = nbin;
        this->size_of_bin = size_of_bin;

        this->chi2_sample.resize(size_of_bin);
        this->accept_radio_sample.resize(size_of_bin);
        this->chi2_bin.resize(nbin);
        this->accept_radio_bin.resize(nbin);
    }

    void Measure::resize(int nbin, int size_of_bin) {
        this->nbin = nbin;
        this->size_of_bin = size_of_bin;
        this->accept_radio_sample.resize(size_of_bin);
        this->chi2_sample.resize(size_of_bin);
        this->chi2_bin.resize(nbin);
        this->accept_radio_bin.resize(nbin);
        this->clear();
    }

    void Measure::clear() {
        // clear previous data
        this->chi2_mean = 0.0;
        this->chi2_err = 0.0;
        this->accept_radio_mean = 0.0;
        this->accept_radio_err = 0.0;

        this->chi2_sample.setZero();
        this->accept_radio_sample.setZero();
        this->chi2_bin.setZero();
        this->accept_radio_bin.setZero();
    }

    void Measure::collect(int s, double chi2, double accept_radio) {
        assert( s >= 0 && s < this->size_of_bin );
        // s labels index of samples within certain bin
        this->accept_radio_sample(s) = accept_radio;
        this->chi2_sample(s) = chi2;
    }

    void Measure::bin_analyse(int n) {
        // compute mean value for one bin of samples
        // n labels index of certain one bin
        this->chi2_bin(n) = this->chi2_sample.sum() / this->size_of_bin;
        this->accept_radio_bin(n) = this->accept_radio_sample.sum() / this->size_of_bin;
    }

    void Measure::analyse() {
        // compute mean values and estimate errors for all bins
        this->chi2_mean = this->chi2_bin.sum() / this->nbin;
        this->accept_radio_mean = this->accept_radio_bin.sum() / this->nbin;
        this->chi2_err = sqrt( (this->chi2_bin.array() - this->chi2_mean).square().sum() / (this->nbin - 1) );
        this->accept_radio_err = sqrt( (this->accept_radio_bin.array() - this->accept_radio_mean).square().sum() / (this->nbin - 1) );
    }

    double Measure::chi2() const{
        return this->chi2_mean;
    }

    double Measure::accept_radio() const {
        return this->accept_radio_mean;
    }

} // namespace Measure
