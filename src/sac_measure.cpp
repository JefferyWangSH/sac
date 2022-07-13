#include "sac_measure.h"

namespace SAC {

    // interface member functions
    int Measure::number_of_bin() const { return this->m_number_of_bin; }
    int Measure::size_of_bin() const { return this->m_size_of_bin; }
    double Measure::chi2() const{ return this->m_chi2_mean; }
    double Measure::chi2_error() const { return this->m_chi2_err; }
    double Measure::accept_ratio() const { return this->m_accept_ratio_mean; }
    double Measure::accept_ratio_error() const { return this->m_accept_ratio_err; }

    double Measure::chi2( int n ) const 
    {
        assert( n >= 0 && n < this->m_number_of_bin );
        return this->m_chi2_bin(n);
    }

    double Measure::accept_ratio( int n ) const 
    {
        assert( n >= 0 && n < this->m_number_of_bin );
        return this->m_accept_ratio_bin(n);
    }

    
    Measure::Measure( int number_of_bin, int size_of_bin ) 
    {
        this->m_number_of_bin = number_of_bin;
        this->m_size_of_bin = size_of_bin;

        this->m_chi2_sample.resize(size_of_bin);
        this->m_accept_ratio_sample.resize(size_of_bin);
        this->m_chi2_bin.resize(number_of_bin);
        this->m_accept_ratio_bin.resize(number_of_bin);
    }


    void Measure::resize( int number_of_bin, int size_of_bin ) 
    {
        this->m_number_of_bin = number_of_bin;
        this->m_size_of_bin = size_of_bin;

        this->m_accept_ratio_sample.resize(size_of_bin);
        this->m_chi2_sample.resize(size_of_bin);
        this->m_chi2_bin.resize(number_of_bin);
        this->m_accept_ratio_bin.resize(number_of_bin);
        
        this->clear();
    }


    void Measure::clear() 
    {
        // clear previous data
        this->m_chi2_mean = 0.0;
        this->m_chi2_err = 0.0;
        this->m_accept_ratio_mean = 0.0;
        this->m_accept_ratio_err = 0.0;

        this->m_chi2_sample.setZero();
        this->m_accept_ratio_sample.setZero();
        this->m_chi2_bin.setZero();
        this->m_accept_ratio_bin.setZero();
    }


    void Measure::collect( int s, double chi2, double accept_ratio )
    {
        assert( s >= 0 && s < this->m_size_of_bin );
        // s labels the index of samples
        this->m_accept_ratio_sample(s) = accept_ratio;
        this->m_chi2_sample(s) = chi2;
    }


    void Measure::bin_analyse( int n )
    {
        assert( n >= 0 && n < this->m_number_of_bin );
        // compute the mean value for one bin of samples
        // n labels the index of bin
        this->m_chi2_bin(n) = this->m_chi2_sample.sum() / this->m_size_of_bin;
        this->m_accept_ratio_bin(n) = this->m_accept_ratio_sample.sum() / this->m_size_of_bin;
    }


    void Measure::analyse()
    {
        // compute mean values and estimate errors for all bins
        this->m_chi2_mean = this->m_chi2_bin.sum() / this->m_number_of_bin;
        this->m_accept_ratio_mean = this->m_accept_ratio_bin.sum() / this->m_number_of_bin;
        this->m_chi2_err = std::sqrt( (this->m_chi2_bin.array() - this->m_chi2_mean).square().sum() / (this->m_number_of_bin - 1) );
        this->m_accept_ratio_err = std::sqrt( (this->m_accept_ratio_bin.array() - this->m_accept_ratio_mean).square().sum() / (this->m_number_of_bin - 1) );
    }

} // namespace SAC
