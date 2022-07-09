#include "freq_grids.h"
#include <cassert>
#include <cmath>

namespace Grids {

    double FreqGrids::FreqInterval() const { return this->m_freq_interval; }
    double FreqGrids::SpecInterval() const { return this->m_spec_interval; }
    int FreqGrids::FreqNum() const { return this->m_num_freq; }
    int FreqGrids::SpecNum() const { return this->m_num_spec; }


    FreqGrids::FreqGrids( double freq_interval, double spec_interval, double freq_min, double freq_max ) 
    {
        this->m_freq_interval = freq_interval;
        this->m_spec_interval = spec_interval;
        this->m_freq_min = freq_min;
        this->m_freq_max = freq_max;
    }


    void FreqGrids::initial() 
    {
        // generate both hyperfine and spectral frequency grids
        this->m_int_freq_min = 0;
        this->m_int_freq_max = std::ceil((this->m_freq_max - this->m_freq_min) / this->m_freq_interval);
        this->m_num_freq = this->m_int_freq_max;
        this->m_num_spec = std::ceil((this->m_freq_max - this->m_freq_min) / this->m_spec_interval);
    }


    double FreqGrids::FreqIndex2Freq( int freq_index ) const 
    {
        assert( freq_index >= 0 );
        assert( freq_index < this->m_num_freq );
        return this->m_freq_min + freq_index * this->m_freq_interval;
    }


    int FreqGrids::Freq2FreqIndex( double freq ) const 
    {
        assert( freq >= this->m_freq_min && freq < this->m_freq_max );
        int index = std::ceil((freq - this->m_freq_min) / this->m_freq_interval);
        assert( index >= this->m_int_freq_min && index < this->m_int_freq_max );
        return index;
    }


    double FreqGrids::SpecIndex2Freq( int spec_index ) const 
    {
        assert( spec_index >= 0 );
        assert( spec_index < this->m_num_spec );
        return this->m_freq_min + spec_index * this->m_spec_interval;
    }


    int FreqGrids::FreqIndex2SpecIndex( int freq_index ) const 
    {
        assert( freq_index >= 0 );
        assert( freq_index < this->m_num_freq );
        return std::floor(freq_index * this->m_freq_interval / this->m_spec_interval);
    }

} // namespace Grids
