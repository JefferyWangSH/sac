#include "freq_grids.h"
#include <cassert>
#include <cmath>

namespace Grids {

    FreqGrids::FreqGrids(double freq_interval, double spec_interval, double freq_min, double freq_max) {
        this->freq_interval = freq_interval;
        this->spec_interval = spec_interval;
        this->freq_min = freq_min;
        this->freq_max = freq_max;
    }

    void FreqGrids::init() {
        // express frequency in unit of griding interval
        this->int_freq_min = 0;
        this->int_freq_max = std::ceil((this->freq_max - this->freq_min) / this->freq_interval);
        this->num_freq = this->int_freq_max;
        this->num_spec = std::ceil((this->freq_max - this->freq_min) / this->spec_interval);
    }

    double FreqGrids::FreqInterval() const {
        return this->freq_interval;
    }

    double FreqGrids::SpecInterval() const {
        return this->spec_interval;
    }

    int FreqGrids::FreqNum() const {
        return this->num_freq;
    }

    int FreqGrids::SpecNum() const {
        return this->num_spec;
    }

    double FreqGrids::FreqIndex2Freq(const int &freq_index) const {
        assert( freq_index >= 0 );
        assert( freq_index < this->FreqNum() );
        return this->freq_min + freq_index * this->freq_interval;
    }

    int FreqGrids::Freq2FreqIndex(const double &freq) const {
        assert( freq >= this->freq_min && freq < this->freq_max );
        int index = std::ceil((freq - freq_min) / freq_interval);
        assert( index >= this->int_freq_min && index < this->int_freq_max );
        return index;
    }

    double FreqGrids::SpecIndex2Freq(const int &spec_index) const {
        assert( spec_index >= 0 );
        assert( spec_index < this->SpecNum() );
        return this->freq_min + spec_index * this->spec_interval;
    }

    int FreqGrids::FreqIndex2SpecIndex(const int &freq_index) const {
        assert( freq_index >= 0 );
        assert( freq_index < this->FreqNum() );
        return std::floor(freq_index * this->freq_interval / this->spec_interval);
    }

} // namespace Grids
