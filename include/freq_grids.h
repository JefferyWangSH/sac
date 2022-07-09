#ifndef SAC_FREQ_GRIDS_H
#define SAC_FREQ_GRIDS_H
#pragma once

/**
  *  This header file defines `Grids::FreqGrids` class
  *  for the construction of (hyperfine) grids in frequency domain.
  */

namespace Grids {

    // ---------------------------------------  Grids::FreqGrids class  -----------------------------------------
    class FreqGrids{

        private:
            // griding params
            double m_freq_min{};          // minimum of the frequency
            double m_freq_max{};          // maximum of the frequency

            int m_int_freq_min{};         // minimum of the frequency indices, in unit of the griding interval
            int m_int_freq_max{};         // maximum of the frequency indices, in unit of the griding interval

            double m_freq_interval{};     // (hyperfine) frequency interval in the sampling space
            double m_spec_interval{};     // frequency interval of spectrum (frequency spacing in the accumulated histogram)

            int m_num_freq{};             // number of (hyperfine) discrete frequency points
            int m_num_spec{};             // number of discrete frequency points in the accumulated spectrum

        public:

            FreqGrids() = default;
            FreqGrids( double freq_interval, double spec_interval, double freq_min, double freq_max );

            void initial();

            // frequency interval of the hyperfine grids
            double FreqInterval() const;

            // frequency interval of the accumulated spectrum
            double SpecInterval() const;

            // number of discrete hyperfine frequency points
            int FreqNum() const;

            // number of discrete frequency points in the accumulated spectrum
            int SpecNum() const;

            // conversion between indices of discrete (hyperfine) frequency points and continuum frequency
            int Freq2FreqIndex( double freq )         const;
            double FreqIndex2Freq( int freq_index )   const;

            // convert the indices of discrete (spectral) frequency points to continuum frequency
            double SpecIndex2Freq( int spec_index )   const;

            // convert the index of hyperfine frequency grids to that of the spectral ones
            int FreqIndex2SpecIndex( int freq_index ) const;

    };

} // namespace Grids

#endif // SAC_FREQ_GRIDS_H
