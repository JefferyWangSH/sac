#ifndef SAC_FREQ_GRIDS_H
#define SAC_FREQ_GRIDS_H
#pragma once

/**
  *  This header file includes FreqGrids class
  *  for the construction of grids in frequency domain.
  */

namespace Grids {

    class FreqGrids{
    private:
        /* griding params */
        double freq_min{};         // minimum of frequency ( in space of continuum )
        double freq_max{};         // maximum of frequency ( in space of continuum )

        int int_freq_min{};        // min and max of frequency index, in unit of griding interval
        int int_freq_max{};

        double freq_interval{};     // frequency interval in sampling space
        double spec_interval{};     // frequency interval of spectrum (frequency spacing in accumulated histogram)

        int num_freq{};             // number of (fine) discrete frequencies in grids
        int num_spec{};             // number of discrete frequencies in accumulated spectrum

    public:
        FreqGrids() = default;

        FreqGrids(double freq_interval, double spec_interval, double freq_min, double freq_max);

        void init();

        /* frequency interval in fine grids */
        double FreqInterval() const;

        /* frequency interval in accumulated spectrum */
        double SpecInterval() const;

        /* number of intervals in fine grids */
        int FreqNum() const;

        /* number of intervals in accumulated spectrum */
        int SpecNum() const;

        /* conversion between discrete index of grids to continuum frequency */
        double FreqIndex2Freq(const int &freq_index) const;
        int Freq2FreqIndex(const double &freq) const;

        /* convert discrete index of spectrum to continuum frequency */
        double SpecIndex2Freq(const int &spec_index) const;

        /* convert index of frequency grids to index of spectrum bins */
        int FreqIndex2SpecIndex(const int &freq_index) const;
    };

} // namespace Grids

#endif //SAC_FREQ_GRIDS_H
