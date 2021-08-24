#ifndef SAC_FREQENCYGRID_H
#define SAC_FREQENCYGRID_H
#pragma once

/**
  * This header file includes FrequencyGrid class
  * for the construction of grids in frequency domain.
  */

namespace Grid {
    class FrequencyGrid{
    private:
        /* griding params */
        double omega_min{};                 // minimum of frequency ( in space of continuum )
        double omega_max{};                 // maximum of frequency ( in space of continuum )

        int int_omega_min{};                // min and max of frequency index, in unit of griding interval
        int int_omega_max{};

        double grid_interval{};             // delta frequency of sampling space
        double spec_interval{};             // delta frequency of spectrum (frequency spacing in accumulated histogram)


    public:

        FrequencyGrid() = default;

        FrequencyGrid(double grid_interval, double spec_interval, double omega_min, double omega_max);

        void init();

        /* lower index of discrete frequency */
        const int lower() const;

        /* upper index of discrete frequency */
        const int upper() const;

        /* frequency interval in fine grids */
        const double interval() const;

        /* convert discrete index to continuum frequency */
        const double index_to_freq(const int &grid_index);

    };
}


#endif //SAC_FREQENCYGRID_H
