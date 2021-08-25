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

        int num_grid_index{};               // number of (fine) discrete frequencies in grids
        int num_spec_index{};               // number of discrete frequencies in accumulated spectrum


    public:

        FrequencyGrid() = default;

        FrequencyGrid(double grid_interval, double spec_interval, double omega_min, double omega_max);

        void init();

        /* frequency interval in fine grids */
        const double GridInterval() const;

        /* frequency interval in accumulated spectrum */
        const double SpecInterval() const;

        /* number of intervals in fine grids */
        const int GridsNum() const;

        /* number of intervals in accumulated spectrum */
        const int SpecNum() const;

        /* convert discrete index of grids to continuum frequency */
        const double GridIndex2Freq(const int &grid_index) const;

        /* convert discrete index of spectrum to continuum frequency */
        const double SpecIndex2Freq(const int &spec_index) const;

        /* convert index of grids to index of spectrum bins */
        const int Grid2Spec(const int &grid_index) const;

    };
}


#endif //SAC_FREQENCYGRID_H
