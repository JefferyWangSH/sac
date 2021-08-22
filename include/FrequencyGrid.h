#ifndef SAC_FREQENCYGRID_H
#define SAC_FREQENCYGRID_H
#pragma once

/**
  * This header file includes FrequencyGrid class
  * for the construction of grids in frequency domain.
  */

class FrequencyGrid{

private:
    /* griding params */
    double omega_min{};                 // minimum of frequency
    double omega_max{};                 // maximum of frequency

    int int_omega_min{};                // min and max of frequency number, in unit of griding interval
    int int_omega_max{};

    double grid_interval{};             // delta frequency of sampling space
    double spec_interval{};             // delta frequency of spectrum (frequency spacing in accumulated histogram)


public:

    FrequencyGrid() = default;

    void set_params(double _grid_interval, double _spec_interval, double _omega_min, double _omega_max);

    void init();

    const int lower() const;

    const int upper() const;

    const double interval() const;

    const double index_to_freq(const int &grid_index);

};


#endif //SAC_FREQENCYGRID_H
