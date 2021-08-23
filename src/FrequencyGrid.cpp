#include "FrequencyGrid.h"
#include <cassert>
#include <cmath>


Grid::FrequencyGrid::FrequencyGrid(double _grid_interval, double _spec_interval, double _omega_min, double _omega_max) {
    this->grid_interval = _grid_interval;
    this->spec_interval = _spec_interval;
    this->omega_min = _omega_min;
    this->omega_max = _omega_max;
}

void Grid::FrequencyGrid::init() {
    // convert frequency into unit of griding interval
    this->int_omega_min = 0;
    this->int_omega_max = ceil((omega_max - omega_min) / grid_interval);
}

const int Grid::FrequencyGrid::lower() const{
    return this->int_omega_min;
}

const int Grid::FrequencyGrid::upper() const{
    return this->int_omega_max;
}

const double Grid::FrequencyGrid::interval() const{
    return this->grid_interval;
}

const double Grid::FrequencyGrid::index_to_freq(const int &grid_index) {
    assert( grid_index >= 0 );
    assert( grid_index < int_omega_max );
    return omega_min + grid_index * grid_interval;
}

