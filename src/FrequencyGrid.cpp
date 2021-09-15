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
    this->num_grid_index = this->int_omega_max;
    this->num_spec_index = ceil((omega_max - omega_min) / spec_interval);
}

const double Grid::FrequencyGrid::GridInterval() const{
    return this->grid_interval;
}

const double Grid::FrequencyGrid::SpecInterval() const {
    return this->spec_interval;
}

const int Grid::FrequencyGrid::GridsNum() const {
    return this->num_grid_index;
}

const int Grid::FrequencyGrid::SpecNum() const {
    return this->num_spec_index;
}

const double Grid::FrequencyGrid::GridIndex2Freq(const int &grid_index) const{
    assert( grid_index >= 0 );
    assert( grid_index < this->GridsNum() );
    return this->omega_min + grid_index * this->grid_interval;
}

const int Grid::FrequencyGrid::Freq2GridIndex(const double &freq) const {
    assert( freq >= this->omega_min && freq < this->omega_max );
    int grid = ceil((freq - omega_min) / grid_interval);
    assert( grid >= this->int_omega_min && grid < this->int_omega_max );
    return grid;
}

const double Grid::FrequencyGrid::SpecIndex2Freq(const int &spec_index) const{
    assert( spec_index >= 0 );
    assert( spec_index < this->SpecNum() );
    return this->omega_min + spec_index * this->spec_interval;
}

const int Grid::FrequencyGrid::Grid2Spec(const int &grid_index) const{
    assert( grid_index >= 0 );
    assert( grid_index < this->GridsNum() );
    return floor(grid_index * this->grid_interval / this->spec_interval);
}
