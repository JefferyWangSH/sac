#include <iostream>
#include "Kernel.h"
#include "SAC.h"
#include "FrequencyGrid.h"

void Kernel::set_params(int _nt, int _nfreq) {
    this->nt = _nt;
    this->nfreq = _nfreq;
    this->kernel.resize(nt, nfreq);
}

void Kernel::init(const SAC &sac, FrequencyGrid grid, const std::string &mode) {
    assert( sac.tau.size() == nt );
    assert( grid.upper() == nfreq );
    assert( mode == "fermion" );

    if ( mode == "fermion" ) {
        for (int i = grid.lower(); i < grid.upper(); ++i) {
            const double freq = grid.index_to_freq(i);
            kernel.col(i) = (- freq * sac.tau.array()).exp().matrix() / ( 1.0 + exp(- sac.beta * freq) );
        }
    }

    if ( mode == "boson" ) {

    }
}

void Kernel::rotate(const Eigen::MatrixXd &rotate_mat) {
    assert( rotate_mat.rows() == nt );
    assert( rotate_mat.cols() == nt );

    for (int i = 0; i < nt; ++i) {
        kernel.col(i) = rotate_mat * kernel.col(i);
    }
}

