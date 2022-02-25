#include "kernel.h"
#include "sac.h"
#include "freq_grids.h"
#include <iostream>

namespace Kernel {

    Kernel::Kernel(int nt, int nfreq) {
        this->nt = nt;
        this->nfreq = nfreq;
        this->kernel.resize(nt, nfreq);
    }

    void Kernel::init(const Simulation::SAC &sac, const Grids::FreqGrids &grids, const std::string &kernel_type) {
        assert( sac.tau_from_qmc.size() == nt );
        assert( grids.FreqNum() == nfreq );

        // kernel for fermion greens function
        if ( kernel_type == "fermion" ) {
            for (int i = 0; i < grids.FreqNum(); ++i) {
                const double freq = grids.FreqIndex2Freq(i);
                kernel.col(i) = (-freq*sac.tau_from_qmc.array()).exp() / (1.0 + exp(-sac.beta*freq));
            }
        }

        // kernel for boson greens function
        if ( kernel_type == "boson" ) {
            for (int i = 0; i < grids.FreqNum(); ++i) {
                const double freq = grids.FreqIndex2Freq(i);
                kernel.col(i) = ( (-freq*sac.tau_from_qmc.array()).exp() + (-freq*(sac.beta-sac.tau_from_qmc.array())).exp() ) / (1.0 + exp(-sac.beta*freq));
            }
        }
    }

    void Kernel::rotate(const Eigen::MatrixXd &rotate_mat) {
        assert( rotate_mat.rows() == nt );
        assert( rotate_mat.cols() == nt );

        kernel = rotate_mat * kernel;
    //    for (int i = 0; i < nt; ++i) {
    //        kernel.col(i) = rotate_mat * kernel.col(i);
    //    }
    }

} // namespace Kernel
