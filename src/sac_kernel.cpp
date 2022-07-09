#include "sac_kernel.h"
#include "sac_core.h"
#include "freq_grids.h"
#include <iostream>

namespace SAC {

    // interfaces
    int Kernel::time_size() const { return this->m_nt; }
    int Kernel::freq_size() const { return this->m_nfreq; }
    const Eigen::MatrixXd& Kernel::kernel() const { return this->m_kernel; }

    double Kernel::kernel( int t, int freq ) const 
    {
        assert( t >= 0 && t < this->m_nt );
        assert( freq >= 0 && freq < this->m_nfreq );
        return this->m_kernel(t, freq);
    }


    Kernel::Kernel( int nt, int nfreq )
    {
        this->m_nt = nt;
        this->m_nfreq = nfreq;
        this->m_kernel.resize(nt, nfreq);
    }


    void Kernel::initial( const SAC::SacCore& sac, 
                          const Grids::FreqGrids& grids, 
                          std::string_view kernel_type ) 
    {
        assert( sac.tau_from_qmc.size() == this->m_nt );
        assert( grids.FreqNum() == this->m_nfreq );

        // kernel for fermionic green's functions
        // e.g. fermionic green's function  ->  usual fermionic spectral function
        if ( kernel_type == "fermion" ) {
            for (auto i = 0; i < grids.FreqNum(); ++i) {
                const double freq = grids.FreqIndex2Freq(i);
                this->m_kernel.col(i) = (-freq*sac.tau_from_qmc.array()).exp() / (1.0 + exp(-sac.beta*freq));
            }
        }

        // kernel for bosonic green's functions
        // e.g. bosonic Matsubara correlation function ->  dynamic susceptibility
        if ( kernel_type == "boson" ) {
            for (auto i = 0; i < grids.FreqNum(); ++i) {
                const double freq = grids.FreqIndex2Freq(i);
                this->m_kernel.col(i) = ( (-freq*sac.tau_from_qmc.array()).exp() + (-freq*(sac.beta-sac.tau_from_qmc.array())).exp() ) / (1.0 + exp(-sac.beta*freq));
            }
        }
    }


    void Kernel::rotate( const Eigen::MatrixXd& rotate_mat ) 
    {
        assert( rotate_mat.rows() == this->m_nt );
        assert( rotate_mat.cols() == this->m_nt );

        this->m_kernel = rotate_mat * this->m_kernel;
        // for (int i = 0; i < this->m_nt; ++i) {
        //     this->m_kernel.col(i) = rotate_mat * this->m_kernel.col(i);
        // }
    }

} // namespace SAC
