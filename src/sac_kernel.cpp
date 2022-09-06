#include "sac_kernel.h"
#include "qmc_reader.h"
#include "freq_grids.h"
#include <iostream>

namespace SAC {

    // interface member functions
    int Kernel::time_size() const { return this->m_time_size; }
    int Kernel::freq_size() const { return this->m_freq_size; }
    const Eigen::MatrixXd& Kernel::kernel() const { return this->m_kernel; }
    
    // double Kernel::kernel( int t, int freq ) const 
    // {
    //     assert( t >= 0 && t < this->m_time_size );
    //     assert( freq >= 0 && freq < this->m_freq_size );
    //     return this->m_kernel(t, freq);
    // }

    bool Kernel::NeedManualNormalize() const { return ( this->m_kernel_type == "boson" ); }
    bool Kernel::OnlyUsePositiveFreqDomain() const { return ( this->m_kernel_type == "boson" ); }


    void Kernel::set_kernel_params( int time_size, int freq_size, const std::string& kernel_type )
    {   
        if ( kernel_type != "fermion" && kernel_type != "boson" ) {
            std::cerr << "SAC::Kernel::set_kernel_params(): "
                      << "undefined kernel type \'" << kernel_type 
                      << "\' ( the default type of the kernel is \'fermion\' )." << std::endl;
            exit(1);
        }
        else { this->m_kernel_type = kernel_type; }
        
        this->m_time_size = time_size;
        this->m_freq_size = freq_size;
        this->m_kernel.resize(time_size, freq_size);
    }


    void Kernel::initial( const Initializer::QmcReader& qmc_reader, const Grids::FreqGrids& grids ) 
    {
        assert( qmc_reader.tgrids_qmc().size() == this->m_time_size );
        assert( grids.FreqNum() == this->m_freq_size );

        // kernel for fermionic green's functions
        // e.g. fermionic green's function  ->  usual fermionic spectral function
        //
        //      G(t) = int_omega{-inf}{+inf}  exp( -omega * t ) / ( 1 + exp( -omega * beta ) ) * A( omega )
        //
        // where A( omega ) is positive definite and satisfies the sum rule sum_omega A = 1

        if ( this->m_kernel_type == "fermion" ) {
            for ( auto i = 0; i < grids.FreqNum(); ++i ) {
                const double freq = grids.FreqIndex2Freq(i);
                this->m_kernel.col(i) = ( -freq * qmc_reader.tgrids_qmc().array() ).exp() 
                                      / ( 1.0 + exp( -qmc_reader.beta() * freq ) );
            }
        }

        // kernel for bosonic correlation functions
        // e.g. dynamic spin structure factor ( susceptibility )
        // 
        //      G(t) = int_omega{0}{+inf}  exp( -omega * t ) * B( omega )
        // 
        // where B( omega ) is positive definite in (0, +inf) and B( -omega ) = exp( -beta*omega ) B( omega ) 
        // the sum rule is guaranteed by rescaling G(t=0) = 1

        if ( this->m_kernel_type == "boson" ) {
            for ( auto i = 0; i < grids.FreqNum(); ++i ) {
                const double freq = grids.FreqIndex2Freq(i);
                this->m_kernel.col(i) = ( -freq * qmc_reader.tgrids_qmc().array() ).exp();
            }
        }

        // rotate to the diagonal representation of the QMC correlation functions
        this->rotate( qmc_reader.rotate_mat() );
    }


    void Kernel::rotate( const Eigen::MatrixXd& rotate_mat ) 
    {
        assert( rotate_mat.rows() == this->m_time_size );
        assert( rotate_mat.cols() == this->m_time_size );

        this->m_kernel = rotate_mat * this->m_kernel;
        // for ( int i = 0; i < this->m_time_size; ++i ) {
        //     this->m_kernel.col(i) = rotate_mat * this->m_kernel.col(i);
        // }
    }

} // namespace SAC
