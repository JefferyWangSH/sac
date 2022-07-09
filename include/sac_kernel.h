#ifndef SAC_KERNEL_H
#define SAC_KERNEL_H
#pragma once

/**
  *  This header file defines `SAC::Kernel` class,
  *  which illustrates the relationships between imaginary-time correlation function G(k,t) 
  *  and the spectral function A(k,omega).
  *  Both fermionic and bosonic transformations are implemented,
  *  and customized kernels for other physical quantities are also supported.
  */

#include <string_view>
#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#include <Eigen/Core>

// forward declaration
namespace Grids { class FreqGrids; }


namespace SAC {
    
    class SacCore;

    // -------------------------------------  SAC::Kernel class  -----------------------------------------
    class Kernel {

        private:
            int m_time_size{};                // dimension of time
            int m_freq_size{};                // dimension of (hyperfine) frequency
            Eigen::MatrixXd m_kernel{};       // kernel matrix

        public:

            Kernel() = default;
            Kernel( int time_size, int freq_size );
            
            // initialize the kernel object from time and frequency grids
            void initial ( const SAC::SacCore& sac, 
                           const Grids::FreqGrids& grids, 
                           std::string_view kernel_type = "fermion" );

            // rotate the kernel to the diagonal representation of the covariance matrix
            void rotate  ( const Eigen::MatrixXd& rotate_mat );

            // interface member functions
            int time_size() const;
            int freq_size() const;
            const Eigen::MatrixXd& kernel()  const;
            // double kernel( int t, int freq ) const;

    };

} // namespace SAC

#endif // SAC_KERNEL_H
