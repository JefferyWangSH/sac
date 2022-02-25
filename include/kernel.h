#ifndef SAC_KERNEL_H
#define SAC_KERNEL_H
#pragma once

/**
  *  This header file includes Kernel class,
  *  which describe relationships between correlation function G(k,t) and spectral function A(k,omega).
  *  Different kinds of kernel for various physical system are supported,
  *  including kernels for fermion, boson system, and ...
  */

#include <string>
#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#include <Eigen/Core>

namespace Simulation { class SAC; }
namespace Grids { class FreqGrids; }

namespace Kernel {

    class Kernel {
    public:
        int nt{};                   // dimension of time
        int nfreq{};                // dimension of frequency
        Eigen::MatrixXd kernel{};   // kernel matrix

    public:
        Kernel() = default;

        Kernel(int nt, int nfreq);

        /* init kernel from grids in frequency and time */
        void init(const Simulation::SAC &sac, const Grids::FreqGrids &grids, const std::string &kernel_type="fermion");

        /* transformation of picture (bases) */
        void rotate(const Eigen::MatrixXd &rotate_mat);
    };

} // namespace Kernel

#endif //SAC_KERNEL_H
