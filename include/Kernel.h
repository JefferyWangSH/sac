#ifndef SAC_KERNEL_H
#define SAC_KERNEL_H
#pragma once

/**
  * This header file includes Kernel class,
  * which describe relationships between time correlations and spectrum.
  * Different kinds of kernel for various physical system are provided in this file,
  * including kernels for fermion, boson system, and ...
  */

#include <string>

#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#include <Eigen/Core>

namespace Simulation { class SAC; }
namespace  Grid { class FrequencyGrid; }

namespace Kernel {
    class Kernel {
    public:
        int nt;
        int nfreq;

        Eigen::MatrixXd kernel;

    public:
        Kernel() = default;

        Kernel(int nt, int nfreq);

        void init(const Simulation::SAC &sac, Grid::FrequencyGrid grid, const std::string &mode="fermion");

        void rotate(const Eigen::MatrixXd &rotate_mat);
    };
}

#endif //SAC_KERNEL_H
