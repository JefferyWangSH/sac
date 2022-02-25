#ifndef SAC_ANNEALINGCHAIN_H
#define SAC_ANNEALINGCHAIN_H
#pragma once

/**
  *  This header file includes `AnnealingData` struct and `AnnealingChain` class
  *  for the storage of simulating information and params during the annealing process of SAC.
  */

#include <vector>
#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#include <Eigen/Core>


namespace Annealing{

    struct AnnealingData{
        double theta{};                   // sampling temperature
        int ndelta{};                     // number of delta functions
        int window_width{};               // width of random move window
        double amplitude{};               // amplitude of delta functions
        double chi2{};                    // average fitting goodness chi2 at current sampling temperature
        Eigen::VectorXi locations{};      // locations of delta functions
    };

    class AnnealingChain {
    public:
        int length{};
        int max_length{};

        std::vector<AnnealingData> chain{};

    public:
        AnnealingChain() = default;

        AnnealingChain(int len);

        void push(const AnnealingData &data);

        int len() const;

        void clear();
    };

} // namespace Annealing

#endif //SAC_ANNEALINGCHAIN_H
