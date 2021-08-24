#ifndef SAC_ANNEALCHAIN_H
#define SAC_ANNEALCHAIN_H
#pragma once

/**
  *  This header file includes `DeltaData` struct and `AnnealChain` class
  *  for the storage of simulating information during the annealing process of SAC.
  */


#include <vector>

#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#include <Eigen/Core>


namespace Annealing{

    struct DeltaData{
        double theta{};                 // sampling temperature
        int ndelta{};                   // number of delta functions
        int window_width{};             // width of random move window
        double amplitude{};             // amplitude of delta functions
        Eigen::VectorXi locations;      // locations of delta functions
    };

    class AnnealChain {
    public:
        int length{};
        int max_length{};

        std::vector<DeltaData> chain;

    public:

        AnnealChain() = default;

        AnnealChain(int len);

        void push(const DeltaData &data);

        const int len() const;
    };
}

#endif //SAC_ANNEALCHAIN_H
