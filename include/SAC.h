#ifndef SAC_SAC_H
#define SAC_SAC_H
#pragma once

/**
 *  This head file includes SAC class for the implementation of stochastic analytic continuation method,
 *  proposed by Anders.W. Sandvik.
 *  A Monte Carlo process and simulated annealing are performed
 *  to extract real-frequency information from imaginary-time correlations,
 *  which are obtained previously by QMC calculations.
 *
 */


#include <random>
#include "ReadInModule.h"

// random engine
static std::default_random_engine rand_engine_sac(time(nullptr));


class SAC {
public:


};

#endif //SAC_SAC_H
