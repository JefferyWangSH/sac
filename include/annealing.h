#ifndef SIMULATED_ANNEALING_H
#define SIMULATED_ANNEALING_H
#pragma once

/**
  *  This header file defines `SimulatedAnnealing::MetaData` structrue 
  *  and `SimulatedAnnealing::Chain` class 
  *  for the storage of simulating information and params 
  *  during the simulated annealing process of SAC.
  */

#include <vector>
#define EIGEN_USE_MKL_ALL
#define EIGEN_VECTORIZE_SSE4_2
#include <Eigen/Core>


namespace SimulatedAnnealing{
    
    // -----------------------------  SimulatedAnnealing::MetaData structure  --------------------------------
    struct MetaData{

        double theta{};                         // sampling temperature
        int ndelta{};                           // number of delta functions
        int window_width{};                     // width of random move window
        double amplitude{};                     // amplitude of delta functions
        double chi2{};                          // averaged fitting goodness chi2 at current sampling temperature
        Eigen::VectorXi locations{};            // locations of delta functions
    
    };


    // ---------------------------------  SimulatedAnnealing::Chain class  -----------------------------------
    class Chain {
        
        private:
            int m_length{};                     // length of the metadata chain
            int m_max_length{};                 // maximum length of the metadata chain

            std::vector<MetaData> m_chain{};    // chain object

        public:

            Chain() = default;
            Chain( int max_length );

            // push back the metadata
            void push( const MetaData& data ); 

            // clear up metadata in storage
            void clear();

            // interface member functions
            int length()     const;
            int max_length() const;
            const MetaData& chain( int index ) const;
 
    };

} // namespace SimulatedAnnealing

#endif // SIMULATED_ANNEALING_H
