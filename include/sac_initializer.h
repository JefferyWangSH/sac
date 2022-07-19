#ifndef SAC_INITIALIZER_H
#define SAC_INITIALIZER_H
#pragma once

/**
  *  This header file defines the interface class `SAC::Initializer`
  *  for the initialization of the relevant modules used in the SAC simulation.
  *  Program options are read from the toml configuration file
  *  using the static member function SAC::Initializer::parse_toml().
  */

#include <string_view>


// forward declaration
namespace Grids { class FreqGrids; }

namespace SAC {

    // forward declarations
    class SacCore;
    class Kernel;
    class Measure;
    class QmcReader;
    namespace Annealing { class Chain; }
    
    // --------------------------------  SAC::Initializer class  -------------------------------
    class Initializer {

        public:

            static void parse_toml( std::string_view toml_config,
                                    QmcReader& reader,
                                    SacCore& core,
                                    Kernel& kernel,
                                    Grids::FreqGrids& grids,
                                    Measure& measure,
                                    Annealing::Chain& chain );

            static void initial();

    };

} // namespace SAC


#endif // SAC_INITIALIZER_H
