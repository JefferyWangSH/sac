#include "sac_initializer.h"
#include "sac_core.h"
#include "sac_kernel.h"
#include "sac_annealing.h"
#include "sac_measure.h"
#include "freq_grids.h"
#include "qmc_reader.h"
#include "utils/toml.hpp"

namespace SAC {

    void Initializer::parse_toml( std::string_view toml_config,
                                  QmcReader& reader,
                                  SacCore& core,
                                  Kernel& kernel,
                                  Grids::FreqGrids& grids,
                                  Measure& measure,
                                  Annealing::Chain& chain )
    {
        // parse the configuration file
        auto config = toml::parse_file( toml_config );



    }

} // namespace SAC
