#ifndef SAC_WRITER_H
#define SAC_WRITER_H
#pragma once

/**
  *  This header file defines `SAC::Writer` class for the output of SAC results,
  *  including the log of the annealing process, 
  *  recovered spectral functions and the quality report.
  */

#include <string>

// forward declaration
namespace Grids { class FreqGrids; }


namespace SAC {

    // forward declarations
    class SacCore;
    class Kernel;
    class Measure;
    namespace Annealing { class Chain; }
    namespace Initializer { class QmcReader; }
    

    // ----------------------------------------  SAC::Writer class  ------------------------------------------
    class Writer {
        
        public:
            
            // log output during the annealing process
            static void write_log             ( const std::string& file, int n,
                                                const SacCore& core,
                                                const Grids::FreqGrids& grids,
                                                const Measure& measure,
                                                const Annealing::Chain& chain );

            // output the recovered spectral functions
            static void write_spectrum        ( const std::string& file,
                                                const SacCore& core );
            
            // output the quality report of the recovery of spectral functions.
            // Both the correlations from QMC and the recovered ones obtained from the SAC spectrum,
            // are written into the output file, together with their differences.
            static void write_quality_report  ( const std::string& file,
                                                const SacCore& core,
                                                const Grids::FreqGrids& grids,
                                                const Kernel& kernel,
                                                const Initializer::QmcReader& qmc_reader );

    };


} // namespace SAC

#endif // SAC_WRITER_H