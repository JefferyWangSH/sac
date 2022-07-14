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
    class Measure;
    namespace Annealing { class Chain; }
    

    // ----------------------------------------  SAC::Writer class  ------------------------------------------
    class Writer {
        
        public:
            
            static void write_log( const std::string& file, int n,
                                   const SacCore& core,
                                   const Grids::FreqGrids& grids,
                                   const Measure& measure,
                                   const Annealing::Chain& chain );

            static void write_spectrum();

            static void write_quality_report();

    };


}

#endif // SAC_WRITER_H