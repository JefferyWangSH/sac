#ifndef SAC_WRITER_H
#define SAC_WRITER_H
#pragma once

/**
  *  This header file defines `SAC::Writer` class for the output of SAC results,
  *  including the log of the annealing process, 
  *  recovered spectral functions and the quality report.
  */


namespace SAC {

    // -------------------------------------  SAC::Writer class  ------------------------------------
    class Writer {
        
        public:
            
            static void write_log();

            static void write_spectrum();

            static void write_quality_report();

    };
}

#endif // SAC_WRITER_H