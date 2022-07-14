#include "sac_writer.h"
#include "sac_core.h"
#include "sac_measure.h"
#include "sac_annealing.h"
#include "freq_grids.h"

#include <iostream>
#include <fstream>
#include <boost/format.hpp>

namespace SAC {

    void Writer::write_log( const std::string& file, int n,
                            const SacCore& core,
                            const Grids::FreqGrids& grids,
                            const Measure& measure,
                            const Annealing::Chain& chain )
    {   
        // n labels the index of bin data in the Measure class
        // first check if the log file is empty
        std::ifstream check_empty( file, std::ios::in|std::ios::app );
        if ( !check_empty.is_open() ) {
            std::cerr << "SAC::Writer::write_log(): "
                      << "fail to open file \'" << file << "\'." << std::endl;
            exit(1);
        }
        const bool is_empty = ( check_empty.peek() == EOF );
        check_empty.close();

        std::ofstream outfile( file, std::ios::out|std::ios::app );
        if ( !outfile.is_open() ) {
            std::cerr << "SAC::Writer::write_log(): "
                      << "fail to open file \'" << file << "\'." << std::endl;
            exit(1);
        }

        if ( is_empty ) {
            // if empty, print the header infomation
            boost::format header_format("%| 15s|%| 13s|%| 15s|%| 15s|%| 18s|%| 15s|%| 15s|%| 15s|");
            outfile << header_format % "AnnealingStep" % "BinIndex" % "Theta" % "MinChi2/nt" % "AveragedChi2/nt" 
                                     % "DeltaChi2" % "AcceptRatio" % "WindowWidth" << std::endl;
        }

        boost::format log_format("%| 15d|%| 13d|%| 15.3e|%| 15.3e|%| 18.3e|%| 15.3e|%| 15.5f|%| 15.5e|");
        outfile << log_format % ( chain.length() + 1 )
                              % ( n + 1 ) 
                              % ( core.Theta() )
                              % ( core.minChi2() / core.TimeSize() )
                              % ( measure.chi2(n) / core.TimeSize() )
                              % ( measure.chi2(n) - core.minChi2() )
                              % ( measure.accept_ratio(n) )
                              % ( core.WindowWidth() * grids.FreqInterval() ) << std::endl;
        outfile.close();
    }

} // namespace SAC