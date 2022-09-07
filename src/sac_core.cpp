#include "sac_core.h"
#include "sac_kernel.h"
#include "sac_measure.h"
#include "sac_writer.h"
#include "qmc_reader.h"
#include "freq_grids.h"
#include "random.h"

#include <iostream>

namespace SAC {

    // interface member functions
    int SacCore::TimeSize() const { return this->m_time_size; }
    int SacCore::WindowWidth() const { return this->m_metadata.window_width; }
    int SacCore::StabilizationPace() const { return this->m_stabilization_pace; }
    int SacCore::NumDeltas() const { return this->m_delta_num; }
    int SacCore::CollectingSteps() const { return this->m_collecting_steps; }
    double SacCore::Theta() const { return this->m_metadata.theta; }
    double SacCore::minChi2() const { return this->m_chi2_min; }
    double SacCore::ScalingFactor() const { return this->m_scaling_factor; }
    double SacCore::AnnealingRate() const { return this->m_annealing_rate; }

    const Eigen::VectorXd& SacCore::tGridsQmc() const { return this->m_tgrids_qmc; }
    const Eigen::VectorXd& SacCore::CorrQmc() const { return this->m_corr_qmc; }
    const Eigen::VectorXd& SacCore::SigmaQmc() const { return this->m_sigma_qmc; }
    const Eigen::VectorXd& SacCore::FrequencyGrids() const { return this->m_freq; }
    const Eigen::VectorXd& SacCore::RecoveredSpectrum() const { return this->m_spec; }
    
    double SacCore::FrequencyGrids( int i ) const {
        assert( i >= 0 && i < this->m_freq.size() );
        return this->m_freq(i);
    }

    double SacCore::RecoveredSpectrum( int i ) const {
        assert( i >= 0 && i < this->m_spec.size() );
        return this->m_spec(i);
    }


    void SacCore::set_sampling_params( int delta_num, int collecting_steps, int stabilization_pace, 
                                       const std::string& update_type )
    {
        assert( delta_num > 0 );
        assert( collecting_steps > 0 );
        assert( stabilization_pace > 0 );

        this->m_delta_num = delta_num;
        this->m_collecting_steps = collecting_steps;
        this->m_stabilization_pace = stabilization_pace;
        
        if ( update_type != "single" && update_type != "pair" ) {
            std::cerr << "SAC::SacCore::set_sampling_params(): "
                      << "undefined update type \'" << update_type 
                      << "\' ( the default type of the MC updates is \'single\' )." << std::endl;
            exit(1);
        }
        else { this->m_update_type = update_type; }
    }


    void SacCore::set_annealing_params( double theta, double annealing_rate, const std::string& log )
    {   
        assert( theta > 0.0 );
        assert( annealing_rate > 0.0 && annealing_rate < 1.0 );
        this->m_metadata.theta = theta;
        this->m_annealing_rate = annealing_rate;
        this->m_log_file = log;
    }


    void SacCore::initial( const Kernel& kernel, 
                           const Initializer::QmcReader& qmc_reader, 
                           const Grids::FreqGrids& grids )
    {
        // ensure that the QmcReader and FreqGrids classes have been initialized before
        // check the validity of parameters
        if ( kernel.OnlyUsePositiveFreqDomain() && grids.FreqIndex2Freq(0) != 0.0 ) {
            std::cerr << "SAC::SacCore::initial: "
                      << "input frequency domain incompatible with the kernel type." << std::endl; 
            exit(1); 
        }


        // -----------------------------------------------------------------------------------------
        //                               Initialize from QmcReader
        // -----------------------------------------------------------------------------------------

        this->m_time_size = qmc_reader.cov_mat_dim();
        this->m_beta = qmc_reader.beta();
        this->m_scaling_factor = qmc_reader.scaling_factor();

        this->m_tgrids_qmc = qmc_reader.tgrids_qmc();
        this->m_corr_qmc = qmc_reader.rotate_mat() * qmc_reader.corr_mean_qmc();
        this->m_sigma_qmc = ( std::sqrt(qmc_reader.bootstrap_num()) / qmc_reader.eig_vec().array().sqrt() ).matrix();


        // -----------------------------------------------------------------------------------------
        //                               Initialize from FreqGrids
        // -----------------------------------------------------------------------------------------
        
        // ----------------  Generate initial configurations of the delta functions  ---------------

        // reading grids info from FreqGrids
        this->m_metadata.locations.resize(this->m_delta_num);

        // one can choose customized strategies for initialization 
        // according to different priori information of the spectral functions
        // e.g. randomly distributed
        std::uniform_int_distribution<> rand_delta(0, grids.FreqNum()-1);
        for ( auto i = 0; i < this->m_delta_num; ++i ) {
            this->m_metadata.locations(i) = rand_delta( Utils::Random::Engine );
        }

        // // delta-like distribution
        // const double delta_freq = 2.0;
        // this->m_metadata.locations = Eigen::VectorXi::Constant( this->m_delta_num, grids.Freq2FreqIndex(delta_freq) );

        // // rectangle-like distribution
        // const double left_edge = -2.0;
        // const double right_edge = 2.0;
        // const int left_edge_index = grids.Freq2FreqIndex(left_edge);
        // const int right_edge_index = grids.Freq2FreqIndex(right_edge);
        // const int rectangle_length = right_edge_index - left_edge_index + 1;
        // for ( auto i = 0; i < this->m_delta_num; ++i ) {
        //     this->m_metadata.locations(i) = left_edge_index + i % rectangle_length;
        // }

        // // TODO: gaussian-like distribution
        // const double gaussian_peak = 1.0;
        // const double gaussian_sigma = 1.0;
        // // filling from peak to edges

        // ------------------  Initialize amplitudes of the delta functions  -----------------------
        
        // all delta functions share the equal amplitude
        // the factor 1.0 comes from the normalization of spectral functions:
        // for fermionic systems, this is an intrinsic property of spectral functions;
        // while for bosonic systems, this is guaranteed by scaling the correlations such that G(t=0) = 1.0
        if ( kernel.NeedManualNormalize() ) {
            this->m_delta_amplitude = 1.0 / this->m_delta_num; 
        }
        else {
            this->m_delta_amplitude = 1.0 / (this->m_scaling_factor * this->m_delta_num);
        }
        
        // --------------  Initialize the window width of delta function moves  --------------------

        // todo: 1/10 of average frequency ?
        const double average_freq = std::abs( std::log( 1.0 / qmc_reader.corr_mean_qmc()[this->m_time_size-1] ) 
            / this->m_tgrids_qmc[this->m_time_size-1] );
        this->m_metadata.window_width = std::ceil( 0.1 * average_freq / grids.FreqInterval() );
        // this->m_metadata.window_width = std::ceil( 0.5 * grids.FreqNum() );


        // -----------------------------------------------------------------------------------------
        //                           Allocate memory for others members
        // -----------------------------------------------------------------------------------------

        // initialize correlation functions
        this->m_corr_now.resize(this->m_time_size);
        this->m_corr_next.resize(this->m_time_size);
        this->compute_corr_from_deltas( kernel );

        // compute the goodness of fitting (chi2) for the initial configurations of delta functions
        this->m_chi2 = this->compute_goodness(this->m_corr_now);
        this->m_chi2_min = this->m_chi2;

        this->m_accept_ratio = 0.0;

        // initialize containers for recovering spectrum
        this->m_freq.resize( grids.SpecNum() );
        this->m_spec.resize( grids.SpecNum() );
    }



    void SacCore::compute_corr_from_deltas( const Kernel& kernel )
    {
        // require Eigen version > 3.3.9
        const Eigen::MatrixXd& tmp_kernel = kernel.kernel()(Eigen::all, this->m_metadata.locations);
        this->m_corr_now = tmp_kernel * Eigen::VectorXd::Constant(this->m_delta_num, this->m_delta_amplitude);
        // todo: test the efficiency
    }


    double SacCore::compute_goodness( const Eigen::VectorXd& corr ) const
    {
        assert( corr.size() == this->m_time_size );
        return ( ( corr - this->m_corr_qmc ).array() * this->m_sigma_qmc.array() ).square().sum();
    }


    void SacCore::update_deltas_1step( const Kernel& kernel, const Grids::FreqGrids& grids )
    {
        if ( this->m_update_type == "single" ) { this->update_deltas_1step_single( kernel, grids ); }
        else if ( this->m_update_type == "pair" ) { this->update_deltas_1step_pair( kernel, grids ); }
    }


    // move a delta function randomly for one single update attempt
    // we define that a MonteCarlo step includes `delta_num` number of such update attempts.
    void SacCore::update_deltas_1step_single( const Kernel& kernel, const Grids::FreqGrids& grids )
    {
        // intermediate parameters
        std::uniform_int_distribution<> rand_delta(0, this->m_delta_num-1);
        std::uniform_int_distribution<> rand_width(1, this->m_metadata.window_width);
        std::uniform_int_distribution<> rand_location(0, grids.FreqNum()-1);
        int select_delta;
        int move_width;
        int location_now;
        int location_next;
        double chi2_next;
        double p;
        int accept_count = 0;

        // attempt to move the delta functions for `delta_num` times
        for ( auto i = 0; i < this->m_delta_num; ++i ) {
            // select one delta function randomly
            select_delta = rand_delta( Utils::Random::Engine );
            location_now = this->m_metadata.locations[select_delta];

            if ( this->m_metadata.window_width > 0 && this->m_metadata.window_width < grids.FreqNum() ) {
                // move the selected delta function randomly within the window
                move_width = rand_width( Utils::Random::Engine );
                if ( std::bernoulli_distribution(0.5)( Utils::Random::Engine ) ) {
                    // move to right
                    location_next = location_now + move_width;
                }
                else {
                    // move to left
                    location_next = location_now - move_width;
                }

                // abort the update if the updated delta function runs out of the frequency domain
                if ( location_next < 0 || location_next >= grids.FreqNum() ) {
                    i -= 1;
                    continue;
                }
            }
            else if ( this->m_metadata.window_width == grids.FreqNum() ) {
                // randomly move over the frequency domain
                location_next = rand_location( Utils::Random::Engine );
            }
            else {
                std::cerr << "SAC::SacCore::update_deltas_1step_single(): "
                          << "window width of delta function moves larger than the frequency domain." << std::endl; 
                exit(1); 
            }

            // compute the updated correlation functions
            this->m_corr_next = this->m_corr_now + this->m_delta_amplitude *
                    ( kernel.kernel().col(location_next) - kernel.kernel().col(location_now) );

            // compute the chi2 and accepting ratio for the updated correlation functions
            chi2_next = this->compute_goodness( this->m_corr_next );
            p = exp( ( this->m_chi2 - chi2_next ) / ( 2.0 * this->m_metadata.theta ) );

            if ( std::bernoulli_distribution(std::min(p,1.0))( Utils::Random::Engine ) ) {
                // accept the new configurations
                this->m_metadata.locations[select_delta] = location_next;
                this->m_corr_now = this->m_corr_next;
                this->m_chi2 = chi2_next;
                if ( this->m_chi2 < this->m_chi2_min ) { this->m_chi2_min = this->m_chi2; }
                accept_count++;
            }
        }

        // compute the averaged accepting ratio
        this->m_accept_ratio = (double)accept_count / this->m_delta_num;
    }


    // move a pair of delta functions randomly for one single update attempt
    // in this case, a MonteCarlo step includes `delta_num/2` number of such update attempts.
    void SacCore::update_deltas_1step_pair( const Kernel& kernel, const Grids::FreqGrids& grids )
    {
        assert( this->m_delta_num >= 2 );

        // intermediate parameters
        std::uniform_int_distribution<> rand_delta(0, this->m_delta_num-1);
        std::uniform_int_distribution<> rand_width(1, this->m_metadata.window_width);
        std::uniform_int_distribution<> rand_location(0, grids.FreqNum()-1);
        int select_delta1, select_delta2;
        int move_width1, move_width2;
        int location_now1, location_now2;
        int location_next1, location_next2;
        bool out_of_domain1, out_of_domain2;
        double chi2_next;
        double p;
        int accept_count = 0;

        // attempt to move pair of delta functions for `delta_num/2` times
        for ( auto i = 0; i < std::ceil(this->m_delta_num/2); i++ ) {
            // select two different delta functions randomly
            select_delta1 = rand_delta( Utils::Random::Engine );
            select_delta2 = select_delta1;
            while ( select_delta1 == select_delta2 ) {
                select_delta2 = rand_delta( Utils::Random::Engine );
            }
            location_now1 = this->m_metadata.locations[select_delta1];
            location_now2 = this->m_metadata.locations[select_delta2];

            if ( this->m_metadata.window_width > 0 && this->m_metadata.window_width < grids.FreqNum() ) {
                // randomly choose a window width for the delta function moves
                move_width1 = rand_width( Utils::Random::Engine );
                move_width2 = rand_width( Utils::Random::Engine );
                if ( std::bernoulli_distribution(0.5)( Utils::Random::Engine ) ) {
                    location_next1 = location_now1 + move_width1;
                    location_next2 = location_now2 - move_width2;
                }
                else {
                    location_next1 = location_now1 - move_width1;
                    location_next2 = location_now2 + move_width2;
                }

                // abort the update if any of the updated delta functions runs out of the frequency domain
                out_of_domain1 = ( location_next1 < 0 || location_next1 >= grids.FreqNum() );
                out_of_domain2 = ( location_next2 < 0 || location_next2 >= grids.FreqNum() );
                if ( out_of_domain1 || out_of_domain2 ) {
                    i -= 1;
                    continue;
                }
            }
            else if (this->m_metadata.window_width == grids.FreqNum()) {
                // randomly move over the frequency domain
                location_next1 = rand_location( Utils::Random::Engine );
                location_next2 = rand_location( Utils::Random::Engine );
            }
            else {
                std::cerr << "SAC::SacCore::update_deltas_1step_pair(): "
                          << "window width of delta function moves larger than the frequency domain." << std::endl; 
                exit(1); 
            }

            // compute the updated correlation functions
            this->m_corr_next = this->m_corr_now + this->m_delta_amplitude *
                    ( kernel.kernel().col(location_next1) - kernel.kernel().col(location_now1)
                    + kernel.kernel().col(location_next2) - kernel.kernel().col(location_now2) );

            // compute the chi2 and accepting ratio for the updated correlation functions
            chi2_next = this->compute_goodness( this->m_corr_next );
            p = exp( ( this->m_chi2 - chi2_next ) / ( 2.0 * this->m_metadata.theta ) );

            if ( std::bernoulli_distribution(std::min(p,1.0))( Utils::Random::Engine ) ) {
                // accept the new configurations
                this->m_metadata.locations[select_delta1] = location_next1;
                this->m_metadata.locations[select_delta2] = location_next2;
                this->m_corr_now = this->m_corr_next;
                this->m_chi2 = chi2_next;
                if ( this->m_chi2 < this->m_chi2_min ) { this->m_chi2_min = this->m_chi2; }
                accept_count++;
            }
        }

        // compute the averaged accepting ratio
        this->m_accept_ratio = (double)accept_count / std::ceil(this->m_delta_num/2);
    }


    void SacCore::update_at_fixed_theta( const Kernel& kernel, 
                                         const Grids::FreqGrids& grids , 
                                         Measure& measure,
                                         const Annealing::Chain& chain )
    {
        // total Monte Carlo steps at a fixed sampling temperature: nbin * sbin
        for ( auto n = 0; n < measure.number_of_bin(); ++n ) {
            // n corresponds to the index of bins
            for ( auto s = 0; s < measure.size_of_bin(); ++s ) {
                // s corresponds to the index of samples in one bin
                // recalculate the goodness of fitting `chi2` every `stabilization_pace` steps
                if ( s % this->m_stabilization_pace == 1 ) {
                    this->m_chi2 = this->compute_goodness( this->m_corr_now );
                }
                this->update_deltas_1step( kernel, grids );
                measure.collect( s, this->m_chi2, this->m_accept_ratio );
            }

            // average over samples in one bin, labeled by n
            measure.bin_analyse(n);

            // write log
            Writer::write_log( this->m_log_file, n, *this, grids, measure, chain );

            // adjust the window width of the delta function moves
            // make sure that the averaged accepting ratio of random moves is around 0.5
            if ( measure.accept_ratio(n) > 0.5 ) {
                this->m_metadata.window_width = std::min( (int)std::ceil( this->m_metadata.window_width * 1.5 ),
                                                          grids.FreqNum() );
            }
            if ( measure.accept_ratio(n) < 0.4 ) {
                this->m_metadata.window_width = std::ceil( this->m_metadata.window_width / 1.5 );
            }
        }
    }

    
    void SacCore::perform_annealing( const Kernel& kernel, 
                                     const Grids::FreqGrids& grids , 
                                     Measure& measure,
                                     Annealing::Chain& chain )
    {
        // annealing process, with maximum number of steps being `max_annealing_steps`
        for ( auto i = 0; i < chain.max_length(); ++i ) {

            // Monte Carlo updates at the fixed temperature
            this->update_at_fixed_theta( kernel, grids, measure, chain );

            // record simulating information for current sampling temperature
            measure.analyse();
            this->m_metadata.chi2 = measure.chi2();
            chain.push(this->m_metadata);

            // exit condition
            if ( measure.chi2() - this->m_chi2_min < 1e-3 ) { break; }

            // lower down the sampling temperature
            this->m_metadata.theta *= this->m_annealing_rate;
        }
    }

    
    void SacCore::decide_sampling_theta( const Kernel& kernel, const Annealing::Chain& chain )
    {
        // determine the sampling temperature after the annealing process
        // by slightly raising the artificial temperature `theta` to avoid overfitting
        for ( auto i = chain.length() - 1; i >= 0; --i ) {
            // raise chi2 by a standard deviation with respect to its minimum
            if ( chain.chain(i).chi2 > this->m_chi2_min + 2.0 * std::sqrt(this->m_chi2_min) ) {
                this->m_metadata = chain.chain(i);
                break;
            }
        }
        this->compute_corr_from_deltas( kernel );
        this->m_chi2 = this->compute_goodness( this->m_corr_now );
    }


    void SacCore::sample_and_collect( const Kernel& kernel, 
                                      const Grids::FreqGrids& grids , 
                                      Measure& measure,
                                      const Annealing::Chain& chain )
    {
        // reach equilibrium at the current sampling temperature
        this->update_at_fixed_theta( kernel, grids, measure, chain );

        // generate frequency grids for the spectral functions
        for ( auto n = 0; n < grids.SpecNum(); ++n ) {
            this->m_freq(n) = grids.SpecIndex2Freq(n);
        }
        this->m_spec.setZero();

        // sampling and collecting the spectrum
        for ( auto i = 0; i < this->m_collecting_steps; ++i ) {
            // recompute the goodness of fitting `chi2` every `stabilization_pace` steps
            if ( i % this->m_stabilization_pace == 1 ) {
                this->m_chi2 = this->compute_goodness( this->m_corr_now );
            }
            // update configurations of the delta functions
            this->update_deltas_1step( kernel, grids );

            // collect the spectral functions
            for ( auto j = 0; j < this->m_delta_num; ++j ) {
                this->m_spec( grids.FreqIndex2SpecIndex(this->m_metadata.locations(j)) ) += this->m_delta_amplitude;
            }
        }

        // re-scale and recover the spectral functions
        if ( kernel.NeedManualNormalize() ) {
            this->m_spec = M_PI * this->m_scaling_factor / ( this->m_collecting_steps * grids.SpecInterval() ) *
                           this->m_spec.array() / ( 1.0 + ( -this->m_beta * this->m_freq.array() ).exp() ) ;
        }
        else {
            this->m_spec *= this->m_scaling_factor / ( this->m_collecting_steps * grids.SpecInterval() );
        }
    }
    

} // namespace SAC
