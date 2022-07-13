#include "sac_annealing.h"
#include <iostream>

namespace SAC::Annealing {

    // interface member functions
    int Chain::length() const { return this->m_length; }
    int Chain::max_length() const { return this->m_max_length; }

    const MetaData& Chain::chain( int index ) const 
    {
        assert( index >= 0 && index < this->m_length );
        return this->m_chain[index];
    }


    Chain::Chain( int max_length )
    {
        this->m_length = 0;
        this->m_max_length = max_length;
        this->m_chain.reserve(max_length);
    }


    void Chain::set_max_length( int max_length )
    {
        this->m_max_length = max_length;
        this->m_length = std::min(this->m_length, this->m_max_length);
        this->m_chain.reserve(this->m_max_length);
    }


    void Chain::push( const MetaData& data )
    {
        if ( this->m_length < this->m_max_length ) {
            this->m_chain.emplace_back(data);
            this->m_length++;
        }
        else {
            std::cerr << "SAC::Annealing::Chain::push(): "
                      << "annealing chain longer than permitted." << std::endl;
            exit(1);
        }
    }


    void Chain::clear()
    {
        // clear the annealing history ( chain data )
        this->m_length = 0;
        this->m_chain.clear();
        this->m_chain.shrink_to_fit();
    }

} // namespace SAC::Annealing
