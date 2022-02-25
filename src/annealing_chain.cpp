#include "annealing_chain.h"
#include <iostream>

namespace Annealing {

    AnnealingChain::AnnealingChain(int len) {
        this->length = 0;
        this->max_length = len;
        this->chain.reserve(len);
    }

    void AnnealingChain::push(const AnnealingData &data) {
        if ( this->len() < this->max_length ) {
            this->chain.emplace_back(data);
            this->length++;
        }
        else {
            std::cerr << " Annealing chain longer than permitted ! \n" << std::endl;
            exit(1);
        }
    }

    int AnnealingChain::len() const {
        return this->length;
    }

    void AnnealingChain::clear() {
    //    this->chain.clear();
    //    this->chain.shrink_to_fit();
        this->length = 0;
    }

} // namespace Annealing
