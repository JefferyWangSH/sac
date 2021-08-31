#include "AnnealChain.h"
#include <iostream>

Annealing::AnnealChain::AnnealChain(int len) {
    this->length = 0;
    this->max_length = len;
    this->chain.reserve(len);
}

void Annealing::AnnealChain::push(const Annealing::AnnealData &data) {
    if ( this->len() < this->max_length ) {
        this->chain.emplace_back(data);
        this->length++;
    }
    else {
        std::cerr << "annealing chain longer than permitted !" << std::endl;
        exit(1);
    }
}

const int Annealing::AnnealChain::len() const {
    return this->length;
}

void Annealing::AnnealChain::clear() {
//    this->chain.clear();
//    this->chain.shrink_to_fit();
    this->length = 0;
}
