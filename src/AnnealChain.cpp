#include "AnnealChain.h"
#include <iostream>

Annealing::AnnealChain::AnnealChain(int len) {
    this->len = len;
    this->chain.reserve(len);
}

void Annealing::AnnealChain::push(const Annealing::DeltaData &data) {
    if (this->chain.size() < len) {
        this->chain.emplace_back(data);
    }
    else {
        std::cerr << " annealing chain longer than permitted ! " << std::endl;
        exit(1);
    }
}

const int Annealing::AnnealChain::maximum_steps() const {
    return this->len;
}
