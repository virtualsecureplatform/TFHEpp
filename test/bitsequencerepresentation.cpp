#include <cassert>
#include <chrono>
#include <iostream>
#include <random>
#include <tfhe++.hpp>

template<class P, uint32_t int_widht>
using TBSR = std::array<TFHEpp::TRLWE<P>,int_width>

int main()
{
    constexpr uint32_t int_width = 8;

    constexpr uint32_t num_test = 10;
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<uint32_t> binary(0, 1);
}