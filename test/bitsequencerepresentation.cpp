#include <cassert>
#include <chrono>
#include <iostream>
#include <random>
#include <tfhe++.hpp>

template<class P, uint32_t int_width>
using TBSR = std::array<TFHEpp::TRLWE<P>,int_width>

template<class P, uint32_t int_width>
inline void TBSRgen(TBSR<P,int_width> &tbsr){
    tbsr = {};
    for(int i = 0; i<int_width;i++) for(int j = 0;j<(1<<int_width);j++) tbsr[i][1][j] = ((j>>i)&1)?P::μ:-P::μ;
    return tbsr;
}

template<class P, uint32_t int_width>
void AddByTBSR(TBSR<P,int_width+1> &res,std::array<TFHEpp::TRGSWFFT<P>,int_width> &A,std::array<TFHEpp::TRGSWFFT<P>,int_width> &B){
    TBSRgen<P,int_width+1>(res);
}

int main()
{
    constexpr uint32_t int_width = 8;
    static_assert(TFHEpp::lvl1param::nbit>=int_width);

    constexpr uint32_t num_test = 10;
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<uint32_t> binary(0, 1);
}