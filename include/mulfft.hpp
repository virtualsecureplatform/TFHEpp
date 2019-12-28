#include<array>

#include<ipp.h>

#include<params.hpp>

namespace TFHEpp{
    template <uint32_t N = DEF_N>
    std::array<double, N> TwistFFT(std::array<double,N> a, std::array<Ipp64f,N> twist);
    template <uint32_t N = DEF_N>
    std::array<double, N> TwistIFFT(std::array<double,N> a, std::array<Ipp64f,N> twist);
    std::array<uint32_t,DEF_N> PolyMullvl1(std::array<uint32_t,DEF_N> a, std::array<uint32_t,DEF_N> b);
}