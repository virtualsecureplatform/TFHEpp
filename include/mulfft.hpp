#include<array>

#include<spqlios-fft.h>

#include<params.hpp>

namespace TFHEpp{
    using namespace std;
    template <uint32_t N = DEF_N>
    array<uint32_t, N> TwistFFT(const array<double,N> &a, FFT_Processor_Spqlios &fftp);
    template <uint32_t N = DEF_N>
    array<double, N> TwistIFFT(const array<uint32_t,N> &a, FFT_Processor_Spqlios &fftp);
    array<uint32_t,DEF_N> PolyMullvl1(const array<uint32_t,DEF_N> &a, const array<uint32_t,DEF_N> &b);
}