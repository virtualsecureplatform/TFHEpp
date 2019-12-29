#include<array>

#include<spqlios-fft.h>

#include<params.hpp>

namespace TFHEpp{
    using namespace std;
    void TwistFFTlvl1(array<uint32_t, DEF_N> &res,const array<double,DEF_N> &a);
    void TwistIFFTlvl1(array<double, DEF_N> &res,const array<uint32_t,DEF_N> &a);
    void PolyMullvl1(array<uint32_t,DEF_N> &res,const array<uint32_t,DEF_N> &a, const array<uint32_t,DEF_N> &b);
}