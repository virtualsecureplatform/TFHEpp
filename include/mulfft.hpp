#pragma once

#include<array>

#include<spqlios-fft.h>

#include<params.hpp>

namespace TFHEpp{
    using namespace std;
    void TwistFFTlvl1(array<uint32_t, DEF_N> &res,const array<double,DEF_N> &a);
    void TwistIFFTlvl1(array<double, DEF_N> &res,const array<uint32_t,DEF_N> &a);
    void PolyMullvl1(array<uint32_t,DEF_N> &res,const array<uint32_t,DEF_N> &a, const array<uint32_t,DEF_N> &b);

    void TwistFFTlvl2(array<uint64_t, DEF_nbar> &res,const array<double,DEF_nbar> &a);
    void TwistIFFTlvl2(array<double, DEF_nbar> &res,const array<uint64_t,DEF_nbar> &a);
    void PolyMullvl2(array<uint64_t,DEF_nbar> &res,const array<uint64_t,DEF_nbar> &a, const array<uint64_t,DEF_nbar> &b);
}