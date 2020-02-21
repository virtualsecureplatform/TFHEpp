#pragma once

#include <spqlios-fft.h>

#include <params.hpp>
#include <utils.hpp>

namespace TFHEpp {
using namespace std;

inline void TwistFFTlvl1(Polynomiallvl1 &res, const PolynomialInFDlvl1 &a)
{
    fftplvl1.execute_direct_torus32(res.data(), a.data());
}

inline void TwistIFFTlvl1(array<double, DEF_N> &res,
                          const array<uint32_t, DEF_N> &a)
{
    fftplvl1.execute_reverse_torus32(res.data(), a.data());
}

inline void PolyMullvl1(Polynomiallvl1 &res, const Polynomiallvl1 &a,
                        const Polynomiallvl1 &b)
{
    alignas(32) PolynomialInFDlvl1 ffta;
    TwistIFFTlvl1(ffta, a);
    alignas(32) PolynomialInFDlvl1 fftb;
    TwistIFFTlvl1(fftb, b);
    MulInFD<DEF_N>(ffta, ffta, fftb);
    TwistFFTlvl1(res, ffta);
}

inline void TwistFFTlvl2(Polynomiallvl2 &res, const PolynomialInFDlvl2 &a)
{
    fftplvl2.execute_direct_torus64(res.data(), a.data());
}

inline void TwistIFFTlvl2(PolynomialInFDlvl2 &res, const Polynomiallvl2 &a)
{
    fftplvl2.execute_reverse_torus64(res.data(), a.data());
}

inline void PolyMulNaievelvl2(Polynomiallvl2 &res, const Polynomiallvl2 &a,
                              const Polynomiallvl2 &b)
{
    for (int i = 0; i < DEF_nbar; i++) {
        uint64_t ri = 0;
        for (int j = 0; j <= i; j++)
            ri += static_cast<int64_t>(a[j]) * b[i - j];
        for (int j = i + 1; j < DEF_nbar; j++)
            ri -= static_cast<int64_t>(a[j]) * b[DEF_nbar + i - j];
        res[i] = ri;
    }
}
}  // namespace TFHEpp