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
    PolynomialInFDlvl1 ffta;
    TwistIFFTlvl1(ffta, a);
    PolynomialInFDlvl1 fftb;
    TwistIFFTlvl1(fftb, b);
    MulInFD<DEF_N>(ffta, ffta, fftb);
    TwistFFTlvl1(res, ffta);
}

inline void TwistFFTlvl2(array<uint64_t, DEF_nbar> &res,
                         const array<double, DEF_nbar> &a)
{
    fftplvl2.execute_direct_torus64(res.data(), a.data());
}

inline void TwistIFFTlvl2(PolynomialInFDlvl2 &res, const Polynomiallvl2 &a)
{
    fftplvl2.execute_reverse_torus64(res.data(), a.data());
}

inline void PolyMullvl2(Polynomiallvl2 &res, const Polynomiallvl2 &a,
                        const Polynomiallvl2 &b)
{
    array<double, DEF_nbar> ffta;
    TwistIFFTlvl2(ffta, a);
    array<double, DEF_nbar> fftb;
    TwistIFFTlvl2(fftb, b);
    MulInFD<DEF_nbar>(ffta, ffta, fftb);
    TwistFFTlvl2(res, ffta);
}
}  // namespace TFHEpp