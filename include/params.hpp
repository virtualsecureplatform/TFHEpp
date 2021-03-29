#pragma once
#include <array>
#include <cmath>
#include <cstdint>

namespace TFHEpp {
using namespace std;

// Use old 80bit security parameters. It is faster, but not recommended.

#ifdef USE_80BIT_SECURITY
#include <params/CGGI16.hpp>
#elif defined(USE_CGGI19)
#include <params/CGGI19.hpp >
#else
#include <params/128bit.hpp>
#endif

struct lvl01param {
    using domainP = lvl0param;
    using targetP = lvl1param;
};
struct lvl02param {
    using domainP = lvl0param;
    using targetP = lvl2param;
};

struct lweParams {
    static constexpr lvl0param lvl0();
    static constexpr lvl1param lvl1();
    static constexpr lvl2param lvl2();
    static constexpr lvl01param lvl01();
    static constexpr lvl02param lvl02();
    static constexpr lvl10param lvl10();
    static constexpr lvl21param lvl21();
    static constexpr lvl22param lvl22();
};

template <class P>
using Key = array<typename P::T, P::n>;

template <class P>
using TLWE = array<typename P::T, P::n + 1>;

template <class P>
using Polynomial = array<typename P::T, P::n>;
template <class P>
using PolynomialInFD = array<double, P::n>;
template <class P>
using DecomposedPolynomial = Polynomial<P>;
template <class P>
using DecomposedPolynomialInFD = PolynomialInFD<P>;

template <class P>
using TRLWE = array<Polynomial<P>, 2>;
template <class P>
using TRLWEInFD = array<PolynomialInFD<P>, 2>;

template <class P>
using TRGSW = array<TRLWE<P>, 2 * P::l>;
template <class P>
using TRGSWFFT = array<TRLWEInFD<P>, 2 * P::l>;

template <class P>
using BootStrappingKeyFFT = array<TRGSWFFT<typename P::targetP>, P::domainP::n>;

template <class P>
using KeySwitchingKey =
    array<array<array<TLWE<typename P::targetP>, (1 << P::basebit) - 1>, P::t>,
          P::domainP::n>;
template <class P>
using PrivKeySwitchKey =
    array<array<array<TRLWE<typename P::targetP>, (1 << P::basebit) - 1>, P::t>,
          P::domainP::n + 1>;

}  // namespace TFHEpp