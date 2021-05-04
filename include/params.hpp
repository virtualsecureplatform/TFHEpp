#pragma once
#include <array>
#include <cmath>
#include <cstdint>

namespace TFHEpp {
using namespace std;

// Use old 80bit security parameters. It is faster, but not recommended.

#ifdef USE_80BIT_SECURITY
#include "./params/CGGI16.hpp"
#elif defined(USE_CGGI19)
#include "./params/CGGI19.hpp"
#else
#include "./params/128bit.hpp"
#endif

struct lvl01param {
    using domainP = lvl0param;
    using targetP = lvl1param;
};
struct lvl02param {
    using domainP = lvl0param;
    using targetP = lvl2param;
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
using BootstrappingKey = array<TRGSW<typename P::targetP>, P::domainP::n>;
template <class P>
using BootstrappingKeyFFT = array<TRGSWFFT<typename P::targetP>, P::domainP::n>;

template <class P>
using KeySwitchingKey =
    array<array<array<TLWE<typename P::targetP>, (1 << P::basebit) - 1>, P::t>,
          P::domainP::n>;
template <class P>
using PrivKeySwitchKey =
    array<array<array<TRLWE<typename P::targetP>, (1 << P::basebit) - 1>, P::t>,
          P::domainP::n + 1>;

#define TFHEPP_EXPLICIT_INSTANTIATION_LVL0_1_2(fun) \
    fun(lvl0param);                            \
    fun(lvl1param);                            \
    fun(lvl2param);
#define TFHEPP_EXPLICIT_INSTANTIATION_LVL1_2(fun) \
    fun(lvl1param);                          \
    fun(lvl2param);
#define TFHEPP_EXPLICIT_INSTANTIATION_LVL01_02(fun) \
    fun(lvl01param);                           \
    fun(lvl02param);
#define TFHEPP_EXPLICIT_INSTANTIATION_KEY_SWITCH(fun) \
    fun(lvl10param);                              \
    fun(lvl20param);                              \
    fun(lvl21param);                              \
    fun(lvl22param);
#define TFHEPP_EXPLICIT_INSTANTIATION_LVL0221_0222(fun) \
    fun(lvl02param, lvl21param);                   \
    fun(lvl02param, lvl22param);
}  // namespace TFHEpp