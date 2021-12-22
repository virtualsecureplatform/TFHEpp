#pragma once
#include <array>
#include <cmath>
#include <cstdint>

#include "cuhe++.hpp"

namespace TFHEpp {

// Use old 80bit security parameters. It is faster, but not recommended.

#if defined(USE_80BIT_SECURITY)
#include "params/CGGI16.hpp"
#elif defined(USE_CGGI19)
#include "params/CGGI19.hpp"
#elif defined(USE_CONCRETE)
#include "params/concrete.hpp"
#else
#include "params/128bit.hpp"
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
using Key = std::array<typename P::T, P::k * P::n>;

template <class P>
using TLWE = std::array<typename P::T, P::k * P::n + 1>;

template <class P>
using Polynomial = std::array<typename P::T, P::n>;
template <class P>
using UnsignedPolynomial = Polynomial<P>;
template <class P>
using PolynomialInFD = std::array<double, P::n>;
template <class P>
using PolynomialNTT = std::array<cuHEpp::INTorus, P::n>;
template <class P>
using DecomposedPolynomial = Polynomial<P>;
template <class P>
using DecomposedPolynomialInFD = PolynomialInFD<P>;
template <class P>
using DecomposedPolynomialNTT = PolynomialNTT<P>;

template <class P>
using TRLWE = std::array<Polynomial<P>, P::k + 1>;
template <class P>
using UnsignedTRLWE = std::array<Polynomial<P>, P::k + 1>;
template <class P>
using TRLWE3 = std::array<Polynomial<P>, 3>;
template <class P>
using TRLWEInFD = std::array<PolynomialInFD<P>, P::k + 1>;
template <class P>
using TRLWENTT = std::array<PolynomialNTT<P>, P::k + 1>;

template <class P>
using TRGSW = std::array<TRLWE<P>, (P::k + 1) * P::l>;
template <class P>
using TRGSWFFT = std::array<TRLWEInFD<P>, (P::k + 1) * P::l>;
template <class P>
using TRGSWNTT = std::array<TRLWENTT<P>, (P::k + 1) * P::l>;

template <class P>
using BootstrappingKey = std::array<TRGSW<typename P::targetP>, P::domainP::n>;
template <class P>
using BootstrappingKeyFFT =
    std::array<TRGSWFFT<typename P::targetP>, P::domainP::n>;
template <class P>
using BootstrappingKeyNTT =
    std::array<TRGSWNTT<typename P::targetP>, P::domainP::n>;

template <class P>
using KeySwitchingKey = std::array<
    std::array<std::array<TLWE<typename P::targetP>, (1 << P::basebit) - 1>,
               P::t>,
    P::domainP::k * P::domainP::n>;
template <class P>
using TLWE2TRLWEIKSKey = std::array<
    std::array<std::array<TRLWE<typename P::targetP>, (1 << P::basebit) - 1>,
               P::t>,
    P::domainP::n>;
template <class P>
using AnnihilateKey = std::array<TRGSWFFT<P>, P::nbit>;
template <class P>
using PrivateKeySwitchingKey = std::array<
    std::array<std::array<TRLWE<typename P::targetP>, (1 << P::basebit) - 1>,
               P::t>,
    P::domainP::n + 1>;
template <class P>
using relinKey = std::array<TRLWE<P>, P::l>;
template <class P>
using relinKeyFFT = std::array<TRLWEInFD<P>, P::l>;

#define TFHEPP_EXPLICIT_INSTANTIATION_TLWE(fun) \
    fun(lvl0param);                             \
    fun(lvl1param);                             \
    fun(lvl2param);
#define TFHEPP_EXPLICIT_INSTANTIATION_TRLWE(fun) \
    fun(lvl1param);                              \
    fun(lvl2param);
#define TFHEPP_EXPLICIT_INSTANTIATION_BLIND_ROTATE(fun) \
    fun(lvl01param);                                    \
    fun(lvl02param);
#define TFHEPP_EXPLICIT_INSTANTIATION_KEY_SWITCH_TO_TLWE(fun) \
    fun(lvl10param);                                          \
    fun(lvl11param);                                          \
    fun(lvl20param);                                          \
    fun(lvl21param);                                          \
    fun(lvl22param);
#define TFHEPP_EXPLICIT_INSTANTIATION_KEY_SWITCH_TO_TRLWE(fun) \
    fun(lvl11param);                                           \
    fun(lvl21param);                                           \
    fun(lvl22param);
#define TFHEPP_EXPLICIT_INSTANTIATION_GATE(fun) \
    fun(lvl10param, lvl01param);                \
    fun(lvl10param, lvl02param);                \
    fun(lvl20param, lvl01param);                \
    fun(lvl20param, lvl02param);
#define TFHEPP_EXPLICIT_INSTANTIATION_CIRCUIT_KEY(fun) \
    fun(lvl02param, lvl21param);                       \
    fun(lvl02param, lvl22param);
#define TFHEPP_EXPLICIT_INSTANTIATION_CIRCUIT_BOOTSTRAPPING(fun) \
    fun(lvl10param, lvl02param, lvl21param);                     \
    fun(lvl10param, lvl02param, lvl22param);
}  // namespace TFHEpp