#pragma once
#include <array>
#include <cmath>
#include <cstdint>

#include "cuhe++.hpp"
#include "raintt.hpp"

namespace TFHEpp {

template <class T, size_t N>
struct alignas(64) aligned_array : public std::array<T, N> {};

enum class ErrorDistribution { ModularGaussian, CenteredBinomial };

// Use old 80bit security parameters. It is faster, but not recommended.
#if defined(USE_80BIT_SECURITY)
#include "params/CGGI16.hpp"
#elif defined(USE_COMPRESS)
#include "params/compress.hpp"
#elif defined(USE_CGGI19)
#include "params/CGGI19.hpp"
#elif defined(USE_CONCRETE)
#include "params/concrete.hpp"
#elif defined(USE_TFHE_RS)
#include "params/tfhe-rs.hpp"
#elif defined(USE_TERNARY)
#include "params/ternary.hpp"
#else
#include "params/128bit.hpp"
#endif

struct lvl01param {
    using domainP = lvl0param;
    using targetP = lvl1param;
#ifdef USE_KEY_BUNDLE
    static constexpr uint32_t Addends = 2;
#else
    static constexpr uint32_t Addends = 1;
#endif
};

struct lvlh1param {
    using domainP = lvlhalfparam;
    using targetP = lvl1param;
#ifdef USE_KEY_BUNDLE
    static constexpr uint32_t Addends = 2;
#else
    static constexpr uint32_t Addends = 1;
#endif
};

struct lvl02param {
    using domainP = lvl0param;
    using targetP = lvl2param;
#ifdef USE_KEY_BUNDLE
    static constexpr uint32_t Addends = 2;
#else
    static constexpr uint32_t Addends = 1;
#endif
};

struct lvlh2param {
    using domainP = lvlhalfparam;
    using targetP = lvl2param;
#ifdef USE_KEY_BUNDLE
    static constexpr uint32_t Addends = 2;
#else
    static constexpr uint32_t Addends = 1;
#endif
};

template <class P>
using Key = std::array<typename P::T, P::k * P::n>;

template <class P>
using TLWE = aligned_array<typename P::T, P::k * P::n + 1>;

template <class P>
using Polynomial = std::array<typename P::T, P::n>;
template <class P>
using PolynomialInFD = std::array<double, P::n>;
template <class P>
using PolynomialNTT = std::array<cuHEpp::INTorus, P::n>;
template <class P>
using PolynomialRAINTT = std::array<raintt::DoubleSWord, P::n>;
template <class P>
using DecomposedPolynomial = std::array<Polynomial<P>, P::l>;
template <class P>
using DecomposedPolynomialNTT = std::array<PolynomialNTT<P>, P::l>;
template <class P>
using DecomposedPolynomialRAINTT = std::array<PolynomialRAINTT<P>, P::l>;

template <class P>
using TRLWE = std::array<Polynomial<P>, P::k + 1>;
template <class P>
using TRLWE3 = std::array<Polynomial<P>, 3>;
template <class P>
using TRLWEInFD = std::array<PolynomialInFD<P>, P::k + 1>;
template <class P>
using TRLWENTT = std::array<PolynomialNTT<P>, P::k + 1>;
template <class P>
using TRLWERAINTT = std::array<PolynomialRAINTT<P>, P::k + 1>;

template <class P>
using TRGSW = std::array<TRLWE<P>, (P::k + 1) * P::l>;
template <class P>
using HalfTRGSW = std::array<TRLWE<P>, P::l>;
template <class P>
using TRGSWFFT = aligned_array<TRLWEInFD<P>, (P::k + 1) * P::l>;
template <class P>
using HalfTRGSWFFT = aligned_array<TRLWEInFD<P>, P::l>;
template <class P>
using TRGSWNTT = std::array<TRLWENTT<P>, (P::k + 1) * P::l>;
template <class P>
using TRGSWRAINTT = std::array<TRLWERAINTT<P>, (P::k + 1) * P::l>;

#ifdef USE_KEY_BUNDLE
template <class P>
using BootstrappingKeyElement =
    std::array<TRGSW<typename P::targetP>, (1 << P::Addends) - 1>;
template <class P>
using BootstrappingKeyElementFFT =
    std::array<TRGSWFFT<typename P::targetP>, 1 << P::Addends>;
#else
template <class P>
using BootstrappingKeyElement =
    std::array<TRGSW<typename P::targetP>, P::domainP::key_value_diff>;
template <class P>
using BootstrappingKeyElementFFT =
    std::array<TRGSWFFT<typename P::targetP>, P::domainP::key_value_diff>;
#endif

template <class P>
using BootstrappingKey =
    std::array<BootstrappingKeyElement<P>, P::domainP::k * P::domainP::n>;
template <class P>
using BootstrappingKeyFFT =
    std::array<BootstrappingKeyElementFFT<P>,
               P::domainP::k * P::domainP::n / P::Addends>;
template <class P>
using BootstrappingKeyNTT =
    std::array<TRGSWNTT<typename P::targetP>, P::domainP::k * P::domainP::n>;
template <class P>
using BootstrappingKeyRAINTT =
    std::array<TRGSWRAINTT<typename P::targetP>, P::domainP::k * P::domainP::n>;

template <class P>
using KeySwitchingKey = std::array<
    std::array<std::array<TLWE<typename P::targetP>, (1 << (P::basebit - 1))>,
               P::t>,
    P::domainP::k * P::domainP::n>;
template <class P>
using SubsetKeySwitchingKey = std::array<
    std::array<std::array<TLWE<typename P::targetP>, (1 << P::basebit) - 1>,
               P::t>,
    P::domainP::k * P::domainP::n - P::targetP::k * P::targetP::n>;
template <class P>
using TLWE2TRLWEIKSKey = std::array<
    std::array<std::array<TRLWE<typename P::targetP>, (1 << P::basebit) - 1>,
               P::t>,
    P::domainP::n>;
template <class P>
using EvalAutoKey = std::array<HalfTRGSWFFT<P>, P::k>;
template <class P>
using AnnihilateKey = std::array<EvalAutoKey<P>, P::nbit>;
template <class P>
using PrivateKeySwitchingKey = std::array<
    std::array<std::array<TRLWE<typename P::targetP>, (1 << (P::basebit - 1))>,
               P::t>,
    P::domainP::k * P::domainP::n + 1>;
template <class P>
using SubsetPrivateKeySwitchingKey = std::array<
    std::array<std::array<TRLWE<typename P::targetP>, (1 << P::basebit) - 1>,
               P::t>,
    P::targetP::k * P::targetP::n + 1>;
template <class P>
using relinKey = std::array<TRLWE<P>, P::l>;
template <class P>
using relinKeyFFT = std::array<TRLWEInFD<P>, P::l>;

#define TFHEPP_EXPLICIT_INSTANTIATION_TLWE(fun) \
    fun(lvl0param);                             \
    fun(lvlhalfparam);                          \
    fun(lvl1param);                             \
    fun(lvl2param);                             \
    fun(lvl3param);
#define TFHEPP_EXPLICIT_INSTANTIATION_TRLWE(fun) \
    fun(lvl1param);                              \
    fun(lvl2param);
#define TFHEPP_EXPLICIT_INSTANTIATION_BLIND_ROTATE(fun) \
    fun(lvl01param);                                    \
    fun(lvl02param);
#define TFHEPP_EXPLICIT_INSTANTIATION_KEY_SWITCH_TO_TLWE(fun) \
    fun(lvl10param);                                          \
    fun(lvl1hparam);                                          \
    fun(lvl20param);                                          \
    fun(lvl2hparam);                                          \
    fun(lvl21param);
#define TFHEPP_EXPLICIT_INSTANTIATION_KEY_SWITCH_TO_TRLWE(fun) \
    fun(lvl11param);                                           \
    fun(lvl21param);                                           \
    fun(lvl22param);
#define TFHEPP_EXPLICIT_INSTANTIATION_SUBSET_KEY_SWITCH_TO_TLWE(fun) \
    fun(lvl21param);
#define TFHEPP_EXPLICIT_INSTANTIATION_SUBSET_KEY_SWITCH_TO_TRLWE(fun) \
    fun(lvl21param);
#define TFHEPP_EXPLICIT_INSTANTIATION_GATE_IKSBR(fun) \
    fun(lvl10param, lvl01param, lvl1param::μ);
#define TFHEPP_EXPLICIT_INSTANTIATION_GATE_BRIKS(fun) \
    fun(lvl01param, lvl1param::μ, lvl10param);
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
#define TFHEPP_EXPLICIT_INSTANTIATION_CIRCUIT_BOOTSTRAPPING_SUBIKS(fun) \
    fun(lvl10param, lvl02param, lvl21param);
}  // namespace TFHEpp