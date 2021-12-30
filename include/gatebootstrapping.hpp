#pragma once

#include <cmath>
#include <limits>

#include "cloudkey.hpp"
#include "detwfa.hpp"
#include "keyswitch.hpp"
#include "params.hpp"
#include "trlwe.hpp"
#include "utils.hpp"

namespace TFHEpp {

template <class P, uint32_t num_out = 1>
void BlindRotate(TRLWE<typename P::targetP> &res,
                 const TLWE<typename P::domainP> &tlwe,
                 const BootstrappingKeyFFT<P> &bkfft,
                 const Polynomial<typename P::targetP> &testvector)
{
    constexpr uint32_t bitwidth = bits_needed<num_out - 1>();
    const uint32_t b̄ = 2 * P::targetP::n -
                       ((tlwe[P::domainP::k * P::domainP::n] >>
                         (std::numeric_limits<typename P::domainP::T>::digits -
                          1 - P::targetP::nbit + bitwidth))
                        << bitwidth);
    res = {};
    PolynomialMulByXai<typename P::targetP>(res[P::targetP::k], testvector, b̄);
    for (int i = 0; i < P::domainP::k * P::domainP::n; i++) {
        constexpr typename P::domainP::T roundoffset =
            1ULL << (std::numeric_limits<typename P::domainP::T>::digits - 2 -
                     P::targetP::nbit + bitwidth);
        const uint32_t ā =
            (tlwe[i] + roundoffset) >>
            (std::numeric_limits<typename P::domainP::T>::digits - 1 -
             P::targetP::nbit + bitwidth)
                << bitwidth;
        if (ā == 0) continue;
        // Do not use CMUXFFT to avoid unnecessary copy.
        CMUXFFTwithPolynomialMulByXaiMinusOne<typename P::targetP>(res,
                                                                   bkfft[i], ā);
    }
}

template <class P, uint32_t num_out = 1>
void BlindRotate(TRLWE<typename P::targetP> &res,
                 const TLWE<typename P::domainP> &tlwe,
                 const BootstrappingKeyFFT<P> &bkfft,
                 const TRLWE<typename P::targetP> &testvector)
{
    constexpr uint32_t bitwidth = bits_needed<num_out - 1>();
    const uint32_t b̄ = 2 * P::targetP::n -
                       ((tlwe[P::domainP::k * P::domainP::n] >>
                         (std::numeric_limits<typename P::domainP::T>::digits -
                          1 - P::targetP::nbit + bitwidth))
                        << bitwidth);
    for (int k = 0; k < P::targetP::k + 1; k++)
        PolynomialMulByXai<typename P::targetP>(res[k], testvector[k], b̄);
    for (int i = 0; i < P::domainP::k * P::domainP::n; i++) {
        constexpr typename P::domainP::T roundoffset =
            1ULL << (std::numeric_limits<typename P::domainP::T>::digits - 2 -
                     P::targetP::nbit + bitwidth);
        const uint32_t ā =
            (tlwe[i] + roundoffset) >>
            (std::numeric_limits<typename P::domainP::T>::digits - 1 -
             P::targetP::nbit + bitwidth)
                << bitwidth;
        if (ā == 0) continue;
        // Do not use CMUXFFT to avoid unnecessary copy.
        CMUXFFTwithPolynomialMulByXaiMinusOne<typename P::targetP>(res,
                                                                   bkfft[i], ā);
    }
}

template <class P, uint32_t num_out = 1>
void BlindRotate(TRLWE<typename P::targetP> &res,
                 const TLWE<typename P::domainP> &tlwe,
                 const BootstrappingKeyNTT<P> &bkntt,
                 const Polynomial<typename P::targetP> &testvector)
{
    constexpr uint32_t bitwidth = bits_needed<num_out - 1>();
    const uint32_t b̄ = 2 * P::targetP::n -
                       ((tlwe[P::domainP::k * P::domainP::n] >>
                         (std::numeric_limits<typename P::domainP::T>::digits -
                          1 - P::targetP::nbit + bitwidth))
                        << bitwidth);
    res = {};
    PolynomialMulByXai<typename P::targetP>(res[P::targetP::k], testvector, b̄);
    for (int i = 0; i < P::domainP::k * P::domainP::n; i++) {
        constexpr typename P::domainP::T roundoffset =
            1ULL << (std::numeric_limits<typename P::domainP::T>::digits - 2 -
                     P::targetP::nbit + bitwidth);
        const uint32_t ā =
            (tlwe[i] + roundoffset) >>
            (std::numeric_limits<typename P::domainP::T>::digits - 1 -
             P::targetP::nbit + bitwidth)
                << bitwidth;
        if (ā == 0) continue;
        // Do not use CMUXNTT to avoid unnecessary copy.
        CMUXNTTwithPolynomialMulByXaiMinusOne<typename P::targetP>(res,
                                                                   bkntt[i], ā);
    }
}

template <class P>
void GateBootstrappingTLWE2TLWEFFT(
    TLWE<typename P::targetP> &res, const TLWE<typename P::domainP> &tlwe,
    const BootstrappingKeyFFT<P> &bkfft,
    const Polynomial<typename P::targetP> &testvector);

template <class P>
void GateBootstrappingTLWE2TLWENTT(
    TLWE<typename P::targetP> &res, const TLWE<typename P::domainP> &tlwe,
    const BootstrappingKeyNTT<P> &bkntt,
    const Polynomial<typename P::targetP> &testvector){
        TRLWE<typename P::targetP> acc;
        BlindRotate<P>(acc, tlwe, bkntt, testvector);
        SampleExtractIndex<typename P::targetP>(res, acc, 0);
    }

template <class P, uint32_t num_out>
void GateBootstrappingManyLUT(
    std::array<TLWE<typename P::targetP>, num_out> &res,
    const TLWE<typename P::domainP> &tlwe, const BootstrappingKeyFFT<P> &bkfft,
    const Polynomial<typename P::targetP> &testvector)
{
    TRLWE<typename P::targetP> acc;
    BlindRotate<P, num_out>(acc, tlwe, bkfft, testvector);
    for (int i = 0; i < num_out; i++)
        SampleExtractIndex<typename P::targetP>(res[i], acc, i);
}

template <class P, typename P::T μ>
constexpr Polynomial<P> μpolygen()
{
    Polynomial<P> poly;
    for (typename P::T &p : poly) p = μ;
    return poly;
}

template <typename lvl1param::T μ = lvl1param::μ>
void GateBootstrapping(TLWE<lvl0param> &res, const TLWE<lvl0param> &tlwe,
                       const EvalKey &ek)
{
    TLWE<lvl1param> tlwelvl1;
    GateBootstrappingTLWE2TLWEFFT<lvl01param>(tlwelvl1, tlwe, *ek.bkfftlvl01,
                                              μpolygen<lvl1param, μ>());
    IdentityKeySwitch<lvl10param>(res, tlwelvl1, *ek.iksklvl10);
}

template <typename lvl1param::T μ = lvl1param::μ>
void GateBootstrapping(TLWE<lvl1param> &res, const TLWE<lvl1param> &tlwe,
                       const EvalKey &ek)
{
    TLWE<lvl0param> tlwelvl0;
    IdentityKeySwitch<lvl10param>(tlwelvl0, tlwe, *ek.iksklvl10);
    GateBootstrappingTLWE2TLWEFFT<lvl01param>(res, tlwelvl0, *ek.bkfftlvl01,
                                              μpolygen<lvl1param, μ>());
}

template <typename lvl1param::T μ = lvl1param::μ>
void GateBootstrappingNTT(TLWE<lvl0param> &res, const TLWE<lvl0param> &tlwe,
                       const EvalKey &ek)
{
    TLWE<lvl1param> tlwelvl1;
    GateBootstrappingTLWE2TLWENTT<lvl01param>(tlwelvl1, tlwe, *ek.bknttlvl01,
                                              μpolygen<lvl1param, μ>());
    IdentityKeySwitch<lvl10param>(res, tlwelvl1, *ek.iksklvl10);
}

template <typename lvl1param::T μ = lvl1param::μ>
void GateBootstrappingNTT(TLWE<lvl1param> &res, const TLWE<lvl1param> &tlwe,
                       const EvalKey &ek)
{
    TLWE<lvl0param> tlwelvl0;
    IdentityKeySwitch<lvl10param>(tlwelvl0, tlwe, *ek.iksklvl10);
    GateBootstrappingTLWE2TLWENTT<lvl01param>(res, tlwelvl0, *ek.bknttlvl01,
                                              μpolygen<lvl1param, μ>());
}
}  // namespace TFHEpp