#pragma once

#include <cmath>
#include <limits>

#include "cloudkey.hpp"
#include "detwfa.hpp"
#include "keyswitch.hpp"
#include "params.hpp"
#include "trlwe.hpp"
#include "utils.hpp"

#ifdef USE_KEY_BUNDLE
#include "keybundle.hpp"
#endif

namespace TFHEpp {

// https://eprint.iacr.org/2025/809
template <class P, uint32_t num_out>
void BRModSwitch(ModswitchTLWE<typename P::domainP> &moded,
                 const TLWE<typename P::domainP> &tlwe)
{
    constexpr uint32_t bitwidth = bits_needed<num_out - 1>();
    std::make_signed_t<typename P::domainP::T> c = 0;
    constexpr typename P::domainP::T roundoffset =
        1ULL << (std::numeric_limits<typename P::domainP::T>::digits - 2 -
                 P::targetP::nbit + bitwidth);
    for (int i = 0; i < P::domainP::k * P::domainP::n; i++) {
        moded[i] = (tlwe[i] + roundoffset) >>
                   (std::numeric_limits<typename P::domainP::T>::digits - 1 -
                    P::targetP::nbit + bitwidth)
                       << bitwidth;
        c += tlwe[i] -
             (moded[i] << (std::numeric_limits<typename P::domainP::T>::digits -
                           1 - P::targetP::nbit));
    }
    moded[P::domainP::k * P::domainP::n] =
        2 * P::targetP::n -
        (static_cast<typename P::domainP::T>(
             tlwe[P::domainP::k * P::domainP::n] - c / 2 + roundoffset) >>
         (std::numeric_limits<typename P::domainP::T>::digits - 1 -
          P::targetP::nbit + bitwidth)
             << bitwidth);
}

template <class P, uint32_t num_out = 1>
void BlindRotate(TRLWE<typename P::targetP> &res,
                 const TLWE<typename P::domainP> &tlwe,
                 const BootstrappingKeyFFT<P> &bkfft,
                 const Polynomial<typename P::targetP> &testvector)
{
    ModswitchTLWE<typename P::domainP> moded;
    BRModSwitch<P, num_out>(moded, tlwe);
    res = {};
    PolynomialMulByXai<typename P::targetP>(
        res[P::targetP::k], testvector, moded[P::domainP::k * P::domainP::n]);
#ifdef USE_KEY_BUNDLE
    for (int i = 0; i < P::domainP::k * P::domainP::n / P::Addends; i++) {
        alignas(64) TRGSWFFT<typename P::targetP> BKadded;
        KeyBundleFFT<P>(BKadded, bkfft[i],
                        std::span(moded)
                            .subspan(P::Addends * i, P::Addends)
                            .template first<P::Addends>());
        ExternalProduct<typename P::targetP>(res, res, BKadded);
    }
#else
    for (int i = 0; i < P::domainP::k * P::domainP::n; i++) {
        if (moded[i] == 0) continue;
        // Do not use CMUXFFT to avoid unnecessary copy.
        CMUXwithPolynomialMulByXaiMinusOne<P>(res, bkfft[i], moded[i]);
    }
#endif
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
#ifdef USE_KEY_BUNDLE
    alignas(64) std::array<TRGSWFFT<typename P::targetP>,
                           P::domainP::k * P::domainP::n / P::Addends>
        BKadded;
#pragma omp parallel for num_threads(4)
    for (int i = 0; i < P::domainP::k * P::domainP::n / P::Addends; i++) {
        constexpr typename P::domainP::T roundoffset =
            1ULL << (std::numeric_limits<typename P::domainP::T>::digits - 2 -
                     P::targetP::nbit + bitwidth);
        std::array<typename P::domainP::T, P::Addends> bara;
        bara[0] = (tlwe[2 * i] + roundoffset) >>
                  (std::numeric_limits<typename P::domainP::T>::digits - 1 -
                   P::targetP::nbit + bitwidth)
                      << bitwidth;
        bara[1] = (tlwe[2 * i + 1] + roundoffset) >>
                  (std::numeric_limits<typename P::domainP::T>::digits - 1 -
                   P::targetP::nbit + bitwidth)
                      << bitwidth;
        KeyBundleFFT<P>(BKadded[i], bkfft[i], bara);
    }
    for (int i = 0; i < P::domainP::k * P::domainP::n / P::Addends; i++) {
        ExternalProduct<typename P::targetP>(res, res, BKadded[i]);
    }
#else
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
        CMUXwithPolynomialMulByXaiMinusOne<P>(res, bkfft[i], ā);
    }
#endif
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
        CMUXwithPolynomialMulByXaiMinusOne<typename P::targetP>(res,
                                                                   bkntt[i], ā);
    }
}

template <class P, uint32_t num_out = 1>
void BlindRotate(TRLWE<typename P::targetP> &res,
                 const TLWE<typename P::domainP> &tlwe,
                 const BootstrappingKeyRAINTT<P> &bkraintt,
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
        CMUXwithPolynomialMulByXaiMinusOne<typename P::targetP>(
            res, bkraintt[i], ā);
    }
}

template <class P>
void GateBootstrappingTLWE2TLWE(
    TLWE<typename P::targetP> &res, const TLWE<typename P::domainP> &tlwe,
    const BootstrappingKeyFFT<P> &bkfft,
    const Polynomial<typename P::targetP> &testvector)
{
    alignas(64) TRLWE<typename P::targetP> acc;
    BlindRotate<P>(acc, tlwe, bkfft, testvector);
    SampleExtractIndex<typename P::targetP>(res, acc, 0);
}

template <class P>
void GateBootstrappingTLWE2TLWE(TLWE<typename P::targetP> &res,
                                const TLWE<typename P::domainP> &tlwe,
                                const BootstrappingKeyFFT<P> &bkfft,
                                const TRLWE<typename P::targetP> &testvector)
{
    alignas(64) TRLWE<typename P::targetP> acc;
    BlindRotate<P>(acc, tlwe, bkfft, testvector);
    SampleExtractIndex<typename P::targetP>(res, acc, 0);
}

template <class P>
void GateBootstrappingTLWE2TLWE(
    TLWE<typename P::targetP> &res, const TLWE<typename P::domainP> &tlwe,
    const BootstrappingKeyNTT<P> &bkntt,
    const Polynomial<typename P::targetP> &testvector)
{
    alignas(64) TRLWE<typename P::targetP> acc;
    BlindRotate<P>(acc, tlwe, bkntt, testvector);
    SampleExtractIndex<typename P::targetP>(res, acc, 0);
}

template <class P>
void GateBootstrappingTLWE2TLWE(
    TLWE<typename P::targetP> &res, const TLWE<typename P::domainP> &tlwe,
    const BootstrappingKeyRAINTT<P> &bkraintt,
    const Polynomial<typename P::targetP> &testvector)
{
    alignas(64) TRLWE<typename P::targetP> acc;
    BlindRotate<P>(acc, tlwe, bkraintt, testvector);
    SampleExtractIndex<typename P::targetP>(res, acc, 0);
}

template <class P, uint32_t num_out>
void GateBootstrappingManyLUT(
    std::array<TLWE<typename P::targetP>, num_out> &res,
    const TLWE<typename P::domainP> &tlwe, const BootstrappingKeyFFT<P> &bkfft,
    const Polynomial<typename P::targetP> &testvector)
{
    alignas(64) TRLWE<typename P::targetP> acc;
    BlindRotate<P, num_out>(acc, tlwe, bkfft, testvector);
    for (int i = 0; i < num_out; i++)
        SampleExtractIndex<typename P::targetP>(res[i], acc, i);
}

template <class P, uint32_t num_out>
void GateBootstrappingManyLUT(
    std::array<TLWE<typename P::targetP>, num_out> &res,
    const TLWE<typename P::domainP> &tlwe, const BootstrappingKeyFFT<P> &bkfft,
    const TRLWE<typename P::targetP> &testvector)
{
    alignas(64) TRLWE<typename P::targetP> acc;
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

// Double Decomposition variants
template <class P, uint32_t num_out = 1>
void BlindRotateDD(TRLWE<typename P::targetP> &res,
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
        CMUXwithPolynomialMulByXaiMinusOneDD<P>(res, bkfft[i], ā);
    }
}

template <class P>
void GateBootstrappingTLWE2TLWEDD(
    TLWE<typename P::targetP> &res, const TLWE<typename P::domainP> &tlwe,
    const BootstrappingKeyFFT<P> &bkfft,
    const Polynomial<typename P::targetP> &testvector)
{
    alignas(64) TRLWE<typename P::targetP> acc;
    BlindRotateDD<P>(acc, tlwe, bkfft, testvector);
    SampleExtractIndex<typename P::targetP>(res, acc, 0);
}

template <class bkP, typename bkP::targetP::T μ, class iksP>
void GateBootstrapping(TLWE<typename iksP::targetP> &res,
                       const TLWE<typename bkP::domainP> &tlwe,
                       const EvalKey &ek)
{
    alignas(64) TLWE<typename bkP::targetP> tlwelvl1;
    GateBootstrappingTLWE2TLWE<bkP>(tlwelvl1, tlwe, ek.getbkfft<bkP>(),
                                    μpolygen<typename bkP::targetP, μ>());
    IdentityKeySwitch<iksP>(res, tlwelvl1, ek.getiksk<iksP>());
}

template <class iksP, class bkP, typename bkP::targetP::T μ>
void GateBootstrapping(TLWE<typename bkP::targetP> &res,
                       const TLWE<typename iksP::domainP> &tlwe,
                       const EvalKey &ek)
{
    alignas(64) TLWE<typename iksP::targetP> tlwelvl0;
    IdentityKeySwitch<iksP>(tlwelvl0, tlwe, ek.getiksk<iksP>());
    GateBootstrappingTLWE2TLWE<bkP>(res, tlwelvl0, ek.getbkfft<bkP>(),
                                    μpolygen<typename bkP::targetP, μ>());
}

template <class bkP, typename bkP::targetP::T μ, class iksP>
void GateBootstrappingNTT(TLWE<typename iksP::tagetP> &res,
                          const TLWE<typename bkP::domainP> &tlwe,
                          const EvalKey &ek)
{
    alignas(64) TLWE<typename bkP::targetP> tlwelvl1;
    GateBootstrappingTLWE2TLWE<bkP>(tlwelvl1, tlwe, ek.getbkntt<bkP>(),
                                    μpolygen<typename bkP::targetP, μ>());
    IdentityKeySwitch<iksP>(res, tlwelvl1, ek.getiksk<iksP>());
}

template <class iksP, class bkP, typename bkP::targetP::T μ>
void GateBootstrappingNTT(TLWE<typename bkP::targetP> &res,
                          const TLWE<typename iksP::domainP> &tlwe,
                          const EvalKey &ek)
{
    alignas(64) TLWE<typename iksP::targetP> tlwelvl0;
    IdentityKeySwitch<iksP>(tlwelvl0, tlwe, ek.getiksk<iksP>());
    GateBootstrappingTLWE2TLWE<bkP>(res, tlwelvl0, ek.getbkntt<bkP>(),
                                    μpolygen<typename bkP::targetP, μ>());
}

}  // namespace TFHEpp