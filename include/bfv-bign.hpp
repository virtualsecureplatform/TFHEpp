#pragma once

#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

#include "gatebootstrapping.hpp"
#include "homdecomp.hpp"
#include "mulfft.hpp"
#include "params.hpp"
#include "trlwe.hpp"
#include "trgsw.hpp"
#include "utils.hpp"

namespace TFHEpp {

template <class P>
inline void TwistFFTrescaleBigNum(Polynomial<P> &res,
                                  const PolynomialInFD<P> &a)
{
    const double q = std::ldexp(1.0, std::numeric_limits<typename P::T>::digits - 1);
    if constexpr (std::is_same_v<typename P::T, uint32_t>)
        fftplvl1.execute_direct_torus32_rescale_bignum(
            res.data(), a.data(), q, P::plain_modulus);
    else if constexpr (std::is_same_v<typename P::T, uint64_t>)
        fftplvl2.execute_direct_torus64_rescale_bignum(
            res.data(), a.data(), P::plain_modulus);
    else
        static_assert(false_v<typename P::T>, "Undefined BigNum FFT rescale");
}

template <class P>
inline void PolyMulRescaleUnsignedBigNum(Polynomial<P> &res,
                                         const Polynomial<P> &a,
                                         const Polynomial<P> &b)
{
    PolynomialInFD<P> ffta, fftb;
    TwistIFFTUInt<P>(ffta, a);
    TwistIFFTUInt<P>(fftb, b);
    MulInFD<P::n>(ffta, fftb);
    TwistFFTrescaleBigNum<P>(res, ffta);
}

template <class P>
inline void TRLWEMultWithoutRelinerizationBigNum(TRLWE3<P> &res,
                                                 const TRLWE<P> &a,
                                                 const TRLWE<P> &b)
{
    PolynomialInFD<P> ffta, fftb, fftc;
    TwistIFFTUInt<P>(ffta, a[0]);
    TwistIFFTUInt<P>(fftb, b[1]);
    MulInFD<P::n>(fftc, ffta, fftb);
    TwistIFFTUInt<P>(ffta, a[1]);
    TwistIFFTUInt<P>(fftb, b[0]);
    FMAInFD<P::n>(fftc, ffta, fftb);
    TwistFFTrescaleBigNum<P>(res[0], fftc);

    PolyMulRescaleUnsignedBigNum<P>(res[1], a[1], b[1]);
    PolyMulRescaleUnsignedBigNum<P>(res[2], a[0], b[0]);
}

template <class P>
inline void BigNumMult(TRLWE<P> &res, const TRLWE<P> &a, const TRLWE<P> &b,
                       const relinKeyFFT<P> &relinkeyfft)
{
    TRLWE3<P> mult;
    TRLWEMultWithoutRelinerizationBigNum<P>(mult, a, b);
    Relinearization<P>(res, mult, relinkeyfft);
}

template <class P>
void TLWES2TRLWEIKS(TRLWE<typename P::targetP> &res,
                    const std::vector<TLWE<typename P::domainP>> &tlwes,
                    const TLWE2TRLWEIKSKey<P> &iksk)
{
    const int tlnum = static_cast<int>(tlwes.size());
    constexpr typename P::domainP::T prec_offset =
        typename P::domainP::T(1)
        << (std::numeric_limits<typename P::domainP::T>::digits -
            (1 + P::basebit * P::t));
    constexpr uint32_t mask = (uint32_t(1) << P::basebit) - 1;
    res = {};
    constexpr uint32_t domain_digit =
        std::numeric_limits<typename P::domainP::T>::digits;
    constexpr uint32_t target_digit =
        std::numeric_limits<typename P::targetP::T>::digits;

    if constexpr (domain_digit == target_digit) {
        for (int i = 0; i < tlnum; i++) res[P::targetP::k][i] = tlwes[i][P::domainP::n];
    }
    else if constexpr (domain_digit > target_digit) {
        for (int i = 0; i < tlnum; i++)
            res[P::targetP::k][i] =
                (tlwes[i][P::domainP::n] +
                 (typename P::domainP::T(1)
                  << (domain_digit - target_digit - 1))) >>
                (domain_digit - target_digit);
    }
    else {
        for (int i = 0; i < tlnum; i++)
            res[P::targetP::k][i] =
                tlwes[i][P::domainP::n] << (target_digit - domain_digit);
    }

    for (int m = 0; m < tlnum; m++) {
        for (int i = 0; i < P::domainP::n; i++) {
            const typename P::domainP::T aibar = tlwes[m][i] + prec_offset;
            for (int j = 0; j < P::t; j++) {
                const uint32_t aij =
                    (aibar >> (std::numeric_limits<typename P::domainP::T>::digits -
                               (j + 1) * P::basebit)) &
                    mask;
                if (aij == 0) continue;
                for (int l = 0; l < P::targetP::k + 1; l++) {
                    for (uint64_t k = 0; k < P::targetP::n; k++) {
                        if (k + m < P::targetP::n)
                            res[l][k + m] -= iksk[i][j][aij - 1][l][k];
                        else
                            res[l][k + m - P::targetP::n] +=
                                iksk[i][j][aij - 1][l][k];
                    }
                }
            }
        }
    }
}

template <class P>
inline typename P::T tlweSymIntDecryptpra(const TLWE<P> &c, const Key<P> &key,
                                          const double delta)
{
    return static_cast<typename P::T>(std::round(tlweSymPhase<P>(c, key) / delta));
}

template <class P>
inline Polynomial<P> trlweSymDecryptpra(const TRLWE<P> &c, const Key<P> &key)
{
    return trlwePhase<P>(c, key);
}

template <class iksP, class bkP, class sskP, int num_multi, int shift = 0>
void TLWES2BigNumIKSezM(TRLWE<typename bkP::targetP> &res,
                       const std::vector<TLWE<typename iksP::domainP>> &tlwes,
                       const AnnihilateKey<typename bkP::targetP> &ahk,
                       const EvalKey &ek, const SecretKey &)
{
    constexpr int t = std::numeric_limits<typename bkP::targetP::T>::digits;
    const uint32_t tlnum = static_cast<uint32_t>(tlwes.size());
    std::vector<TLWE<typename sskP::domainP>> temp1(tlnum + t - 1);
    for (auto &item : temp1) item = {};

    constexpr int jcir = t / (num_multi * (shift + 1));
    std::array<Polynomial<typename bkP::targetP>, jcir> test_vectors{};
    int jpra = 0;
    for (int j = t - 1; j >= 0; j -= num_multi * (shift + 1)) {
        for (int k = 0; k < bkP::targetP::n; k += num_multi) {
            for (int l = 0; l < num_multi; l++) {
                test_vectors[jpra][k + l] = static_cast<typename bkP::targetP::T>(
                    std::ldexp(1.0,
                               std::numeric_limits<typename bkP::targetP::T>::digits -
                                   j + l * (shift + 1) - 3));
            }
        }
        jpra++;
    }

    for (uint32_t i = 0; i < tlnum; i++) {
        std::array<TLWE<typename bkP::targetP>, num_multi> temps{};
        jpra = 0;
        TLWE<typename iksP::targetP> tlwelvl0;
        IdentityKeySwitch<iksP>(tlwelvl0, tlwes[i], ek.getiksk<iksP>());
        for (int j = t - 1; j >= 0; j -= num_multi * (shift + 1)) {
            GateBootstrappingManyLUT<bkP, num_multi>(
                temps, tlwelvl0, ek.getbkfft<bkP>(), test_vectors[jpra]);
            for (int l = 0; l < num_multi; l++)
                temps[l][bkP::targetP::n] += test_vectors[jpra][l];

            for (int l = 0; l < num_multi; l++) {
                for (int m = 0; m < (shift + 1); m++) {
                    for (int k = 0; k < bkP::targetP::n + 1; k++) {
                        temp1[i + j - l * (shift + 1) - m][k] -=
                            temps[l][k] << m;
                    }
                }
            }
            jpra++;
        }
    }
    TLWE2TRLWEPacking<typename bkP::targetP>(res, temp1, ahk);
}

template <class iksP, class bkP, class sskP, int num_multi, int shift = 0>
inline void TLWES2BIGNUMIKSezM(TRLWE<typename bkP::targetP> &res,
                               const std::vector<TLWE<typename iksP::domainP>> &tlwes,
                               const AnnihilateKey<typename bkP::targetP> &ahk,
                               const EvalKey &ek, const SecretKey &sk)
{
    TLWES2BigNumIKSezM<iksP, bkP, sskP, num_multi, shift>(res, tlwes, ahk, ek,
                                                          sk);
}

template <class P>
std::array<typename P::T, P::n> generateDelbM(
    const std::array<typename P::T, P::n> &p)
{
    std::array<typename P::T, P::n> delb{};
    std::array<typename P::T, P::n> delbM{};
    const double q = std::ldexp(1.0, std::numeric_limits<typename P::T>::digits - 1);

    for (int i = 0; i < P::n; i++) delb[i] = -q / std::pow(P::plain_modulus, i + 1);
    for (int i = 0; i < P::n; i++) {
        for (int j = 0; j < P::n; j++) {
            if (i + j < P::n)
                delbM[i + j] = static_cast<typename P::T>(delbM[i + j] + delb[i] * p[j]);
            else
                delbM[i + j - P::n] = static_cast<typename P::T>(
                    delbM[i + j - P::n] - delb[i] * p[j]);
        }
    }
    return delbM;
}

template <class P>
inline TRLWE<P> bigNumSymIntEncrypt(const std::array<typename P::T, P::n> &p,
                                    const Key<P> &key)
{
    const auto hatp = generateDelbM<P>(p);
    TRLWE<P> c;
    trlweSymEncrypt<P>(c, hatp, key);
    return c;
}

template <class P>
std::array<int, P::n> bigNumSymIntDecrypt(const TRLWE<P> &c, const Key<P> &key)
{
    const auto phase = trlwePhase<P>(c, key);
    std::array<double, P::n> signed_phase{};
    std::array<int, P::n> digits{};
    const double q = std::ldexp(1.0, std::numeric_limits<typename P::T>::digits);
    const double q_half = std::ldexp(1.0, std::numeric_limits<typename P::T>::digits - 1);

    for (int i = 0; i < P::n; i++) {
        if (phase[i] > q_half)
            signed_phase[i] = static_cast<double>(phase[i]) - q;
        else
            signed_phase[i] = static_cast<double>(phase[i]);
    }

    for (int i = 0; i < P::n; i++) {
        const double value = (i == 0)
                                 ? -signed_phase[P::n - 1] - P::plain_modulus * signed_phase[0]
                                 : signed_phase[i - 1] - P::plain_modulus * signed_phase[i];
        digits[i] = static_cast<int>(std::round(value / q_half));
    }
    return digits;
}

template <class P>
inline Polynomial<P> EncodeHatEncoderP(const typename P::T &p)
{
    typename P::T value = p;
    Polynomial<P> ans{};
    for (int i = 0; i < P::n; i++) {
        ans[i] = value % P::plain_modulus;
        value -= ans[i];
        value /= P::plain_modulus;
    }
    return ans;
}

template <class P, int validbit = 32>
inline std::vector<uint8_t> EncodeHatEncoderInt8(const typename P::T &p)
{
    typename P::T value = p;
    std::vector<uint8_t> ans(validbit);
    for (int i = 0; i < validbit; i++) {
        ans[i] = value % P::plain_modulus;
        value -= ans[i];
        value /= P::plain_modulus;
    }
    return ans;
}

template <class P, int validbit = 32>
inline typename P::T decodeHatEncoderP(const Polynomial<P> &p)
{
    int64_t ans = 0;
    for (int i = P::n - 1; i >= 0; i--) {
        ans *= P::plain_modulus;
        ans += p[i];
    }
    return static_cast<typename P::T>(ans);
}

template <class P, int validbit>
inline typename P::T decodeHatEncoderInt82T(const std::vector<uint8_t> &p)
{
    typename P::T ans = 0;
    for (int i = validbit - 1; i >= 0; i--) {
        ans *= P::plain_modulus;
        ans += static_cast<int>(p[i]);
    }
    return ans;
}

template <class P, int validbit, int numcell>
inline typename P::T decodeHatEncoderInt2T(const std::array<int, numcell> &p)
{
    int64_t ans = 0;
    for (int i = validbit - 1; i >= 0; i--) {
        ans *= P::plain_modulus;
        ans += p[i];
    }
    return static_cast<typename P::T>(ans);
}

}  // namespace TFHEpp
