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
std::array<int, P::n> clpxSymIntDecrypt(const TRLWE<P> &c, const Key<P> &key);

template <class P, int validbit, int numcell>
inline typename P::T decodeHatEncoderInt2T(const std::array<int, numcell> &p);

template <class P>
inline void TwistFFTrescaleCLPX(Polynomial<P> &res,
                                  const PolynomialInFD<P> &a)
{
    const double q = std::ldexp(1.0, std::numeric_limits<typename P::T>::digits - 1);
    if constexpr (std::is_same_v<typename P::T, uint32_t>)
        fftplvl1.execute_direct_torus32_rescale_clpx(
            res.data(), a.data(), q, P::plain_modulus);
    else if constexpr (std::is_same_v<typename P::T, uint64_t>)
        fftplvl2.execute_direct_torus64_rescale_clpx(
            res.data(), a.data(), P::plain_modulus);
    else
        static_assert(false_v<typename P::T>, "Undefined CLPX FFT rescale");
}

template <class P>
inline void PolyMulRescaleUnsignedCLPX(Polynomial<P> &res,
                                         const Polynomial<P> &a,
                                         const Polynomial<P> &b)
{
    PolynomialInFD<P> ffta, fftb;
    TwistIFFTUInt<P>(ffta, a);
    TwistIFFTUInt<P>(fftb, b);
    MulInFD<P::n>(ffta, fftb);
    TwistFFTrescaleCLPX<P>(res, ffta);
}

template <class P>
inline void TRLWEMultWithoutRelinerizationCLPX(TRLWE3<P> &res,
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
    TwistFFTrescaleCLPX<P>(res[0], fftc);

    PolyMulRescaleUnsignedCLPX<P>(res[1], a[1], b[1]);
    PolyMulRescaleUnsignedCLPX<P>(res[2], a[0], b[0]);
}

template <class P>
inline void CLPXMult(TRLWE<P> &res, const TRLWE<P> &a, const TRLWE<P> &b,
                       const relinKeyFFT<P> &relinkeyfft)
{
    TRLWE3<P> mult;
    TRLWEMultWithoutRelinerizationCLPX<P>(mult, a, b);
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
void TLWES2CLPXIKSezM(TRLWE<typename bkP::targetP> &res,
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

template <class P>
inline void GateBootstrappingTLWE2TLWEFFT(
    TLWE<typename P::targetP> &res, const TLWE<typename P::domainP> &tlwe,
    const BootstrappingKeyFFT<P> &bkfft,
    const Polynomial<typename P::targetP> &testvector)
{
    GateBootstrappingTLWE2TLWE<P>(res, tlwe, bkfft, testvector);
}

template <class P, uint32_t num_out>
inline void GateBootstrappingTLWE2TLWEFFTManyLut(
    std::array<TLWE<typename P::targetP>, num_out> &res,
    const TLWE<typename P::domainP> &tlwe, const BootstrappingKeyFFT<P> &bkfft,
    const Polynomial<typename P::targetP> &testvector)
{
    GateBootstrappingManyLUT<P, num_out>(res, tlwe, bkfft, testvector);
}

template <class iksP10, class iksP21, class bkP01, class bkP02, class iksP20,
          uint numdigit, uint basebit>
void CLPX2TLWESIKSAnyBit(
    std::vector<TLWE<typename bkP01::targetP>> &res,
    const TRLWE<typename iksP20::domainP> &trlwe, const EvalKey &ek,
    const SecretKey &sk)
{
    uint32_t batch_num = 16;
    const uint32_t tlnum = res.size();
    const uint32_t epoch_num = (tlnum + batch_num - 1) / batch_num;
    TLWE<typename bkP02::targetP> sumpra = {};
    std::array<TLWE<typename bkP02::targetP>, 2> temps;
    Polynomial<typename bkP02::targetP> testVectora = {};
    Polynomial<typename bkP02::targetP> testVectorb = {};

    for (int i = 0; i < bkP02::targetP::n; i++) {
        testVectora[i] = static_cast<typename bkP02::targetP::T>(
            static_cast<double>(i)
            * static_cast<double>(typename bkP02::targetP::T(1)
                                  << (std::numeric_limits<
                                          typename bkP02::targetP::T>::digits -
                                      (bkP02::targetP::nbit + 1))));
    }

    for (int k = 0; k < bkP02::targetP::n; k += 2) {
        for (int l = 0; l < 2; l++) {
            testVectorb[k + l] = static_cast<typename bkP02::targetP::T>(
                static_cast<double>(k)
                * static_cast<double>(typename bkP02::targetP::T(1)
                                      << (std::numeric_limits<
                                              typename bkP02::targetP::T>::digits -
                                          (bkP02::targetP::nbit + 2 + l))));
        }
    }

    for (uint32_t epoch = 0; epoch < epoch_num; epoch++) {
        if ((tlnum - epoch * batch_num) < batch_num)
            batch_num = tlnum - epoch * batch_num;

        for (uint32_t i = 0; i < batch_num; i++) {
            TLWE<typename iksP20::targetP> temp10 = {};

            if (epoch * batch_num + i > 0) {
                TLWE<typename bkP02::targetP> temp3 = {};
                TLWE<typename bkP02::targetP> temp31 = {};
                TLWE<typename bkP02::domainP> temp = {};

                for (int j = 0; j < bkP02::targetP::n + 1; j++) temp3[j] = temps[1][j];
                SampleExtractIndex<typename iksP20::domainP>(temp31, trlwe,
                                                             epoch * batch_num + i);

                IdentityKeySwitch<iksP20>(temp, temp31, ek.getiksk<iksP20>());
                GateBootstrappingTLWE2TLWEFFT<bkP02>(temp31, temp,
                                                     ek.getbkfft<bkP02>(),
                                                     testVectora);

                IdentityKeySwitch<iksP20>(temp, temp31, ek.getiksk<iksP20>());
                GateBootstrappingTLWE2TLWEFFTManyLut<bkP02, 2>(
                    temps, temp, ek.getbkfft<bkP02>(), testVectorb);

                for (int j = 0; j < bkP02::targetP::n + 1; j++)
                    temp31[j] = temp3[j] - temps[0][j];
                IdentityKeySwitch<iksP20>(temp10, temp31, ek.getiksk<iksP20>());
            }
            else {
                TLWE<typename bkP02::targetP> temp31 = {};
                TLWE<typename bkP02::domainP> temp = {};

                SampleExtractIndex<typename iksP20::domainP>(temp31, trlwe,
                                                             epoch * batch_num + i);

                IdentityKeySwitch<iksP20>(temp, temp31, ek.getiksk<iksP20>());
                GateBootstrappingTLWE2TLWEFFT<bkP02>(temp31, temp,
                                                     ek.getbkfft<bkP02>(),
                                                     testVectora);

                IdentityKeySwitch<iksP20>(temp, temp31, ek.getiksk<iksP20>());
                GateBootstrappingTLWE2TLWEFFTManyLut<bkP02, 2>(
                    temps, temp, ek.getbkfft<bkP02>(), testVectorb);

                for (int j = 0; j < bkP02::targetP::n + 1; j++) temp31[j] = -temps[0][j];
                IdentityKeySwitch<iksP20>(temp10, temp31, ek.getiksk<iksP20>());
            }

            std::array<Polynomial<typename bkP02::targetP>, 3> testVectors = {};
            TLWE<typename bkP01::domainP> temp1 = {};
            TLWE<typename bkP02::targetP> temp3 = {};
            TLWE<typename bkP02::targetP> temp31 = {};
            std::array<TLWE<typename bkP02::targetP>, 2> temp32;

            for (int k = 0; k < 2; k++)
                for (int j = 0; j < bkP02::targetP::n; j++)
                    testVectors[k][j] = static_cast<typename bkP02::targetP::T>(
                        static_cast<double>(typename bkP02::targetP::T(1)
                                            << (std::numeric_limits<
                                                    typename bkP02::targetP::T>::digits -
                                                batch_num + i - k - 3)));

            for (int k = 0; k < bkP02::targetP::n; k += 2) {
                testVectors[2][k] = static_cast<typename bkP02::targetP::T>(
                    static_cast<double>(typename bkP02::targetP::T(1)
                                        << (std::numeric_limits<
                                                typename bkP02::targetP::T>::digits -
                                            batch_num + i - 5)));
                testVectors[2][k + 1] = static_cast<typename bkP02::targetP::T>(
                    static_cast<double>(typename bkP02::targetP::T(1)
                                        << (std::numeric_limits<
                                                typename bkP02::targetP::T>::digits -
                                            4)));
            }

            for (int j = 0; j < bkP01::domainP::n + 1; j++) temp1[j] = temp10[j] << 2;
            temp1[iksP20::targetP::n] += static_cast<typename iksP20::targetP::T>(
                typename iksP20::targetP::T(1)
                << (std::numeric_limits<typename iksP20::targetP::T>::digits - 3));

            GateBootstrappingTLWE2TLWEFFTManyLut<bkP02, 2>(
                temp32, temp1, ek.getbkfft<bkP02>(), testVectors[2]);
            IdentityKeySwitch<iksP20>(temp1, temp32[1], ek.getiksk<iksP20>());
            for (int j = 0; j < iksP20::targetP::n + 1; j++) temp10[j] += temp1[j];
            temp10[iksP20::targetP::n] += static_cast<typename iksP20::targetP::T>(
                typename iksP20::targetP::T(1)
                << (std::numeric_limits<typename iksP20::targetP::T>::digits - 4));

            for (int j = 0; j < bkP01::domainP::n + 1; j++) temp1[j] = temp10[j] << 1;

            GateBootstrappingTLWE2TLWEFFT<bkP02>(temp31, temp1,
                                                 ek.getbkfft<bkP02>(),
                                                 testVectors[1]);
            GateBootstrappingTLWE2TLWEFFT<bkP02>(temp3, temp10,
                                                 ek.getbkfft<bkP02>(),
                                                 testVectors[0]);
            temp3[bkP02::targetP::n] -= static_cast<typename bkP02::targetP::T>(
                typename bkP02::targetP::T(1)
                << (std::numeric_limits<typename bkP02::targetP::T>::digits -
                    batch_num + i - 5));

            for (int j = 0; j < bkP02::targetP::n + 1; j++)
                sumpra[j] = temp3[j] - temp31[j] - temp32[0][j] + sumpra[j];
        }

        sumpra[bkP02::targetP::n] += static_cast<typename bkP02::targetP::T>(
            typename bkP02::targetP::T(1)
            << (std::numeric_limits<typename bkP02::targetP::T>::digits -
                batch_num - 6));

        std::array<TLWE<typename bkP01::targetP>, numdigit + 2> sums;
        HomDecomp<iksP21, iksP10, bkP01, basebit, numdigit + 2>(
            sums, sumpra, ek.getiksk<iksP21>(), ek.getiksk<iksP10>(),
            ek.getbkfft<bkP01>());

        Polynomial<typename bkP01::targetP> testVector3 = {};
        for (int i = 0; i < bkP01::targetP::n; i++)
            testVector3[i] = -static_cast<typename bkP01::targetP::T>(
                typename bkP01::targetP::T(1)
                << (std::numeric_limits<typename bkP01::targetP::T>::digits - 3));

        for (int i = 1; i < static_cast<int>(numdigit); i++) {
            TLWE<typename iksP20::targetP> temp1 = {};
            TLWE<typename iksP20::targetP> temp2 = {};

            sums[i][iksP10::domainP::n] += static_cast<typename iksP10::domainP::T>(
                typename iksP10::domainP::T(1)
                << (std::numeric_limits<typename iksP10::domainP::T>::digits - 3));
            IdentityKeySwitch<iksP10>(temp1, sums[i], ek.getiksk<iksP10>());

            for (int k = 0; k < static_cast<int>(basebit) - 1; k++) {
                for (int j = 0; j < iksP10::targetP::n + 1; j++)
                    temp2[j] = temp1[j] << (basebit - k - 1);

                GateBootstrappingTLWE2TLWEFFT<bkP01>(
                    res[epoch * batch_num + (i - 1) * basebit + k], temp2,
                    ek.getbkfft<bkP01>(), testVector3);
            }

            GateBootstrappingTLWE2TLWEFFT<bkP01>(
                res[epoch * batch_num + i * basebit - 1], temp1,
                ek.getbkfft<bkP01>(), testVector3);
        }

        if (epoch < epoch_num - 1) {
            std::array<Polynomial<typename bkP02::targetP>, 3> testVectors = {};
            TLWE<typename iksP20::targetP> temp1 = {};
            std::array<TLWE<typename iksP20::targetP>, 2> temp11;
            TLWE<typename bkP02::targetP> temp3 = {};
            TLWE<typename bkP02::targetP> temp31 = {};
            TLWE<typename bkP02::targetP> temp32 = {};

            for (int i = 0; i < 2; i++) {
                sums[numdigit + i][iksP10::domainP::n] +=
                    static_cast<typename iksP10::domainP::T>(
                        typename iksP10::domainP::T(1)
                        << (std::numeric_limits<typename iksP10::domainP::T>::digits - 3));
                IdentityKeySwitch<iksP10>(temp11[i], sums[numdigit + i],
                                          ek.getiksk<iksP10>());
            }

            for (int k = 0; k < 3; k++)
                for (int j = 0; j < bkP02::targetP::n; j++)
                    testVectors[k][j] = static_cast<typename bkP02::targetP::T>(
                        typename bkP02::targetP::T(1)
                        << (std::numeric_limits<typename bkP02::targetP::T>::digits -
                            batch_num - k - 3));

            GateBootstrappingTLWE2TLWEFFT<bkP02>(temp3, temp11[1],
                                                 ek.getbkfft<bkP02>(),
                                                 testVectors[0]);
            temp3[bkP02::targetP::n] -= static_cast<typename bkP02::targetP::T>(
                typename bkP02::targetP::T(1)
                << (std::numeric_limits<typename bkP02::targetP::T>::digits -
                    batch_num - 5));

            GateBootstrappingTLWE2TLWEFFT<bkP02>(temp31, temp11[0],
                                                 ek.getbkfft<bkP02>(),
                                                 testVectors[1]);

            for (int j = 0; j < iksP20::targetP::n + 1; j++) temp1[j] = temp11[0][j] << 1;
            GateBootstrappingTLWE2TLWEFFT<bkP02>(temp32, temp1,
                                                 ek.getbkfft<bkP02>(),
                                                 testVectors[2]);
            for (int j = 0; j < bkP02::targetP::n + 1; j++)
                sumpra[j] = temp3[j] - temp31[j] - temp32[j];
        }
    }
}

template <class iksP10, class iksP21, class bkP01, class bkP02, class iksP20,
          uint numdigit, uint basebit>
inline void CLPX2TLWESIKSanybit(
    std::vector<TLWE<typename bkP01::targetP>> &res,
    const TRLWE<typename iksP20::domainP> &trlwe, const EvalKey &ek,
    const SecretKey &sk)
{
    CLPX2TLWESIKSAnyBit<iksP10, iksP21, bkP01, bkP02, iksP20, numdigit,
                          basebit>(res, trlwe, ek, sk);
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
inline TRLWE<P> clpxSymIntEncrypt(const std::array<typename P::T, P::n> &p,
                                    const Key<P> &key)
{
    const auto hatp = generateDelbM<P>(p);
    TRLWE<P> c;
    trlweSymEncrypt<P>(c, hatp, key);
    return c;
}

template <class P>
std::array<int, P::n> clpxSymIntDecrypt(const TRLWE<P> &c, const Key<P> &key)
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
