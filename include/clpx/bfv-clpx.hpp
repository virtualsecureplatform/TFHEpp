#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <type_traits>
#include <vector>

#include "tfhe/gatebootstrapping.hpp"
#include "tfhe/homdecomp.hpp"
#include "mulfft.hpp"
#include "params.hpp"
#include "tfhe/trlwe.hpp"
#include "tfhe/trgsw.hpp"
#include "utils.hpp"

namespace TFHEpp {

template <class T>
inline T clpxSignedDoubleToTorus(const double value)
{
    if constexpr (std::is_same_v<T, uint32_t>)
        return static_cast<uint32_t>(static_cast<int32_t>(std::lround(value)));
    else if constexpr (std::is_same_v<T, uint64_t>)
        return static_cast<uint64_t>(static_cast<int64_t>(std::llround(value)));
    else
        static_assert(false_v<T>, "Undefined CLPX torus conversion");
}

template <class P>
std::array<int, P::n> clpxSymIntDecrypt(const TRLWE<P> &c, const Key<P> &key);

template <class P, int validbit, int numcell>
inline typename P::T decodeHatEncoderInt2T(const std::array<int, numcell> &p);

template <class P, int validbit, int numcell>
inline unsigned __int128 decodeHatEncoderInt2U128(
    const std::array<int, numcell> &p);

template <class P, int validbit, int numcell>
inline __int128 decodeHatEncoderInt2I128(const std::array<int, numcell> &p);

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
    TwistIFFT<P>(ffta, a);
    TwistIFFT<P>(fftb, b);
    MulInFD<P::n>(ffta, fftb);
    TwistFFTrescaleCLPX<P>(res, ffta);
}

template <class P>
inline void TRLWEMultWithoutRelinerizationCLPX(TRLWE3<P> &res,
                                                 const TRLWE<P> &a,
                                                 const TRLWE<P> &b)
{
    PolynomialInFD<P> ffta, fftb, fftc;
    TwistIFFT<P>(ffta, a[0]);
    TwistIFFT<P>(fftb, b[1]);
    MulInFD<P::n>(fftc, ffta, fftb);
    TwistIFFT<P>(ffta, a[1]);
    TwistIFFT<P>(fftb, b[0]);
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

template <class iksP, class bkP, class sskP, int num_multi, int shift = 0,
          int w = std::numeric_limits<typename bkP::targetP::T>::digits>
void TLWES2CLPXIKS(TRLWE<typename bkP::targetP> &res,
                   const std::vector<TLWE<typename iksP::domainP>> &tlwes,
                   const AnnihilateKey<typename bkP::targetP> &ahk,
                   const EvalKey &ek)
{
    constexpr int t = std::numeric_limits<typename bkP::targetP::T>::digits;
    static_assert(w > 0, "TLWES2CLPXIKS requires a positive truncation width");
    static_assert(w <= t, "TLWES2CLPXIKS truncation width exceeds torus width");
    const uint32_t tlnum = static_cast<uint32_t>(tlwes.size());
    std::vector<TLWE<typename sskP::domainP>> temp1(tlnum + w - 1);
    for (auto &item : temp1) item = {};

    // Nagai et al. denote shift + 1 as shiftnum and use w as the
    // Delta_b truncation width. The old default keeps the full torus width.
    constexpr int shiftnum = shift + 1;
    constexpr int step = num_multi * shiftnum;
    constexpr int jcir = (w + step - 1) / step;
    std::array<Polynomial<typename bkP::targetP>, jcir> test_vectors{};
    int jpra = 0;
    for (int j = w - 1; j >= 0; j -= step) {
        for (int k = 0; k < bkP::targetP::n; k += num_multi) {
            for (int l = 0; l < num_multi; l++) {
                test_vectors[jpra][k + l] = static_cast<typename bkP::targetP::T>(
                    std::ldexp(1.0,
                               std::numeric_limits<typename bkP::targetP::T>::digits -
                                   j + l * shiftnum - 3));
            }
        }
        jpra++;
    }

    using BootstrapOutputs =
        std::array<TLWE<typename bkP::targetP>, num_multi>;
    std::vector<BootstrapOutputs> bootstrap_outputs(tlnum * jcir);
    const auto &iksk = ek.getiksk<iksP>();
    const auto &bkfft = ek.getbkfft<bkP>();

    // The key switch and blind rotation for one input TLWE do not depend on
    // any other input.  Keep their comparatively large outputs separate so
    // the overlapping diagonal accumulation into temp1 remains race-free.
#pragma omp parallel for schedule(static)
    for (int64_t i = 0; i < static_cast<int64_t>(tlnum); i++) {
        TLWE<typename iksP::targetP> tlwelvl0;
        IdentityKeySwitch<iksP>(tlwelvl0, tlwes[i], iksk);
        int local_jpra = 0;
        for (int j = w - 1; j >= 0; j -= step) {
            auto &temps = bootstrap_outputs[i * jcir + local_jpra];
            GateBootstrappingManyLUT<bkP, num_multi>(
                temps, tlwelvl0, bkfft, test_vectors[local_jpra]);
            for (int l = 0; l < num_multi; l++)
                temps[l][bkP::targetP::n] += test_vectors[local_jpra][l];

            local_jpra++;
        }
    }

    for (uint32_t i = 0; i < tlnum; i++) {
        jpra = 0;
        for (int j = w - 1; j >= 0; j -= step) {
            const auto &temps = bootstrap_outputs[i * jcir + jpra];
            for (int l = 0; l < num_multi; l++) {
                for (int m = 0; m < shiftnum; m++) {
                    const int bit_index = j - l * shiftnum - m;
                    if (bit_index < 0) continue;
                    for (int k = 0; k < bkP::targetP::n + 1; k++) {
                        temp1[i + bit_index][k] -= temps[l][k] << m;
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
    const SecretKey &)
{
    static_assert(numdigit > 1, "CLPX2TLWESIKSAnyBit needs a data digit");
    static_assert(basebit > 0, "CLPX2TLWESIKSAnyBit needs a positive basebit");
    static_assert(basebit == 2 || basebit == 4,
                  "the reverse switch supports basebit 2 or 4");
    constexpr bool single_carry = basebit >= 4;
    constexpr uint32_t carry_digits = single_carry ? 1 : 2;
    constexpr uint32_t carry_width = carry_digits * basebit;
    // HomDecomp layout: one noise digit, numdigit-1 data digits, followed by
    // either the legacy two base-4 carry digits or one base-16 carry digit.
    constexpr uint32_t block_size = (numdigit - 1) * basebit;
    const uint32_t tlnum = res.size();
    const uint32_t epoch_num = (tlnum + block_size - 1) / block_size;
    TLWE<typename bkP02::targetP> sumpra = {};
    TLWE<typename bkP02::targetP> previous_scaled = {};
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

    for (int i = 0; i < bkP02::targetP::n; i++)
        testVectorb[i] = static_cast<typename bkP02::targetP::T>(
            static_cast<double>(i)
            * static_cast<double>(typename bkP02::targetP::T(1)
                                  << (std::numeric_limits<
                                          typename bkP02::targetP::T>::digits -
                                      (bkP02::targetP::nbit + 4))));

    for (uint32_t epoch = 0; epoch < epoch_num; epoch++) {
        const uint32_t block_offset = epoch * block_size;
        const uint32_t batch_num =
            std::min(block_size, tlnum - block_offset);

        for (uint32_t i = 0; i < batch_num; i++) {
            TLWE<typename iksP20::targetP> temp10 = {};
            TLWE<typename bkP02::targetP> fid_stage1 = {};
            TLWE<typename bkP02::targetP> fid_stage2 = {};
            TLWE<typename bkP02::targetP> x_minus_b = {};
            TLWE<typename bkP02::domainP> temp = {};

            SampleExtractIndex<typename iksP20::domainP>(
                fid_stage1, trlwe, block_offset + i);
            IdentityKeySwitch<iksP20>(temp, fid_stage1,
                                      ek.getiksk<iksP20>());
            GateBootstrappingTLWE2TLWEFFT<bkP02>(
                fid_stage1, temp, ek.getbkfft<bkP02>(), testVectora);

            IdentityKeySwitch<iksP20>(temp, fid_stage1,
                                      ek.getiksk<iksP20>());
            GateBootstrappingTLWE2TLWEFFT<bkP02>(
                fid_stage2, temp, ek.getbkfft<bkP02>(), testVectorb);

            for (int j = 0; j < bkP02::targetP::n + 1; j++)
                x_minus_b[j] = previous_scaled[j] - (fid_stage2[j] << 1);
            IdentityKeySwitch<iksP20>(temp10, x_minus_b,
                                      ek.getiksk<iksP20>());
            previous_scaled = fid_stage2;

            // Algorithm Bign2TFHEanyB: round the centered (x-b) digit and
            // apply its block weight in one PBS.  For CLPX radix b=2 the
            // input is split into 2^(2+log2(b))=8 centered intervals.
            constexpr uint32_t num_intervals = 8;
            static_assert(bkP02::targetP::n % num_intervals == 0);
            constexpr uint32_t interval_size =
                bkP02::targetP::n / num_intervals;
            Polynomial<typename bkP02::targetP> rounded_digit_tv = {};
            const auto unit = typename bkP02::targetP::T(1)
                              << (std::numeric_limits<
                                  typename bkP02::targetP::T>::digits -
                                  batch_num + i - carry_width);
            for (uint32_t interval = 0; interval < num_intervals;
                 interval++) {
                // Blind rotation addresses 2N torus positions, so an N/8
                // table interval represents 1/16 of the torus.  The +1/32
                // input bias centers digit d in its interval; negative phases
                // address the upper intervals with an implicit sign flip.
                const uint32_t magnitude =
                    interval <= num_intervals / 2
                        ? interval
                        : num_intervals - interval;
                const auto value =
                    static_cast<typename bkP02::targetP::T>(magnitude) * unit;
                for (uint32_t j = 0; j < interval_size; j++)
                    rounded_digit_tv[interval * interval_size + j] = value;
            }

            temp10[iksP20::targetP::n] +=
                typename iksP20::targetP::T(1)
                << (std::numeric_limits<typename iksP20::targetP::T>::digits -
                    5);
            TLWE<typename bkP02::targetP> rounded_digit = {};
            GateBootstrappingTLWE2TLWEFFT<bkP02>(
                rounded_digit, temp10, ek.getbkfft<bkP02>(),
                rounded_digit_tv);
            for (int j = 0; j < bkP02::targetP::n + 1; j++)
                sumpra[j] += rounded_digit[j];
        }

        sumpra[bkP02::targetP::n] += static_cast<typename bkP02::targetP::T>(
            typename bkP02::targetP::T(1)
            << (std::numeric_limits<typename bkP02::targetP::T>::digits -
                batch_num - 6));
        std::array<TLWE<typename bkP01::targetP>, numdigit + carry_digits>
            sums;
        HomDecomp<iksP21, iksP10, bkP01, basebit,
                  numdigit + carry_digits,
                  basebit * (numdigit + carry_digits)>(
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

            // Center the base-2^basebit digit in its PBS decision interval.
            sums[i][iksP10::domainP::n] += static_cast<typename iksP10::domainP::T>(
                typename iksP10::domainP::T(1)
                << (std::numeric_limits<typename iksP10::domainP::T>::digits -
                    basebit - 1));
            IdentityKeySwitch<iksP10>(temp1, sums[i], ek.getiksk<iksP10>());

            for (int k = 0; k < static_cast<int>(basebit) - 1; k++) {
                for (int j = 0; j < iksP10::targetP::n + 1; j++)
                    temp2[j] = temp1[j] << (basebit - k - 1);

                const uint32_t output_index =
                    block_offset + (i - 1) * basebit + k;
                if (output_index < tlnum)
                    GateBootstrappingTLWE2TLWEFFT<bkP01>(
                        res[output_index], temp2, ek.getbkfft<bkP01>(),
                        testVector3);
            }

            const uint32_t output_index = block_offset + i * basebit - 1;
            if (output_index < tlnum)
                GateBootstrappingTLWE2TLWEFFT<bkP01>(
                    res[output_index], temp1, ek.getbkfft<bkP01>(),
                    testVector3);
        }

        if (epoch < epoch_num - 1) {
            if constexpr (single_carry) {
                sums[numdigit][iksP10::domainP::n] +=
                    typename iksP10::domainP::T(1)
                    << (std::numeric_limits<
                            typename iksP10::domainP::T>::digits -
                        basebit - 1);
                TLWE<typename iksP20::targetP> carry_input = {};
                IdentityKeySwitch<iksP10>(carry_input, sums[numdigit],
                                          ek.getiksk<iksP10>());

                constexpr uint32_t carry_intervals = 8;
                static_assert(bkP02::targetP::n % carry_intervals == 0);
                constexpr uint32_t carry_interval_size =
                    bkP02::targetP::n / carry_intervals;
                const auto carry_unit = typename bkP02::targetP::T(1)
                                        << (std::numeric_limits<
                                                typename bkP02::targetP::T>::digits -
                                            batch_num - basebit);
                Polynomial<typename bkP02::targetP> carry_tv = {};
                // The triangular table recovers the centered base-16 carry;
                // carry_unit applies the attenuation required by the next
                // block's accumulator.
                for (uint32_t interval = 0; interval < carry_intervals;
                     interval++) {
                    const uint32_t magnitude =
                        interval <= carry_intervals / 2
                            ? interval
                            : carry_intervals - interval;
                    const auto value =
                        static_cast<typename bkP02::targetP::T>(magnitude) *
                        carry_unit;
                    for (uint32_t j = 0; j < carry_interval_size; j++)
                        carry_tv[interval * carry_interval_size + j] = value;
                }
                GateBootstrappingTLWE2TLWEFFT<bkP02>(
                    sumpra, carry_input, ek.getbkfft<bkP02>(), carry_tv);
            }
            else {
                std::array<Polynomial<typename bkP02::targetP>, 3> carry_tvs =
                    {};
                TLWE<typename iksP20::targetP> shifted_carry = {};
                std::array<TLWE<typename iksP20::targetP>, 2> carry_inputs{};
                TLWE<typename bkP02::targetP> carry_high = {};
                TLWE<typename bkP02::targetP> carry_low = {};
                TLWE<typename bkP02::targetP> carry_shifted = {};

                for (int i = 0; i < 2; i++) {
                    sums[numdigit + i][iksP10::domainP::n] +=
                        typename iksP10::domainP::T(1)
                        << (std::numeric_limits<
                                typename iksP10::domainP::T>::digits -
                            3);
                    IdentityKeySwitch<iksP10>(carry_inputs[i],
                                              sums[numdigit + i],
                                              ek.getiksk<iksP10>());
                }

                for (int k = 0; k < 3; k++)
                    for (int j = 0; j < bkP02::targetP::n; j++)
                        carry_tvs[k][j] =
                            typename bkP02::targetP::T(1)
                            << (std::numeric_limits<
                                    typename bkP02::targetP::T>::digits -
                                batch_num - k - 3);

                GateBootstrappingTLWE2TLWEFFT<bkP02>(
                    carry_high, carry_inputs[1], ek.getbkfft<bkP02>(),
                    carry_tvs[0]);
                carry_high[bkP02::targetP::n] -=
                    typename bkP02::targetP::T(1)
                    << (std::numeric_limits<
                            typename bkP02::targetP::T>::digits -
                        batch_num - 5);
                GateBootstrappingTLWE2TLWEFFT<bkP02>(
                    carry_low, carry_inputs[0], ek.getbkfft<bkP02>(),
                    carry_tvs[1]);
                for (int j = 0; j < iksP20::targetP::n + 1; j++)
                    shifted_carry[j] = carry_inputs[0][j] << 1;
                GateBootstrappingTLWE2TLWEFFT<bkP02>(
                    carry_shifted, shifted_carry, ek.getbkfft<bkP02>(),
                    carry_tvs[2]);
                for (int j = 0; j < bkP02::targetP::n + 1; j++)
                    sumpra[j] =
                        carry_high[j] - carry_low[j] - carry_shifted[j];
            }
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
    std::array<double, P::n> delb{};
    std::array<double, P::n> delbM{};
    std::array<typename P::T, P::n> res{};
    const double q = std::ldexp(1.0, std::numeric_limits<typename P::T>::digits - 1);

    for (int i = 0; i < P::n; i++) delb[i] = -q / std::pow(P::plain_modulus, i + 1);
    for (int i = 0; i < P::n; i++) {
        for (int j = 0; j < P::n; j++) {
            if (i + j < P::n)
                delbM[i + j] += delb[i] * p[j];
            else
                delbM[i + j - P::n] -= delb[i] * p[j];
        }
    }
    for (int i = 0; i < P::n; i++)
        res[i] = clpxSignedDoubleToTorus<typename P::T>(delbM[i]);
    return res;
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
    return static_cast<typename P::T>(
        decodeHatEncoderInt2I128<P, validbit, numcell>(p));
}

template <class P, int validbit, int numcell>
inline unsigned __int128 decodeHatEncoderInt2U128(
    const std::array<int, numcell> &p)
{
    return static_cast<unsigned __int128>(
        decodeHatEncoderInt2I128<P, validbit, numcell>(p));
}

template <class P, int validbit, int numcell>
inline __int128 decodeHatEncoderInt2I128(const std::array<int, numcell> &p)
{
    __int128 ans = 0;
    for (int i = validbit - 1; i >= 0; i--) {
        ans *= P::plain_modulus;
        ans += p[i];
    }
    return ans;
}

}  // namespace TFHEpp
