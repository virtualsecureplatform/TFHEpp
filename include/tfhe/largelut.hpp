#pragma once

#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <type_traits>
#include <vector>

#include "gatebootstrapping.hpp"
#include "keyswitch.hpp"
#include "trlwe.hpp"

namespace TFHEpp {

namespace largelut_detail {

constexpr uint32_t floor_log2(uint32_t value)
{
    uint32_t res = 0;
    while (value >>= 1) res++;
    return res;
}

template <uint32_t W, uint32_t K>
constexpr uint32_t DigitCount()
{
    static_assert(W > 0);
    static_assert(K > 0);
    static_assert(K <= ((uint64_t{1} << W) - 1),
                  "LargeLUT requires K <= 2^W - 1");

    uint32_t remaining = K;
    uint32_t count = 0;
    for (uint32_t i = 0; remaining > 0; i++) {
        const uint32_t width = i == 0 ? W : W - 1 - floor_log2(i);
        if (width == 0) break;
        remaining -= width < remaining ? width : remaining;
        count++;
    }
    return count;
}

template <uint32_t W, uint32_t K>
constexpr auto DigitWidths()
{
    constexpr uint32_t digits = DigitCount<W, K>();
    std::array<uint32_t, digits> widths{};
    uint32_t remaining = K;
    for (uint32_t i = 0; i < digits; i++) {
        const uint32_t width = i == 0 ? W : W - 1 - floor_log2(i);
        widths[i] = width < remaining ? width : remaining;
        remaining -= widths[i];
    }
    return widths;
}

template <uint32_t W, uint32_t K>
constexpr auto DigitOffsets()
{
    constexpr uint32_t digits = DigitCount<W, K>();
    constexpr auto widths = DigitWidths<W, K>();
    std::array<uint32_t, digits> offsets{};
    uint32_t suffix = 0;
    for (uint32_t i = digits; i-- > 0;) {
        offsets[i] = suffix;
        suffix += widths[i];
    }
    return offsets;
}

template <class P, uint32_t MessageBits>
constexpr typename P::T EncodeMessage(const int64_t value)
{
    using T = typename P::T;
    static_assert(std::is_integral_v<T>);
    static_assert(MessageBits < std::numeric_limits<T>::digits,
                  "message modulus is too large for this torus type");

    constexpr uint32_t shift = std::numeric_limits<T>::digits - MessageBits;
    const uint64_t magnitude = value < 0 ? static_cast<uint64_t>(-value)
                                         : static_cast<uint64_t>(value);
    const T encoded = static_cast<T>(magnitude) << shift;
    return value < 0 ? T{} - encoded : encoded;
}

template <class P, uint32_t MessageBits>
void AddPlainMessage(TLWE<P> &tlwe, const int64_t value)
{
    tlwe[P::k * P::n] += EncodeMessage<P, MessageBits>(value);
}

inline int64_t CenteredMod(const int64_t value, const uint32_t modulus)
{
    assert(modulus > 0);
    int64_t res = value % static_cast<int64_t>(modulus);
    if (res < 0) res += modulus;
    if (res >= static_cast<int64_t>(modulus / 2)) res -= modulus;
    return res;
}

inline int64_t RoundedDiv(const int64_t value, const uint32_t divisor)
{
    assert(divisor > 0);
    if (value >= 0)
        return (value + static_cast<int64_t>(divisor / 2)) / divisor;
    return -((-value + static_cast<int64_t>(divisor / 2)) / divisor);
}

template <class P, uint32_t W>
int64_t CalibrationTableValue(const uint32_t x, const uint32_t amplitude_bits,
                              const uint32_t precision_bits)
{
    static_assert(W < P::nbit, "LargeLUT word size must be smaller than log N");
    constexpr uint32_t B = P::n >> W;
    constexpr uint32_t logB = P::nbit - W;
    assert(amplitude_bits < W);
    assert(precision_bits + amplitude_bits <= logB);

    const uint32_t period = B >> precision_bits;
    const uint32_t unit = B >> (precision_bits + amplitude_bits);
    const int64_t centered = CenteredMod(
        static_cast<int64_t>(x) - static_cast<int64_t>(B / 2), period);
    return RoundedDiv(centered, unit);
}

template <class P, uint32_t W>
Polynomial<P> CalibrationTable(const uint32_t amplitude_bits,
                               const uint32_t precision_bits)
{
    Polynomial<P> table{};
    for (uint32_t x = 0; x < P::n; x++)
        table[x] = EncodeMessage<P, W + 1>(
            CalibrationTableValue<P, W>(x, amplitude_bits, precision_bits));
    return table;
}

template <class P, uint32_t W, uint32_t K>
Polynomial<P> BatchCalibrationTable(const uint32_t digit_index)
{
    static_assert(W < P::nbit, "LargeLUT word size must be smaller than log N");
    constexpr uint32_t L = DigitCount<W, K>();
    constexpr uint32_t B = P::n >> W;
    constexpr auto widths = DigitWidths<W, K>();
    constexpr auto offsets = DigitOffsets<W, K>();
    assert(digit_index + 1 < L);

    Polynomial<P> table{};
    for (uint32_t x = 0; x < P::n; x++) {
        const uint32_t block = x / B;
        if (block >= L - digit_index - 1) continue;

        const uint32_t j = digit_index + 1 + block;
        const uint32_t precision = offsets[digit_index] - offsets[j - 1];
        table[x] = EncodeMessage<P, W + 1>(
            CalibrationTableValue<P, W>(x, widths[j], precision));
    }
    return table;
}

template <class P>
void TrivialTRLWEFromPolynomial(TRLWE<P> &res, const Polynomial<P> &poly)
{
    res = {};
    res[P::k] = poly;
}

template <class P>
void RotateTRLWEByMonomial(TRLWE<P> &res, const TRLWE<P> &src,
                           const uint32_t exponent)
{
    for (uint32_t k = 0; k < P::k + 1; k++)
        PolynomialMulByXai<P>(res[k], src[k], exponent);
}

template <class T, uint32_t t, uint32_t basebit>
constexpr T R2ROffset()
{
    static_assert(std::is_integral_v<T>);
    static_assert(basebit > 0);
    static_assert(basebit < std::numeric_limits<T>::digits);
    static_assert(t * basebit <= std::numeric_limits<T>::digits);

    T offset = 0;
    for (uint32_t i = 1; i <= t; i++)
        offset += (T{1} << (basebit - 1)) *
                  (T{1} << (std::numeric_limits<T>::digits - i * basebit));
    return offset;
}

template <uint32_t W, uint32_t K, class P>
uint32_t StepExtractIndex(const uint32_t i, const uint32_t p)
{
    constexpr uint32_t B = P::n >> W;
    constexpr auto offsets = DigitOffsets<W, K>();

    const uint64_t Bprev = B >> offsets[i - 1];
    const uint64_t index = (static_cast<uint64_t>(2 * p + 1) * Bprev) / 2;
    assert(index < P::n);
    return static_cast<uint32_t>(index);
}

template <uint32_t W, uint32_t K, class P>
uint32_t FinalExtractIndex()
{
    constexpr uint32_t B = P::n >> W;
    return B / 2;
}

constexpr uint32_t BitCeil(uint32_t value)
{
    uint32_t res = 1;
    while (res < value) res <<= 1;
    return res;
}

template <class brP, class ahP, uint32_t W, uint32_t K, uint32_t I>
void LargeLUTAHPackingStep(
    TRLWE<typename brP::targetP> &current,
    const std::array<TLWE<typename brP::domainP>, DigitCount<W, K>()>
        &calibrated,
    const BootstrappingKeyFFT<brP> &bkfft, const AnnihilateKey<ahP> &ahk)
{
    using targetP = typename brP::targetP;
    constexpr uint32_t B = targetP::n >> W;
    constexpr auto offsets = DigitOffsets<W, K>();
    constexpr uint32_t Si = (I + 1) * (uint32_t{1} << offsets[I - 1]);
    constexpr uint32_t Qi = uint32_t{1} << (W + offsets[I]);
    constexpr uint32_t Bi = B >> offsets[I];
    static_assert(Qi >= Si);
    static_assert(targetP::n / Qi == Bi);

    TRLWE<targetP> rotated;
    BlindRotate<brP>(rotated, calibrated[I - 1], bkfft, current);

    auto packed = std::make_unique<std::array<TLWE<targetP>, Qi>>();
    for (uint32_t p = 0; p < Si; p++) {
        const uint32_t index = StepExtractIndex<W, K, targetP>(I, p);
        SampleExtractIndex<targetP>((*packed)[p], rotated, index);
    }
    TLWE2TablePacking<ahP, Qi>(current, *packed, ahk);
}

template <class brP, class ahP, uint32_t W, uint32_t K, uint32_t I = 1>
void LargeLUTAHPackingSteps(
    TRLWE<typename brP::targetP> &current,
    const std::array<TLWE<typename brP::domainP>, DigitCount<W, K>()>
        &calibrated,
    const BootstrappingKeyFFT<brP> &bkfft, const AnnihilateKey<ahP> &ahk)
{
    if constexpr (I < DigitCount<W, K>()) {
        LargeLUTAHPackingStep<brP, ahP, W, K, I>(current, calibrated, bkfft,
                                                 ahk);
        LargeLUTAHPackingSteps<brP, ahP, W, K, I + 1>(current, calibrated,
                                                      bkfft, ahk);
    }
}

}  // namespace largelut_detail

template <uint32_t W, uint32_t K>
constexpr uint32_t LargeLUTDigitCount()
{
    return largelut_detail::DigitCount<W, K>();
}

template <uint32_t W, uint32_t K>
constexpr auto LargeLUTDigitWidths()
{
    return largelut_detail::DigitWidths<W, K>();
}

template <uint32_t W, uint32_t K>
constexpr auto LargeLUTDigitOffsets()
{
    return largelut_detail::DigitOffsets<W, K>();
}

template <uint32_t W, uint32_t K>
constexpr auto LargeLUTRadixDecompose(const uint32_t value)
{
    constexpr uint32_t digits = LargeLUTDigitCount<W, K>();
    constexpr auto widths = LargeLUTDigitWidths<W, K>();
    constexpr auto offsets = LargeLUTDigitOffsets<W, K>();
    std::array<uint32_t, digits> res{};
    for (uint32_t i = 0; i < digits; i++)
        res[i] = (value >> offsets[i]) & ((uint32_t{1} << widths[i]) - 1);
    return res;
}

template <class P, uint32_t W>
std::vector<uint32_t> LargeLUTBatchR2RPositions()
{
    static_assert(W < P::nbit, "LargeLUT word size must be smaller than log N");
    constexpr uint32_t B = P::n >> W;
    std::vector<uint32_t> positions(uint32_t{1} << W);
    for (uint32_t j = 0; j < positions.size(); j++)
        positions[j] = ((2 * j + 1) * B) / 2;
    return positions;
}

template <class P, uint32_t W, uint32_t K>
std::vector<uint32_t> LargeLUTStepR2RPositions(const uint32_t i)
{
    static_assert(W < P::nbit, "LargeLUT word size must be smaller than log N");
    constexpr uint32_t L = LargeLUTDigitCount<W, K>();
    constexpr uint32_t B = P::n >> W;
    constexpr auto offsets = LargeLUTDigitOffsets<W, K>();
    assert(i > 0);
    assert(i < L);

    const uint32_t Bprev = B >> offsets[i - 1];
    const uint32_t Si = (i + 1) * (uint32_t{1} << offsets[i - 1]);
    std::vector<uint32_t> positions(Si);
    for (uint32_t j = 0; j < Si; j++)
        positions[j] = ((2 * j + 1) * Bprev) / 2;
    return positions;
}

template <class P, uint32_t W, uint32_t K>
uint32_t LargeLUTStepBlockSize(const uint32_t i)
{
    static_assert(W < P::nbit, "LargeLUT word size must be smaller than log N");
    constexpr uint32_t L = LargeLUTDigitCount<W, K>();
    constexpr uint32_t B = P::n >> W;
    constexpr auto offsets = LargeLUTDigitOffsets<W, K>();
    assert(i > 0);
    assert(i < L);
    return B >> offsets[i];
}

template <class P, uint32_t W, uint32_t K, class Table>
Polynomial<P> LargeLUTPolynomial(const Table &table)
{
    static_assert(W < P::nbit, "LargeLUT word size must be smaller than log N");
    static_assert(K <= P::nbit, "LargeLUT table size must be at most N");

    constexpr uint32_t table_size = uint32_t{1} << K;
    constexpr uint32_t block_size = P::n >> K;
    Polynomial<P> poly{};
    for (uint32_t i = 0; i < table_size; i++) {
        const auto encoded = largelut_detail::EncodeMessage<P, W + 1>(
            static_cast<int64_t>(table[i]));
        for (uint32_t j = 0; j < block_size; j++)
            poly[i * block_size + j] = encoded;
    }
    return poly;
}

template <class P>
void TrivialTRLWEFromPolynomial(TRLWE<P> &res, const Polynomial<P> &poly)
{
    largelut_detail::TrivialTRLWEFromPolynomial<P>(res, poly);
}

template <class P, uint32_t t, uint32_t basebit>
using R2RKey =
    std::vector<
        std::array<std::array<TRLWE<P>, (std::size_t{1} << (basebit - 1))>,
                   t>>;

template <class P, uint32_t t, uint32_t basebit>
void R2RKeyGen(R2RKey<P, t, basebit> &r2rk,
               const std::vector<uint32_t> &positions, const uint32_t R,
               const Key<P> &key)
{
    using T = typename P::T;
    static_assert(P::k == 1, "R2RPKS currently supports k=1 TRLWE keys");
    static_assert(std::is_integral_v<T>, "R2RPKS requires an integral torus");
    static_assert(basebit > 0);
    static_assert(basebit < std::numeric_limits<T>::digits);
    static_assert(t * basebit <= std::numeric_limits<T>::digits);

    assert(!positions.empty());
    assert(R > 0);
    assert(positions.size() * static_cast<std::size_t>(R) <= P::n);
    for (const uint32_t position : positions) assert(position < P::n);

    constexpr uint32_t halfbase = uint32_t{1} << (basebit - 1);
    r2rk.resize(P::n);
    for (uint32_t i = 0; i < P::n; i++) {
        for (uint32_t j = 0; j < t; j++) {
            constexpr uint32_t digits = std::numeric_limits<T>::digits;
            const uint32_t shift = digits - (j + 1) * basebit;
            for (uint32_t u = 0; u < halfbase; u++) {
                trlweSymEncryptZero<P>(r2rk[i][j][u], key);
                for (uint32_t p = 0; p < positions.size(); p++) {
                    const uint32_t position = positions[p];
                    const uint32_t gamma =
                        position >= i ? position - i : P::n + position - i;
                    T term = key[gamma] * static_cast<T>(u + 1);
                    term <<= shift;

                    for (uint32_t ell = 0; ell < R; ell++) {
                        T &coeff = r2rk[i][j][u][P::k][p * R + ell];
                        if (position >= i)
                            coeff += term;
                        else
                            coeff -= term;
                    }
                }
            }
        }
    }
}

template <class P, uint32_t t, uint32_t basebit>
void R2RKeyGen(R2RKey<P, t, basebit> &r2rk,
               const std::vector<uint32_t> &positions, const uint32_t R,
               const SecretKey &sk)
{
    R2RKeyGen<P, t, basebit>(r2rk, positions, R, sk.key.get<P>());
}

template <class P, uint32_t t, uint32_t basebit, uint32_t W, uint32_t K>
void LargeLUTR2RKeyGen(
    R2RKey<P, t, basebit> &batch_r2rk,
    std::array<R2RKey<P, t, basebit>, LargeLUTDigitCount<W, K>() - 1>
        &step_r2rks,
    const Key<P> &key)
{
    constexpr uint32_t L = LargeLUTDigitCount<W, K>();
    constexpr uint32_t B = P::n >> W;

    const auto batch_positions = LargeLUTBatchR2RPositions<P, W>();
    R2RKeyGen<P, t, basebit>(batch_r2rk, batch_positions, B, key);

    for (uint32_t i = 1; i < L; i++) {
        const auto positions = LargeLUTStepR2RPositions<P, W, K>(i);
        const uint32_t Bi = LargeLUTStepBlockSize<P, W, K>(i);
        R2RKeyGen<P, t, basebit>(step_r2rks[i - 1], positions, Bi, key);
    }
}

template <class P, uint32_t t, uint32_t basebit, uint32_t W, uint32_t K>
void LargeLUTR2RKeyGen(
    R2RKey<P, t, basebit> &batch_r2rk,
    std::array<R2RKey<P, t, basebit>, LargeLUTDigitCount<W, K>() - 1>
        &step_r2rks,
    const SecretKey &sk)
{
    LargeLUTR2RKeyGen<P, t, basebit, W, K>(batch_r2rk, step_r2rks,
                                           sk.key.get<P>());
}

template <class P, uint32_t t, uint32_t basebit>
void R2RPKS(TRLWE<P> &res, const TRLWE<P> &input,
            const std::vector<uint32_t> &positions, const uint32_t R,
            const R2RKey<P, t, basebit> &r2rk)
{
    using T = typename P::T;
    static_assert(P::k == 1, "R2RPKS currently supports k=1 TRLWE keys");
    static_assert(std::is_integral_v<T>, "R2RPKS requires an integral torus");
    static_assert(basebit > 0);
    static_assert(basebit < std::numeric_limits<T>::digits);
    static_assert(t * basebit <= std::numeric_limits<T>::digits);

    assert(!positions.empty());
    assert(R > 0);
    assert(positions.size() * static_cast<std::size_t>(R) <= P::n);
    assert(r2rk.size() == P::n);
    for (const uint32_t position : positions) assert(position < P::n);

    res = {};
    for (uint32_t p = 0; p < positions.size(); p++)
        for (uint32_t ell = 0; ell < R; ell++)
            res[P::k][p * R + ell] = input[P::k][positions[p]];

    constexpr uint32_t digits = std::numeric_limits<T>::digits;
    constexpr T roundoffset =
        (basebit * t) < digits ? T{1} << (digits - (1 + basebit * t)) : 0;
    constexpr T offset = largelut_detail::R2ROffset<T, t, basebit>();
    constexpr T mask = (T{1} << basebit) - 1;
    constexpr T halfbase = T{1} << (basebit - 1);

    for (uint32_t i = 0; i < P::n; i++) {
        const T aibar = input[0][i] + offset + roundoffset;
        for (uint32_t j = 0; j < t; j++) {
            const int32_t aij =
                static_cast<int32_t>((aibar >> (digits - (j + 1) * basebit)) &
                                     mask) -
                static_cast<int32_t>(halfbase);
            if (aij > 0)
                TRLWESub<P>(res, res, r2rk[i][j][aij - 1]);
            else if (aij < 0)
                TRLWEAdd<P>(res, res, r2rk[i][j][-aij - 1]);
        }
    }
}

template <class brP, class iksP, uint32_t W, uint32_t K>
void LargeLUTCalibrate(
    std::array<TLWE<typename brP::domainP>, LargeLUTDigitCount<W, K>()> &res,
    const std::array<TLWE<typename brP::domainP>, LargeLUTDigitCount<W, K>()>
        &digits,
    const BootstrappingKeyFFT<brP> &bkfft, const KeySwitchingKey<iksP> &iksk)
{
    using domainP = typename brP::domainP;
    using targetP = typename brP::targetP;
    static_assert(std::is_same_v<typename iksP::domainP, targetP>);
    static_assert(std::is_same_v<typename iksP::targetP, domainP>);
    static_assert(W < targetP::nbit);

    constexpr uint32_t L = LargeLUTDigitCount<W, K>();
    constexpr uint32_t B = targetP::n >> W;
    constexpr auto widths = LargeLUTDigitWidths<W, K>();
    constexpr auto offsets = LargeLUTDigitOffsets<W, K>();

    std::array<TLWE<domainP>, L> corrections{};
    res[0] = digits[0];
    for (uint32_t i = 0; i + 1 < L; i++) {
        TLWE<domainP> shifted = res[i];
        largelut_detail::AddPlainMessage<domainP, W + 1>(
            shifted, static_cast<int64_t>(i) * (int64_t{1} << (widths[i] - 1)));

        for (uint32_t j = i + 1; j < L; j++) {
            const uint32_t precision = offsets[i] - offsets[j - 1];
            const auto caltab = largelut_detail::CalibrationTable<targetP, W>(
                widths[j], precision);

            TRLWE<targetP> acc;
            BlindRotate<brP>(acc, shifted, bkfft, caltab);

            TLWE<targetP> extracted;
            SampleExtractIndex<targetP>(extracted, acc, B / 2);

            TLWE<domainP> switched;
            IdentityKeySwitch<iksP>(switched, extracted, iksk);
            for (uint32_t a = 0; a <= domainP::k * domainP::n; a++)
                corrections[j][a] += switched[a];
        }

        for (uint32_t a = 0; a <= domainP::k * domainP::n; a++)
            res[i + 1][a] = digits[i + 1][a] - corrections[i + 1][a];
    }
}

template <class brP, class iksP, uint32_t t, uint32_t basebit, uint32_t W,
          uint32_t K>
void LargeLUTBatchCalibrate(
    std::array<TLWE<typename brP::domainP>, LargeLUTDigitCount<W, K>()> &res,
    const std::array<TLWE<typename brP::domainP>, LargeLUTDigitCount<W, K>()>
        &digits,
    const BootstrappingKeyFFT<brP> &bkfft, const KeySwitchingKey<iksP> &iksk,
    const R2RKey<typename brP::targetP, t, basebit> &batch_r2rk)
{
    using domainP = typename brP::domainP;
    using targetP = typename brP::targetP;
    static_assert(std::is_same_v<typename iksP::domainP, targetP>);
    static_assert(std::is_same_v<typename iksP::targetP, domainP>);
    static_assert(W < targetP::nbit);

    constexpr uint32_t L = LargeLUTDigitCount<W, K>();
    constexpr uint32_t B = targetP::n >> W;
    constexpr auto widths = LargeLUTDigitWidths<W, K>();

    const auto batch_positions = LargeLUTBatchR2RPositions<targetP, W>();
    std::array<TLWE<domainP>, L> corrections{};
    res[0] = digits[0];
    for (uint32_t i = 0; i + 1 < L; i++) {
        TLWE<domainP> shifted = res[i];
        largelut_detail::AddPlainMessage<domainP, W + 1>(
            shifted, static_cast<int64_t>(i) * (int64_t{1} << (widths[i] - 1)));

        const auto batch_table =
            largelut_detail::BatchCalibrationTable<targetP, W, K>(i);
        TRLWE<targetP> first_acc;
        BlindRotate<brP>(first_acc, shifted, bkfft, batch_table);

        TRLWE<targetP> packed_table;
        R2RPKS<targetP, t, basebit>(packed_table, first_acc, batch_positions, B,
                                    batch_r2rk);

        TLWE<domainP> neg_shifted;
        for (uint32_t a = 0; a <= domainP::k * domainP::n; a++)
            neg_shifted[a] = typename domainP::T{} - shifted[a];

        TRLWE<targetP> second_acc;
        BlindRotate<brP>(second_acc, neg_shifted, bkfft, packed_table);

        for (uint32_t j = i + 1; j < L; j++) {
            const uint32_t block = j - (i + 1);
            const uint32_t index = ((2 * block + 1) * B) / 2;

            TLWE<targetP> extracted;
            SampleExtractIndex<targetP>(extracted, second_acc, index);

            TLWE<domainP> switched;
            IdentityKeySwitch<iksP>(switched, extracted, iksk);
            for (uint32_t a = 0; a <= domainP::k * domainP::n; a++)
                corrections[j][a] += switched[a];
        }

        for (uint32_t a = 0; a <= domainP::k * domainP::n; a++)
            res[i + 1][a] = digits[i + 1][a] - corrections[i + 1][a];
    }
}

template <class brP, class ahP, uint32_t W, uint32_t K>
void LargeLUTWithCalibrated(
    TLWE<typename brP::targetP> &res, const TRLWE<typename brP::targetP> &table,
    const std::array<TLWE<typename brP::domainP>, LargeLUTDigitCount<W, K>()>
        &calibrated,
    const BootstrappingKeyFFT<brP> &bkfft, const AnnihilateKey<ahP> &ahk)
{
    using targetP = typename brP::targetP;
    static_assert(std::is_same_v<TRLWE<ahP>, TRLWE<targetP>>);
    static_assert(std::is_same_v<TLWE<ahP>, TLWE<targetP>>);
    static_assert(W < targetP::nbit);
    static_assert(K <= targetP::nbit);

    constexpr uint32_t L = LargeLUTDigitCount<W, K>();

    TRLWE<targetP> current = table;
    largelut_detail::LargeLUTAHPackingSteps<brP, ahP, W, K>(
        current, calibrated, bkfft, ahk);

    TRLWE<targetP> final_acc;
    BlindRotate<brP>(final_acc, calibrated[L - 1], bkfft, current);
    SampleExtractIndex<targetP>(
        res, final_acc, largelut_detail::FinalExtractIndex<W, K, targetP>());
}

template <class brP, uint32_t t, uint32_t basebit, uint32_t W, uint32_t K>
void LargeLUTOptimizedWithCalibrated(
    TLWE<typename brP::targetP> &res, const TRLWE<typename brP::targetP> &table,
    const std::array<TLWE<typename brP::domainP>, LargeLUTDigitCount<W, K>()>
        &calibrated,
    const BootstrappingKeyFFT<brP> &bkfft,
    const std::array<R2RKey<typename brP::targetP, t, basebit>,
                     LargeLUTDigitCount<W, K>() - 1> &step_r2rks)
{
    using targetP = typename brP::targetP;
    static_assert(W < targetP::nbit);
    static_assert(K <= targetP::nbit);

    constexpr uint32_t L = LargeLUTDigitCount<W, K>();
    constexpr uint32_t B = targetP::n >> W;

    TRLWE<targetP> current = table;
    for (uint32_t i = 1; i < L; i++) {
        TRLWE<targetP> rotated;
        BlindRotate<brP>(rotated, calibrated[i - 1], bkfft, current);

        const auto positions = LargeLUTStepR2RPositions<targetP, W, K>(i);
        const uint32_t Bi = LargeLUTStepBlockSize<targetP, W, K>(i);
        R2RPKS<targetP, t, basebit>(current, rotated, positions, Bi,
                                    step_r2rks[i - 1]);
    }

    TRLWE<targetP> final_acc;
    BlindRotate<brP>(final_acc, calibrated[L - 1], bkfft, current);
    SampleExtractIndex<targetP>(res, final_acc, B / 2);
}

template <class brP, class iksP, class ahP, uint32_t W, uint32_t K>
void LargeLUT(TLWE<typename brP::targetP> &res,
              const TRLWE<typename brP::targetP> &table,
              const std::array<TLWE<typename brP::domainP>,
                               LargeLUTDigitCount<W, K>()> &digits,
              const BootstrappingKeyFFT<brP> &bkfft,
              const KeySwitchingKey<iksP> &iksk, const AnnihilateKey<ahP> &ahk)
{
    std::array<TLWE<typename brP::domainP>, LargeLUTDigitCount<W, K>()>
        calibrated;
    LargeLUTCalibrate<brP, iksP, W, K>(calibrated, digits, bkfft, iksk);
    LargeLUTWithCalibrated<brP, ahP, W, K>(res, table, calibrated, bkfft, ahk);
}

template <class brP, class iksP, class ahP, uint32_t W, uint32_t K>
void LargeLUT(TLWE<typename brP::targetP> &res,
              const Polynomial<typename brP::targetP> &table,
              const std::array<TLWE<typename brP::domainP>,
                               LargeLUTDigitCount<W, K>()> &digits,
              const BootstrappingKeyFFT<brP> &bkfft,
              const KeySwitchingKey<iksP> &iksk, const AnnihilateKey<ahP> &ahk)
{
    TRLWE<typename brP::targetP> trlwe;
    TrivialTRLWEFromPolynomial<typename brP::targetP>(trlwe, table);
    LargeLUT<brP, iksP, ahP, W, K>(res, trlwe, digits, bkfft, iksk, ahk);
}

template <class brP, class iksP, uint32_t t, uint32_t basebit, uint32_t W,
          uint32_t K>
void LargeLUTOptimized(
    TLWE<typename brP::targetP> &res, const TRLWE<typename brP::targetP> &table,
    const std::array<TLWE<typename brP::domainP>, LargeLUTDigitCount<W, K>()>
        &digits,
    const BootstrappingKeyFFT<brP> &bkfft, const KeySwitchingKey<iksP> &iksk,
    const R2RKey<typename brP::targetP, t, basebit> &batch_r2rk,
    const std::array<R2RKey<typename brP::targetP, t, basebit>,
                     LargeLUTDigitCount<W, K>() - 1> &step_r2rks)
{
    std::array<TLWE<typename brP::domainP>, LargeLUTDigitCount<W, K>()>
        calibrated;
    LargeLUTBatchCalibrate<brP, iksP, t, basebit, W, K>(
        calibrated, digits, bkfft, iksk, batch_r2rk);
    LargeLUTOptimizedWithCalibrated<brP, t, basebit, W, K>(
        res, table, calibrated, bkfft, step_r2rks);
}

template <class brP, class iksP, uint32_t t, uint32_t basebit, uint32_t W,
          uint32_t K>
void LargeLUTOptimized(
    TLWE<typename brP::targetP> &res,
    const Polynomial<typename brP::targetP> &table,
    const std::array<TLWE<typename brP::domainP>, LargeLUTDigitCount<W, K>()>
        &digits,
    const BootstrappingKeyFFT<brP> &bkfft, const KeySwitchingKey<iksP> &iksk,
    const R2RKey<typename brP::targetP, t, basebit> &batch_r2rk,
    const std::array<R2RKey<typename brP::targetP, t, basebit>,
                     LargeLUTDigitCount<W, K>() - 1> &step_r2rks)
{
    TRLWE<typename brP::targetP> trlwe;
    TrivialTRLWEFromPolynomial<typename brP::targetP>(trlwe, table);
    LargeLUTOptimized<brP, iksP, t, basebit, W, K>(
        res, trlwe, digits, bkfft, iksk, batch_r2rk, step_r2rks);
}

}  // namespace TFHEpp
