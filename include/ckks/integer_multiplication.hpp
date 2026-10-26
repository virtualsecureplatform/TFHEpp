#pragma once

#include "scheme_switching.hpp"

namespace TFHEpp {

template <class L = lvl2param>
inline void TFHEEncryptInteger4(TLWE<L> &out, unsigned value, const Key<L> &key)
{
    static_assert(std::is_same_v<typename L::T, std::uint64_t>);
    if (value > 15) throw std::out_of_range("integer multiplier input must be in [0,15]");
    tlweSymEncrypt<L>(out, std::uint64_t{value} << 56, key);
}

template <class L = lvl2param>
inline std::uint8_t TFHEDecryptInteger8(const TLWE<L> &in, const Key<L> &key)
{
    static_assert(std::is_same_v<typename L::T, std::uint64_t>);
    const auto phase = tlweSymPhase<L>(in, key);
    return static_cast<std::uint8_t>((phase + (std::uint64_t{1} << 55)) >> 56);
}

// Width-calibrated integer encoding; Bits=8 gives exact 16-bit products.
template <unsigned Bits, class L = lvl2param>
inline void TFHEEncryptInteger(TLWE<L> &out, unsigned value, const Key<L> &key)
{
    static_assert(Bits >= 1 && Bits <= 8);
    static_assert(std::is_same_v<typename L::T, std::uint64_t>);
    if (value >= (1U << Bits)) throw std::out_of_range("integer input exceeds configured width");
    tlweSymEncrypt<L>(out, std::uint64_t{value} << (64-2*Bits), key);
}

template <unsigned Bits, class L = lvl2param>
inline unsigned TFHEDecryptIntegerProduct(const TLWE<L> &in, const Key<L> &key)
{
    static_assert(Bits >= 1 && Bits <= 8);
    static_assert(std::is_same_v<typename L::T, std::uint64_t>);
    const auto phase = tlweSymPhase<L>(in, key);
    return static_cast<unsigned>((phase + (std::uint64_t{1} << (63-2*Bits))) >> (64-2*Bits));
}

// One integer per ciphertext (constant CKKS plaintext polynomial). This is
// a separate integer-only path; it does not produce canonical real lifts.
struct CKKSInteger4Param {
    static constexpr std::int32_t key_value_min = -1, key_value_max = 1;
    static constexpr std::uint32_t nbit = 14, n = 1U << nbit, k = 1;
    using T = MultiLimbUInt<4>;
    static constexpr std::uint32_t l = 4, lₐ = 4, Bgbit = 8, Bgₐbit = 8;
    static constexpr T Bg = T{1} << Bgbit, Bgₐ = Bg;
    static constexpr T μ = T{1} << 253;
    static constexpr std::uint32_t l̅ = 32, l̅ₐ = 32;
    static constexpr std::uint32_t B̅gbit = 8, B̅gₐbit = 8;
    static constexpr ErrorDistribution errordist = ErrorDistribution::ModularGaussian;
    static inline const double α = std::ldexp(3.2, -256);
    static constexpr double ckks_min_noise_stddev = 3.2;
};

namespace ckks_integer_detail {

// Exact product modulo 2^LogQ. Each NTT inverse combines at most chunks
// convolutions of 16-bit unsigned digits. The static bound prevents modular
// aliasing of their signed integer coefficients. No floating FFT is used.
template <class P, std::uint32_t LogQ>
inline void exactProduct(Polynomial<P> &out, const Polynomial<P> &a,
                         const Polynomial<P> &b)
{
    static_assert(is_multilimb_uint_v<typename P::T> && P::nbit <= 15);
    static_assert(LogQ <= std::numeric_limits<typename P::T>::digits);
    constexpr std::uint32_t chunks = (LogQ + 15) / 16;
    constexpr std::uint64_t prime = 4611686018427322369ULL;
    static_assert(static_cast<unsigned __int128>(chunks) * P::n * 65535 * 65535 < prime / 2);
    static const modular_ntt::NegacyclicNTTPlan plan(P::n, {prime, 7});
    std::vector<std::vector<std::uint64_t>> af(chunks, std::vector<std::uint64_t>(P::n));
    auto bf = af;
    for (std::uint32_t j = 0; j < chunks; j++) {
        for (std::uint32_t i = 0; i < P::n; i++) {
            af[j][i] = static_cast<std::uint64_t>((a[i] >> (16*j)) & typename P::T{65535});
            bf[j][i] = static_cast<std::uint64_t>((b[i] >> (16*j)) & typename P::T{65535});
        }
        plan.forward(af[j]);
        plan.forward(bf[j]);
    }
    out.fill(0);
    std::vector<std::uint64_t> acc(P::n);
    for (std::uint32_t weight = 0; weight < chunks; weight++) {
        std::fill(acc.begin(), acc.end(), 0);
        for (std::uint32_t j = 0; j <= weight; j++)
            modular_ntt::multiplyAddVector(acc.data(), af[j].data(), bf[weight-j].data(), P::n, prime);
        plan.inverse(acc);
        for (std::uint32_t i = 0; i < P::n; i++) {
            const auto v = acc[i] <= prime/2 ? static_cast<std::int64_t>(acc[i])
                                            : -static_cast<std::int64_t>(prime-acc[i]);
            out[i] += P::T::from_signed_i64(v) << (16*weight);
        }
    }
    for (auto &v : out) v = ckks_detail::reduceToLevel<P, LogQ>(v);
}

}  // namespace ckks_integer_detail

template <class P = CKKSInteger4Param, class L = lvl2param>
struct TFHECKKSInteger4Key {
    static_assert(P::k == 1 && L::k == 1);
    static_assert(std::is_same_v<typename L::T, std::uint64_t>);
    static_assert(std::numeric_limits<typename P::T>::digits == 256);
    // Each row encrypts s_TFHE[j] as a constant under the full CKKS key.
    // Store full masks, independently of seeded-mask archive formats.
    std::vector<CKKSCiphertext<P, 256, 64>> to_ckks;
    CKKSRelinKey<P, 128> relin;
    std::unique_ptr<CKKSToTFHEKeySwitchKey<P, L, 2, 24>> to_tfhe;
};

template <class P, class L>
inline void TFHECKKSInteger4KeyGen(TFHECKKSInteger4Key<P, L> &out,
                                  const Key<P> &ckks_key, const Key<L> &tfhe_key)
{
    const CKKSNoise noise{std::ldexp(3.2, -256), 0};
    for (const auto &v : ckks_key)
        (void)ckks_scheme_switch_detail::ternarySecretToInt<P>(v);
    std::vector<std::int64_t> source_coefficients(L::n);
    for (std::uint32_t j = 0; j < L::n; j++)
        source_coefficients[j] = ckks_scheme_switch_detail::ternarySecretToInt<L>(tfhe_key[j]);
    out.to_ckks.resize(L::n);
#pragma omp parallel for
    for (std::uint32_t j = 0; j < L::n; j++) {
        auto message = std::make_unique<Polynomial<P>>();
        message->fill(0);
        (*message)[0] = P::T::from_signed_i64(source_coefficients[j]) << 64;
        ckks_detail::encryptPolynomialAtLevel<P, 256>(
            out.to_ckks[j].ct, *message, ckks_key, noise);
    }
    auto square = std::make_unique<Polynomial<P>>();
    ckks_detail::exactActiveTorusMulBySmall<P, 128>(*square, ckks_key, ckks_key);
    for (std::uint32_t row = 0; row < out.relin.key_switch_rows; row++) {
        auto message = std::make_unique<Polynomial<P>>();
        const auto shift = 128 - (row + 1) * P::B̅gbit;
        for (std::uint32_t i = 0; i < P::n; i++)
            (*message)[i] = ckks_detail::reduceToLevel<P, 128>((*square)[i] << shift);
        ckks_detail::encryptPolynomialAtLevel<P, 128>(
            out.relin.data[row], *message, ckks_key, noise);
    }
    out.to_tfhe = std::make_unique<CKKSToTFHEKeySwitchKey<P, L, 2, 24>>();
    CKKSToTFHEKeyGen<P, L, 2, 24>(*out.to_tfhe, ckks_key, tfhe_key);
}

// Exact public scalar combinations followed by rounded 64-bit rescaling.
// The ideal plaintext is the constant u=b/2^64-sum(a_j*s_j)/2^64, which
// equals the encoded integer fraction plus an integer and TFHE encryption error.
template <class P, class L>
inline void TFHEToCKKSInteger4Phase(CKKSCiphertext<P, 192, 64> &out,
                                   const TLWE<L> &in,
                                   const TFHECKKSInteger4Key<P, L> &key)
{
    if (key.to_ckks.size() != L::n)
        throw std::invalid_argument("integer switching key is not initialized");
    auto raw = std::make_unique<TRLWE<P>>();
    for (auto &p : *raw) p.fill(0);
    (*raw)[1][0] = typename P::T{in[L::n]} << 64;
    for (std::uint32_t j = 0; j < L::n; j++) {
        const typename P::T scalar{in[j]};
        for (std::uint32_t c = 0; c < 2; c++)
            for (std::uint32_t i = 0; i < P::n; i++)
                (*raw)[c][i] -= key.to_ckks[j].ct[c][i] * scalar;
    }
    for (std::uint32_t c = 0; c < 2; c++)
        for (std::uint32_t i = 0; i < P::n; i++)
            out.ct[c][i] = ckks_detail::roundedLevelRightShift<P, 256, 64>((*raw)[c][i]);
}

template <class P>
inline void CKKSInteger4Multiply(CKKSCiphertext<P, 128, 64> &out,
                                 const CKKSCiphertext<P, 192, 64> &lhs,
                                 const CKKSCiphertext<P, 192, 64> &rhs,
                                 const CKKSRelinKey<P, 128> &key)
{
    auto tensor = std::make_unique<TRLWE3<P>>();
    auto temporary = std::make_unique<Polynomial<P>>();
    ckks_integer_detail::exactProduct<P, 192>((*tensor)[2], lhs.ct[0], rhs.ct[0]);
    ckks_integer_detail::exactProduct<P, 192>((*tensor)[1], lhs.ct[1], rhs.ct[1]);
    ckks_integer_detail::exactProduct<P, 192>((*tensor)[0], lhs.ct[0], rhs.ct[1]);
    ckks_integer_detail::exactProduct<P, 192>(*temporary, lhs.ct[1], rhs.ct[0]);
    for (std::uint32_t i = 0; i < P::n; i++) (*tensor)[0][i] += (*temporary)[i];
    for (auto &p : *tensor)
        for (auto &v : p) v = ckks_detail::roundedLevelRightShift<P, 192, 64>(v);
    auto switched = std::make_unique<TRLWE<P>>();
    for (auto &p : *switched) p.fill(0);
    constexpr auto first = CKKSKeySwitchFirstRowForLevel<P, 128>();
    constexpr auto rows = CKKSKeySwitchRowCountForLevel<P, 128>();
    auto centered = std::make_unique<Polynomial<P>>();
    auto digits = std::make_unique<std::array<Polynomial<P>, rows>>();
    ckks_detail::centeredPolynomialAtLevel<P, 128>(*centered, (*tensor)[2]);
    ckks_detail::baseBbarDecomposePolynomialRows<P, first>(*digits, *centered);
    for (std::uint32_t row = 0; row < rows; row++)
        for (std::uint32_t c = 0; c < 2; c++) {
            ckks_detail::exactActiveTorusMulBySmall<P, 128>(
                *temporary, key.data[row][c], (*digits)[row]);
            for (std::uint32_t i = 0; i < P::n; i++)
                (*switched)[c][i] += (*temporary)[i];
        }
    for (std::uint32_t c = 0; c < 2; c++)
        for (std::uint32_t i = 0; i < P::n; i++)
            out.ct[c][i] = ckks_detail::reduceToLevel<P, 128>(
                (*tensor)[c][i] + (*switched)[c][i]);
}

// Input phases encode x/M and y/M, M=2^(2*Bits). The output encodes the
// exact 2*Bits-bit product /M. No secret key is used during evaluation.
// M*(x/M+k)*(y/M+l) == x*y/M (mod 1).
template <unsigned Bits, class P, class L>
inline void TFHECKKSIntegerMultiply(TLWE<L> &out, const TLWE<L> &lhs,
                                    const TLWE<L> &rhs,
                                    const TFHECKKSInteger4Key<P, L> &key)
{
    static_assert(Bits >= 1 && Bits <= 8);
    if (!key.to_tfhe) throw std::invalid_argument("missing output switching key");
    auto a = std::make_unique<CKKSCiphertext<P, 192, 64>>();
    auto b = std::make_unique<CKKSCiphertext<P, 192, 64>>();
    TFHEToCKKSInteger4Phase<P, L>(*a, lhs, key);
    TFHEToCKKSInteger4Phase<P, L>(*b, rhs, key);
    auto product = std::make_unique<CKKSCiphertext<P, 128, 64>>();
    CKKSInteger4Multiply<P>(*product, *a, *b, key.relin);
    // Extraction at q=48, scale=40 divides the represented value by 256.
    // Multiplying by 256*M gives output phase M*u*v, M=2^(2*Bits).
    CKKSMulIntegerInPlace<P, 128, 64>(*product, std::uint64_t{1} << (8+2*Bits));
    auto down = std::make_unique<CKKSCiphertext<P, 104, 40>>();
    CKKSScaleDown<P, 128, 64, 24>(*down, *product);
    auto reduced = std::make_unique<CKKSCiphertext<P, 48, 40>>();
    CKKSLevelReduce<P, 104, 48, 40>(*reduced, *down);
    CKKSCoeffToTFHE<P, L, 48, 40, 2, 24>(out, *reduced, 0, *key.to_tfhe);
}

template <class P, class L>
inline void TFHECKKSInteger4Multiply(TLWE<L> &out, const TLWE<L> &lhs,
                                    const TLWE<L> &rhs,
                                    const TFHECKKSInteger4Key<P, L> &key)
{
    TFHECKKSIntegerMultiply<4, P, L>(out, lhs, rhs, key);
}

}  // namespace TFHEpp
