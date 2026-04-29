#pragma once

#include <array>
#include <cassert>
#include <cstdint>
#include <limits>
#include <memory>
#include <type_traits>
#include <vector>

#include "bfv-c2s.hpp"
#include "bfv-digitext.hpp"

namespace TFHEpp {
namespace bfvboot {

namespace detail {

constexpr uint64_t mulmod(uint64_t a, uint64_t b, uint64_t mod)
{
    return static_cast<uint64_t>(
        (static_cast<unsigned __int128>(a) * b) % mod);
}

constexpr uint64_t powmod(uint64_t base, uint64_t exp, uint64_t mod)
{
    uint64_t res = 1 % mod;
    base %= mod;
    while (exp != 0) {
        if ((exp & 1) != 0) res = mulmod(res, base, mod);
        base = mulmod(base, base, mod);
        exp >>= 1;
    }
    return res;
}

constexpr uint64_t modinv(uint64_t a, uint64_t mod)
{
    __int128 old_r = static_cast<__int128>(a % mod);
    __int128 r = static_cast<__int128>(mod);
    __int128 old_s = 1;
    __int128 s = 0;
    while (r != 0) {
        const __int128 q = old_r / r;
        const __int128 next_r = old_r - q * r;
        old_r = r;
        r = next_r;
        const __int128 next_s = old_s - q * s;
        old_s = s;
        s = next_s;
    }
    if (old_r < 0) old_s = -old_s;
    old_s %= static_cast<__int128>(mod);
    if (old_s < 0) old_s += mod;
    return static_cast<uint64_t>(old_s);
}

constexpr uint32_t ceil_log2(uint64_t x)
{
    uint32_t bits = 0;
    uint64_t v = x - 1;
    while (v != 0) {
        bits++;
        v >>= 1;
    }
    return bits;
}

constexpr uint64_t hensel_lift_root_to_square_modulus(uint64_t root_mod_p,
                                                       uint64_t order,
                                                       uint64_t p)
{
    const uint64_t p2 = p * p;
    const uint64_t root_pow = powmod(root_mod_p, order, p2);
    const uint64_t c = (root_pow + p2 - 1) / p;
    const uint64_t root_inv = modinv(root_mod_p, p);
    const uint64_t denom = mulmod(order % p, root_inv, p);
    const uint64_t k = mulmod((p - c % p) % p, modinv(denom, p), p);
    return root_mod_p + k * p;
}

template <class T>
constexpr double torus_delta_double(uint32_t plain_modulusbit)
{
    if constexpr (is_multilimb_uint_v<T>)
        return 0.0;
    else
        return static_cast<double>(
            static_cast<T>(1) << (std::numeric_limits<T>::digits - plain_modulusbit));
}

}  // namespace detail

// Parameter adapter for the p^2 plaintext ring used by the first BFV
// bootstrapping implementation. It keeps the ciphertext/security parameters of
// BaseP but changes the BFV plaintext modulus and lifts the SIMD root from p to
// p^2.
template <class BaseP>
struct PrimePower2Param {
    using T = typename BaseP::T;
    static_assert(std::is_same_v<T, __uint128_t> || is_multilimb_uint_v<T>,
                  "PrimePower2Param currently expects 128-bit or multi-limb BFV params");

    static constexpr int32_t key_value_max = BaseP::key_value_max;
    static constexpr int32_t key_value_min = BaseP::key_value_min;
    static constexpr std::uint32_t nbit = BaseP::nbit;
    static constexpr std::uint32_t n = BaseP::n;
    static constexpr std::uint32_t k = BaseP::k;
    static constexpr std::uint32_t l = BaseP::l;
    static constexpr std::uint32_t lₐ = BaseP::lₐ;
    static constexpr std::uint32_t Bgbit = BaseP::Bgbit;
    static constexpr std::uint32_t Bgₐbit = BaseP::Bgₐbit;
    static constexpr T Bg = BaseP::Bg;
    static constexpr T Bgₐ = BaseP::Bgₐ;
    static constexpr ErrorDistribution errordist = BaseP::errordist;
    static const inline double α = BaseP::α;
    static constexpr T μ = BaseP::μ;

    static constexpr uint64_t base_plain_modulus =
        static_cast<uint64_t>(BaseP::plain_modulus);
    static constexpr uint64_t bfv_bootstrap_digit_error_bound =
        BaseP::bfv_bootstrap_digit_error_bound;
    static_assert(bfv_bootstrap_digit_error_bound <=
                      (base_plain_modulus - 1) / 2,
                  "BFV bootstrap digit-error bound must fit in one p-digit");
    static constexpr int bfv_bootstrap_linear_bsgs_step =
        BaseP::bfv_bootstrap_linear_bsgs_step;
    static_assert(base_plain_modulus <=
                      std::numeric_limits<uint64_t>::max() / base_plain_modulus,
                  "p^2 must fit in uint64_t");
    static constexpr uint64_t plain_modulus_u64 =
        base_plain_modulus * base_plain_modulus;
    static constexpr T plain_modulus = static_cast<T>(plain_modulus_u64);
    static constexpr uint32_t plain_modulusbit =
        detail::ceil_log2(plain_modulus_u64);
    static constexpr double Δ = detail::torus_delta_double<T>(plain_modulusbit);
    static constexpr T delta_int =
        std::numeric_limits<T>::max() / plain_modulus_u64;
    static constexpr uint64_t Q_mod_t =
        (std::numeric_limits<T>::max() % plain_modulus_u64) + 1;

    static constexpr std::uint32_t l̅ = BaseP::l̅;
    static constexpr std::uint32_t l̅ₐ = BaseP::l̅ₐ;
    static constexpr std::uint32_t B̅gbit = BaseP::B̅gbit;
    static constexpr std::uint32_t B̅gₐbit = BaseP::B̅gₐbit;

    static constexpr uint64_t simd_modulus = plain_modulus_u64;
    static constexpr uint64_t simd_psi =
        detail::hensel_lift_root_to_square_modulus(
            SIMDConstants<BaseP>::PSI, 2 * BaseP::n, base_plain_modulus);
    static constexpr uint64_t simd_psi_inv =
        detail::modinv(simd_psi, simd_modulus);
    static constexpr uint64_t simd_n_inv =
        detail::modinv(BaseP::n % simd_modulus, simd_modulus);
};

template <class P>
void BfvPolyEncrypt(TRLWE<P> &ct, const std::array<uint64_t, P::n> &poly,
                    const Key<P> &key)
{
    Polynomial<P> scaled;
    for (uint32_t i = 0; i < P::n; i++)
        scaled[i] = bfvEncodeCoeff<P>(poly[i]);
    trlweSymEncrypt<P>(ct, scaled, key);
}

template <class P>
void BfvPolyDecrypt(std::array<uint64_t, P::n> &poly, const TRLWE<P> &ct,
                    const Key<P> &key)
{
    const Polynomial<P> phase = trlwePhase<P>(ct, key);
    for (uint32_t i = 0; i < P::n; i++)
        poly[i] = bfvDecodeCoeff<P>(phase[i]);
}

template <class FromP, class ToP>
void ConvertSecretKey(Key<ToP> &out, const Key<FromP> &in)
{
    static_assert(FromP::n == ToP::n && FromP::k == ToP::k,
                  "key conversion requires matching dimensions");
    for (uint32_t i = 0; i < FromP::k * FromP::n; i++)
        out[i] = static_cast<typename ToP::T>(in[i]);
}

template <class PlainP, class KeyP>
void SecretKeyAsPlaintext(std::array<uint64_t, PlainP::n> &out,
                          const Key<KeyP> &key)
{
    static_assert(PlainP::n == KeyP::n && PlainP::k == KeyP::k,
                  "key conversion requires matching dimensions");
    constexpr uint64_t q = static_cast<uint64_t>(PlainP::plain_modulus);
    for (uint32_t i = 0; i < PlainP::n; i++) {
        if (key[i] == static_cast<typename KeyP::T>(-1))
            out[i] = q - 1;
        else
            out[i] = static_cast<uint64_t>(key[i]) % q;
    }
}

template <class P>
void PlainPolynomialMul(TRLWE<P> &res, const TRLWE<P> &ct,
                        const std::array<uint64_t, P::n> &plain)
{
    Polynomial<P> p;
    for (uint32_t i = 0; i < P::n; i++)
        p[i] = static_cast<typename P::T>(
            plain[i] % static_cast<uint64_t>(P::plain_modulus));
    for (int c = 0; c <= static_cast<int>(P::k); c++) {
        if constexpr (is_multilimb_uint_v<typename P::T>)
            PolyMulTorusByUnsigned<P>(res[c], ct[c], plain);
        else
            PolyMul<P>(res[c], ct[c], p);
    }
}

template <class InP, class BootP>
void ModSwitchCiphertextToPlainPolys(std::array<uint64_t, BootP::n> &a,
                                     std::array<uint64_t, BootP::n> &b,
                                     const TRLWE<InP> &ct)
{
    static_assert(InP::n == BootP::n && InP::k == BootP::k,
                  "mod-switch requires matching dimensions");
    constexpr uint64_t q = static_cast<uint64_t>(BootP::plain_modulus);
    for (uint32_t i = 0; i < InP::n; i++) {
        a[i] = static_cast<uint64_t>(modSwitchCoeff<InP>(ct[0][i], q));
        b[i] = static_cast<uint64_t>(modSwitchCoeff<InP>(ct[InP::k][i], q));
    }
}

// Computes an encryption under BootP of the mod-switched noisy decryption
// b' - a' * s, where (a', b') are the input ciphertext components interpreted
// in Z/(BootP::plain_modulus) after BFV modulus switching.
template <class InP, class BootP>
void NoisyDecrypt(TRLWE<BootP> &res, const TRLWE<InP> &ct,
                  const TRLWE<BootP> &enc_sk)
{
    std::array<uint64_t, BootP::n> a{}, b{};
    ModSwitchCiphertextToPlainPolys<InP, BootP>(a, b, ct);

    auto prod = std::make_unique<TRLWE<BootP>>();
    PlainPolynomialMul<BootP>(*prod, enc_sk, a);

    for (int c = 0; c < static_cast<int>(BootP::k); c++)
        for (uint32_t i = 0; i < BootP::n; i++)
            res[c][i] = -(*prod)[c][i];
    for (uint32_t i = 0; i < BootP::n; i++)
        res[BootP::k][i] =
            bfvEncodeCoeff<BootP>(b[i]) - (*prod)[BootP::k][i];
}

template <class BaseP>
struct BootstrapKey {
    using BootP = PrimePower2Param<BaseP>;
    static constexpr uint64_t digit_error_bound =
        BaseP::bfv_bootstrap_digit_error_bound;
    static constexpr int default_linear_bsgs_step =
        BaseP::bfv_bootstrap_linear_bsgs_step;

    std::unique_ptr<GaloisKey<BaseP>> base_galois;
    std::unique_ptr<GaloisKey<BootP>> boot_galois;
    std::unique_ptr<relinKeyFFT<BootP>> relin;
    std::unique_ptr<TRLWE<BootP>> enc_sk;
    std::vector<uint64_t> digit_removal_polynomial;

    std::unique_ptr<std::vector<std::array<uint64_t, BaseP::n>>> stc_same;
    std::unique_ptr<std::vector<std::array<uint64_t, BaseP::n>>> stc_cross;
    std::unique_ptr<std::vector<std::array<uint64_t, BootP::n>>> cts_same;
    std::unique_ptr<std::vector<std::array<uint64_t, BootP::n>>> cts_cross;

    int linear_bsgs_step = default_linear_bsgs_step;
};

template <class BaseP>
BootstrapKey<BaseP> MakeBootstrapKey(const Key<BaseP> &base_key,
                                     bool build_linear_maps = true)
{
    using BootP = typename BootstrapKey<BaseP>::BootP;

    BootstrapKey<BaseP> bk;
    Key<BootP> boot_key;
    ConvertSecretKey<BaseP, BootP>(boot_key, base_key);

    bk.relin = makeRelinKeyFFT<BootP>(boot_key);

    std::array<uint64_t, BootP::n> sk_plain{};
    SecretKeyAsPlaintext<BootP, BaseP>(sk_plain, base_key);
    bk.enc_sk = std::make_unique<TRLWE<BootP>>();
    BfvPolyEncrypt<BootP>(*bk.enc_sk, sk_plain, boot_key);
    bk.digit_removal_polynomial =
        digitext::GetLowestDigitRemovalPolynomialOverRange(
            static_cast<uint64_t>(BaseP::plain_modulus),
            BootstrapKey<BaseP>::digit_error_bound);

    if (build_linear_maps) {
        bk.base_galois = std::make_unique<GaloisKey<BaseP>>();
        GaloisKeyGen<BaseP>(*bk.base_galois, base_key);

        bk.boot_galois = std::make_unique<GaloisKey<BootP>>();
        GaloisKeyGen<BootP>(*bk.boot_galois, boot_key);

        bk.stc_same =
            std::make_unique<std::vector<std::array<uint64_t, BaseP::n>>>();
        bk.stc_cross =
            std::make_unique<std::vector<std::array<uint64_t, BaseP::n>>>();
        c2s::BuildSlotToCoeffKeys<BaseP>(*bk.stc_same, *bk.stc_cross);

        bk.cts_same =
            std::make_unique<std::vector<std::array<uint64_t, BootP::n>>>();
        bk.cts_cross =
            std::make_unique<std::vector<std::array<uint64_t, BootP::n>>>();
        c2s::BuildCoeffToSlotKeys<BootP>(*bk.cts_same, *bk.cts_cross);
    }

    return bk;
}

// First online bootstrapping stage:
//   Enc_p(slots) -> SlotToCoeff -> mod-switched noisy decryption -> CoeffToSlot
//
// The result is a BootP ciphertext whose slots contain p*m + low-order noise
// in Z/p^2Z. The remaining digit-extraction stage will turn this into a
// refreshed BaseP ciphertext.
template <class BaseP>
void BootstrapNoisySlots(TRLWE<typename BootstrapKey<BaseP>::BootP> &res,
                         const TRLWE<BaseP> &ct,
                         const BootstrapKey<BaseP> &bk)
{
    using BootP = typename BootstrapKey<BaseP>::BootP;
    assert(bk.base_galois && bk.boot_galois && bk.stc_same && bk.stc_cross &&
           bk.cts_same && bk.cts_cross);

    auto coeff_ct = std::make_unique<TRLWE<BaseP>>();
    c2s::SlotToCoeff<BaseP>(*coeff_ct, ct, *bk.stc_same, *bk.stc_cross,
                            bk.linear_bsgs_step, *bk.base_galois);

    auto noisy_coeffs = std::make_unique<TRLWE<BootP>>();
    NoisyDecrypt<BaseP, BootP>(*noisy_coeffs, *coeff_ct, *bk.enc_sk);
    coeff_ct.reset();

    c2s::CoeffToSlot<BootP>(res, *noisy_coeffs, *bk.cts_same, *bk.cts_cross,
                            bk.linear_bsgs_step, *bk.boot_galois);
}

template <class BaseP>
void ProjectToBase(TRLWE<BaseP> &res,
                   const TRLWE<typename BootstrapKey<BaseP>::BootP> &ct)
{
    using BootP = typename BootstrapKey<BaseP>::BootP;
    static_assert(std::is_same_v<typename BaseP::T, typename BootP::T>,
                  "projection requires matching torus types");
    static_assert(BaseP::n == BootP::n && BaseP::k == BootP::k,
                  "projection requires matching ciphertext dimensions");

    for (int c = 0; c <= static_cast<int>(BaseP::k); c++)
        for (uint32_t i = 0; i < BaseP::n; i++)
            res[c][i] = static_cast<typename BaseP::T>(ct[c][i]);
}

// Final online bootstrapping stage:
//   Enc_{p^2}(p*m + e in slots) -> Enc_{p^2}(p*m) -> Enc_p(m)
//
// For BFV scaling, encoding p*m modulo p^2 uses the same torus phase as
// encoding m modulo p, so the final projection is a component-wise copy.
template <class BaseP>
void FinalizeBootstrap(TRLWE<BaseP> &res,
                       const TRLWE<typename BootstrapKey<BaseP>::BootP> &noisy,
                       const BootstrapKey<BaseP> &bk)
{
    using BootP = typename BootstrapKey<BaseP>::BootP;
    assert(bk.relin);
    assert(!bk.digit_removal_polynomial.empty());

    auto removed = std::make_unique<TRLWE<BootP>>();
    PolyEval<BootP>(*removed, bk.digit_removal_polynomial, noisy, *bk.relin);
    ProjectToBase<BaseP>(res, *removed);
}

template <class BaseP>
void Bootstrap(TRLWE<BaseP> &res, const TRLWE<BaseP> &ct,
               const BootstrapKey<BaseP> &bk)
{
    using BootP = typename BootstrapKey<BaseP>::BootP;
    auto noisy_slots = std::make_unique<TRLWE<BootP>>();
    BootstrapNoisySlots<BaseP>(*noisy_slots, ct, bk);
    FinalizeBootstrap<BaseP>(res, *noisy_slots, bk);
}

}  // namespace bfvboot
}  // namespace TFHEpp
