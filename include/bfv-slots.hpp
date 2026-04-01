#pragma once

#include <array>
#include <cassert>
#include <cstdint>
#include <limits>
#include <memory>
#include <type_traits>
#include <vector>

#include "bfv++.hpp"
#include "evalkeygens.hpp"
#include "keyswitch.hpp"
#include "params.hpp"
#include "trlwe.hpp"

// bfv-slots.hpp: SIMD slot operations for BFV in TFHEpp
//
// SIMD slots exploit the CRT factorization of Z_t[x]/(x^n+1) when t is prime
// and t ≡ 1 (mod 2n).  Under this condition x^n+1 splits into n linear factors,
// giving n independent "slots" automatically operated on in parallel by poly
// arithmetic.  Slot rotation requires Galois automorphisms (x → x^d) handled by
// the existing EvalAuto infrastructure.
//
// BFV scaling:
//   Δ = floor(Q/t) where Q = 2^128, t = 114689.
//   Encrypt: floor(m·Q/t) per coefficient (≤1 rounding error, independent of m).
//   Decrypt: round(phase·t/Q) (exact BFV decoding).
//   Multiply: accumulate digit products in 384-bit, divide by Δ at the end.

namespace TFHEpp {

// ---------------------------------------------------------------------------
// Modular arithmetic helpers (all operations mod 64-bit prime)
// ---------------------------------------------------------------------------

inline uint64_t mulmod64(uint64_t a, uint64_t b, uint64_t mod)
{
    return static_cast<uint64_t>(
        static_cast<unsigned __int128>(a) * b % mod);
}

inline uint64_t powmod64(uint64_t base, uint64_t exp, uint64_t mod)
{
    uint64_t result = 1;
    base %= mod;
    while (exp > 0) {
        if (exp & 1) result = mulmod64(result, base, mod);
        base = mulmod64(base, base, mod);
        exp >>= 1;
    }
    return result;
}

// ---------------------------------------------------------------------------
// Per-parameter SIMD constants
//   PSI      — primitive 2n-th root of unity mod t  (NTTmod pre-twist factor)
//   PSI_INV  — modular inverse of PSI
//   N_INV    — modular inverse of n mod t  (used in INTT normalization)
// ---------------------------------------------------------------------------

template <class P> struct SIMDConstants;  // primary template intentionally undefined

template <>
struct SIMDConstants<lvl3simdparam> {
    // t = 114689,  n = 4096,  2n = 8192
    // g = 3 (primitive root of t),  ψ = g^((t-1)/(2n)) = 3^14 mod 114689
    static constexpr uint64_t t     = 114689;
    static constexpr uint64_t PSI     = 80720;   // ψ: primitive 2n-th root mod t
    static constexpr uint64_t PSI_INV = 7887;    // ψ^{-1} mod t
    static constexpr uint64_t N_INV   = 114661;  // n^{-1} mod t  (n=4096)
};

// ---------------------------------------------------------------------------
// Galois permutation tables for slot rotation.
//
// The Galois group for x^n+1 with n=2^k is (Z/2n)* ≅ Z_2 × Z_{n/2}.
// Generator g=5 generates the Z_{n/2} factor.  The n NTT evaluation points
// split into two "rows" of n/2, indexed by powers of 5:
//   Row 0: NTT positions (5^i - 1)/2 mod n   for i=0..n/2-1
//   Row 1: NTT positions (2n - 5^i - 1)/2 mod n
//
// For slot rotation to work, the user's slot[i] must correspond to the
// i-th Galois orbit position (5^i mod 2n), NOT the i-th NTT position.
// The automorphism x→x^5 then cyclically rotates slot indices by 1.
//
// galois_ntt_index[i] = (5^i - 1)/2 mod n  (row 0)
// galois_ntt_index[i + n/2] = (2n - 5^i - 1)/2 mod n  (row 1)
// ---------------------------------------------------------------------------

template <class P>
struct GaloisPermutation {
    static const std::array<uint32_t, P::n> &ntt_index() {
        static const auto table = []() {
            std::array<uint32_t, P::n> tbl;
            constexpr uint64_t twon = 2 * P::n;
            uint64_t g = 1;  // 5^0
            for (uint32_t i = 0; i < P::n / 2; i++) {
                tbl[i]         = static_cast<uint32_t>((g - 1) / 2);
                tbl[i + P::n/2] = static_cast<uint32_t>((twon - g - 1) / 2);
                g = g * 5 % twon;
            }
            return tbl;
        }();
        return table;
    }
    // Inverse: given NTT position, return slot index
    static const std::array<uint32_t, P::n> &slot_index() {
        static const auto table = []() {
            std::array<uint32_t, P::n> tbl;
            const auto &fwd = ntt_index();
            for (uint32_t i = 0; i < P::n; i++)
                tbl[fwd[i]] = i;
            return tbl;
        }();
        return table;
    }
};

// ---------------------------------------------------------------------------
// In-place bit-reversal permutation (standard for Cooley-Tukey NTT)
// ---------------------------------------------------------------------------

template <uint32_t N>
static void bit_reverse(std::array<uint64_t, N> &a)
{
    static_assert((N & (N - 1)) == 0, "N must be a power of 2");
    constexpr int LOG = __builtin_ctz(N);
    for (uint32_t i = 0; i < N; i++) {
        uint32_t j = 0;
        uint32_t x = i;
        for (int b = 0; b < LOG; b++) {
            j = (j << 1) | (x & 1);
            x >>= 1;
        }
        if (j > i) std::swap(a[i], a[j]);
    }
}

// ---------------------------------------------------------------------------
// NTTmod<P> — forward twisted NTT mod t
//   Input:  array of n polynomial coefficients in [0, t)
//   Output: array of n slot evaluations at ψ^{2i+1} for i=0..n-1
//
// Algorithm:
//   1. Pre-twist:  a[i] *= ψ^i   (makes the ring negacyclic)
//   2. Bit-reversal permutation
//   3. Cooley-Tukey DIT butterfly with ω = ψ² as primitive n-th root
// ---------------------------------------------------------------------------

template <class P>
void NTTmod(std::array<uint64_t, P::n> &a)
{
    using C = SIMDConstants<P>;
    constexpr uint64_t t   = C::t;
    constexpr uint64_t psi = C::PSI;
    constexpr uint32_t n   = P::n;

    // Step 1: pre-twist a[i] *= ψ^i
    uint64_t psi_pow = 1;
    for (uint32_t i = 0; i < n; i++) {
        a[i] = mulmod64(a[i], psi_pow, t);
        psi_pow = mulmod64(psi_pow, psi, t);
    }

    // Step 2: bit-reversal
    bit_reverse<P::n>(a);

    // Step 3: Cooley-Tukey DIT butterflies with ω = ψ²
    // ω is the primitive n-th root of unity (ψ^2)
    const uint64_t omega = mulmod64(psi, psi, t);

    for (uint32_t len = 2; len <= n; len <<= 1) {
        // w_len = ω^(n/len): principal len-th root of unity
        const uint64_t w_len = powmod64(omega, n / len, t);
        for (uint32_t i = 0; i < n; i += len) {
            uint64_t w = 1;
            for (uint32_t j = 0; j < len / 2; j++) {
                const uint64_t u = a[i + j];
                const uint64_t v = mulmod64(a[i + j + len / 2], w, t);
                a[i + j]            = (u + v >= t) ? u + v - t : u + v;
                a[i + j + len / 2]  = (u >= v) ? u - v : u - v + t;
                w = mulmod64(w, w_len, t);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// INTTmod<P> — inverse twisted NTT mod t
//   Input:  array of n slot values in [0, t)
//   Output: array of n polynomial coefficients in [0, t)
//
// Algorithm:
//   1. Bit-reversal permutation
//   2. Cooley-Tukey DIT butterfly with ω_inv = ω^{-1} = (ψ²)^{-1} = ψ^{-2}
//   3. Normalize: multiply each element by n^{-1} mod t
//   4. Post-untwist: a[i] *= ψ^{-i}
// ---------------------------------------------------------------------------

template <class P>
void INTTmod(std::array<uint64_t, P::n> &a)
{
    using C = SIMDConstants<P>;
    constexpr uint64_t t       = C::t;
    constexpr uint64_t psi_inv = C::PSI_INV;
    constexpr uint64_t n_inv   = C::N_INV;
    constexpr uint32_t n       = P::n;

    // ω_inv = ψ^{-2}
    const uint64_t omega_inv = mulmod64(psi_inv, psi_inv, t);

    // Step 1: bit-reversal
    bit_reverse<P::n>(a);

    // Step 2: DIT butterfly with ω_inv (same structure as forward NTT)
    for (uint32_t len = 2; len <= n; len <<= 1) {
        const uint64_t w_len = powmod64(omega_inv, n / len, t);
        for (uint32_t i = 0; i < n; i += len) {
            uint64_t w = 1;
            for (uint32_t j = 0; j < len / 2; j++) {
                const uint64_t u = a[i + j];
                const uint64_t v = mulmod64(a[i + j + len / 2], w, t);
                a[i + j]            = (u + v >= t) ? u + v - t : u + v;
                a[i + j + len / 2]  = (u >= v) ? u - v : u - v + t;
                w = mulmod64(w, w_len, t);
            }
        }
    }

    // Step 3: normalize by n^{-1}
    for (uint32_t i = 0; i < n; i++)
        a[i] = mulmod64(a[i], n_inv, t);

    // Step 4: post-untwist a[i] *= ψ^{-i}
    uint64_t psi_inv_pow = 1;
    for (uint32_t i = 0; i < n; i++) {
        a[i] = mulmod64(a[i], psi_inv_pow, t);
        psi_inv_pow = mulmod64(psi_inv_pow, psi_inv, t);
    }
}

// ---------------------------------------------------------------------------
// SlotEncode<P> — slot vector → polynomial coefficients
//   slots[i] is a value in Z_t; output polynomial has coefficients in [0,t)
// ---------------------------------------------------------------------------

template <class P>
void SlotEncode(Polynomial<P> &poly, const std::array<uint64_t, P::n> &slots)
{
    // Apply Galois permutation: user slot[i] → NTT position galois_ntt_index[i]
    const auto &perm = GaloisPermutation<P>::ntt_index();
    std::array<uint64_t, P::n> ntt_order;
    for (uint32_t i = 0; i < P::n; i++)
        ntt_order[perm[i]] = slots[i];

    INTTmod<P>(ntt_order);  // inverse NTT: slot evaluations → polynomial coefficients
    for (uint32_t i = 0; i < P::n; i++)
        poly[i] = static_cast<typename P::T>(ntt_order[i]);
}

// ---------------------------------------------------------------------------
// SlotDecode<P> — polynomial coefficients → slot vector
// ---------------------------------------------------------------------------

template <class P>
void SlotDecode(std::array<uint64_t, P::n> &slots, const Polynomial<P> &poly)
{
    constexpr uint64_t t = static_cast<uint64_t>(P::plain_modulus);
    std::array<uint64_t, P::n> ntt_order;
    for (uint32_t i = 0; i < P::n; i++)
        ntt_order[i] = static_cast<uint64_t>(poly[i]) % t;
    NTTmod<P>(ntt_order);  // forward NTT: polynomial coefficients → slot evaluations

    // Apply inverse Galois permutation: NTT position → user slot[i]
    const auto &inv_perm = GaloisPermutation<P>::slot_index();
    for (uint32_t i = 0; i < P::n; i++)
        slots[inv_perm[i]] = ntt_order[i];
}

// ---------------------------------------------------------------------------
// trlweSlotEncrypt<P> — encode SIMD slots then encrypt as TRLWE
//
// Uses integer Δ (not double) to avoid precision loss when T=__uint128_t and
// t is an odd prime:  Δ_int = floor((2^width - 1) / t)
// The error introduced by floor vs exact division is < 1, far below the
// noise budget.
// ---------------------------------------------------------------------------

template <class P>
void trlweSlotEncrypt(TRLWE<P> &ct, const std::array<uint64_t, P::n> &slots,
                      const Key<P> &key)
{
    Polynomial<P> poly;
    SlotEncode<P>(poly, slots);

    // BFV encoding: floor(m · Q / t) per coefficient.
    // Decompose as: m·delta_int + floor(m·Q_mod_t / t)
    // where delta_int = floor(Q/t) and Q_mod_t = Q mod t.
    // This gives at most 1 unit of rounding error regardless of m.
    constexpr typename P::T delta = P::delta_int;
    constexpr uint64_t r = P::Q_mod_t;
    constexpr uint64_t t_val = static_cast<uint64_t>(P::plain_modulus);

    Polynomial<P> scaled;
    for (uint32_t i = 0; i < P::n; i++) {
        uint64_t m = static_cast<uint64_t>(poly[i]);
        // floor(m·Q/t) = m·delta_int + floor(m·r/t)
        // m < t < 2^17, r < t < 2^17, so m*r < 2^34 — fits in uint64_t
        scaled[i] = static_cast<typename P::T>(m) * delta
                  + static_cast<typename P::T>(m * r / t_val);
    }

    trlweSymEncrypt<P>(ct, scaled, key);
}

// ---------------------------------------------------------------------------
// trlweSlotDecrypt<P> — decrypt TRLWE then decode SIMD slots
//
// Uses integer arithmetic throughout to avoid double-precision loss for
// 128-bit torus values with a prime modulus.
// ---------------------------------------------------------------------------

// Compute round(phase · t / Q) mod t where Q = 2^128.
// This is the standard BFV decoding that works for any Δ = floor(Q/t).
// Returns a value in [0, t-1].
template <class P>
inline uint64_t bfvDecodeCoeff(typename P::T phase_u)
{
    constexpr uint64_t t_val = static_cast<uint64_t>(P::plain_modulus);
    // Compute round(phase · t / 2^128) = floor((phase · t + 2^127) / 2^128)
    //
    // phase · t is a 145-bit value (128 + 17 bits).
    // We compute the upper 64 bits of this 145-bit product.
    const uint64_t lo = static_cast<uint64_t>(phase_u);
    const uint64_t hi = static_cast<uint64_t>(phase_u >> 64);

    // phase · t = hi·t·2^64 + lo·t
    __uint128_t prod_lo = static_cast<__uint128_t>(lo) * t_val;
    __uint128_t prod_hi = static_cast<__uint128_t>(hi) * t_val;
    // Combine: add carry from prod_lo into prod_hi
    prod_hi += (prod_lo >> 64);
    // Add rounding bias: 2^127 >> 64 = 2^63
    prod_hi += static_cast<__uint128_t>(1) << 63;
    // Result = top 64 bits = (phase·t + 2^127) >> 128
    uint64_t result = static_cast<uint64_t>(prod_hi >> 64);
    return result % t_val;
}

template <class P>
void trlweSlotDecrypt(std::array<uint64_t, P::n> &slots, const TRLWE<P> &ct,
                      const Key<P> &key)
{
    Polynomial<P> phase = trlwePhase<P>(ct, key);

    Polynomial<P> poly;
    for (uint32_t i = 0; i < P::n; i++)
        poly[i] = static_cast<typename P::T>(bfvDecodeCoeff<P>(phase[i]));

    SlotDecode<P>(slots, poly);
}

// ---------------------------------------------------------------------------
// GaloisKey<P> — evaluation keys for all log2(n) power-of-2 slot rotations
//                plus one conjugation key.
//
//   index i  (0 ≤ i < nbit)  : key for rotation by 2^i slots
//                                Galois element = 5^{2^i} mod 2n
//   index nbit               : key for conjugation  (Galois element = 2n-1)
//
// Arbitrary rotation by k steps is applied via binary decomposition:
//   RotateSlots uses at most nbit sequential EvalAuto calls.
// ---------------------------------------------------------------------------

template <class P>
using GaloisKey = std::array<EvalAutoKey<P>, P::nbit + 1>;

// ---------------------------------------------------------------------------
// GaloisKeyGen<P> — generate the full GaloisKey from the secret key
// ---------------------------------------------------------------------------

template <class P>
void GaloisKeyGen(GaloisKey<P> &gk, const Key<P> &key)
{
    // Galois group generator g=5 for power-of-2 cyclotomics.
    // Rotation by 2^i uses element 5^{2^i} mod 2n.
    uint64_t d = 5;
    for (int i = 0; i < static_cast<int>(P::nbit); i++) {
        evalautokeygen<P>(gk[i], static_cast<uint>(d), key);
        d = d * d % (2 * P::n);  // 5^{2^{i+1}} mod 2n
    }
    // Conjugation: automorphism x → x^{-1} ≡ x^{2n-1}
    evalautokeygen<P>(gk[P::nbit], 2 * P::n - 1, key);
}

// ---------------------------------------------------------------------------
// RotateSlots<P> — rotate SIMD slots cyclically by 'steps' positions
//
// Uses binary decomposition of steps: applies at most nbit sequential
// EvalAuto calls, each using the corresponding power-of-2 rotation key.
// EvalAuto internally calls ExternalProduct which auto-selects the DD path
// when P::l̅ > 1.
// ---------------------------------------------------------------------------

template <class P>
void RotateSlots(TRLWE<P> &res, const TRLWE<P> &ct, int steps,
                 const GaloisKey<P> &gk)
{
    // The Galois group Z_{n/2} acts on each row of n/2 slots.
    // Rotation by k applies automorphism 5^k mod 2n.
    // Normalize steps to [0, n/2).
    constexpr int half = static_cast<int>(P::n) / 2;
    steps = ((steps % half) + half) % half;

    if (steps == 0) {
        res = ct;
        return;
    }

    // Binary decomposition of steps in [0, n/2):
    // bit i of steps → apply automorphism 5^{2^i} mod 2n (key gk[i]).
    // Only nbit-1 bits needed (since steps < n/2 = 2^{nbit-1}).
    TRLWE<P> cur = ct, tmp;
    uint64_t d = 5;
    for (int i = 0; i < static_cast<int>(P::nbit) - 1; i++) {
        if ((steps >> i) & 1) {
            EvalAuto<P>(tmp, cur, static_cast<uint>(d), gk[i]);
            cur = tmp;
        }
        d = d * d % (2 * P::n);
    }
    res = cur;
}

// ---------------------------------------------------------------------------
// ConjugateSlots<P> — apply the conjugation automorphism (x → x^{-1})
// ---------------------------------------------------------------------------

template <class P>
void ConjugateSlots(TRLWE<P> &res, const TRLWE<P> &ct, const GaloisKey<P> &gk)
{
    EvalAuto<P>(res, ct, 2 * P::n - 1, gk[P::nbit]);
}

// ---------------------------------------------------------------------------
// makeRelinKeyFFT<P> — heap-safe relinearization key generation
//
// For large parameter sets (e.g. lvl3simdparam) the standard relinKeyFFTgen
// creates a 4 MB relinKey + 2 MB relinKeyFFT on the stack, which overflows.
// This function avoids that by:
//   1. Allocating relinKeyFFT on the heap.
//   2. Calling halftrgswSymEncrypt(HalfTRGSWFFT<P>&, ...) which already uses
//      heap allocation internally for its HalfTRGSW<P> intermediate.
//   3. Returning the result as a unique_ptr.
//
// Note: HalfTRGSWFFT<P> and relinKeyFFT<P> are the same underlying type
// (aligned_array<TRLWEInFD<P>, P::l * P::l̅>), so the cast is valid.
// ---------------------------------------------------------------------------

template <class P>
std::unique_ptr<relinKeyFFT<P>> makeRelinKeyFFT(const Key<P> &key)
{
    // Compute s^2 (key polynomial squared) — same as in relinKeygen
    Polynomial<P> keysquare, partkey;
    for (int i = 0; i < static_cast<int>(P::n); i++) partkey[i] = key[i];
    PolyMulNaive<P>(keysquare, partkey, partkey);

    // relinKeyFFT<P> == HalfTRGSWFFT<P> (same type alias), allocate on heap.
    // halftrgswSymEncrypt(HalfTRGSWFFT<P>&, ...) uses make_unique<HalfTRGSW<P>>
    // internally, keeping the large intermediate off the stack.
    auto relinkeyfft = std::make_unique<relinKeyFFT<P>>();
    halftrgswSymEncrypt<P>(*relinkeyfft, keysquare, key);
    return relinkeyfft;
}

// ---------------------------------------------------------------------------
// makeBootKeyFFT<P> — heap-safe boot key generation for BFV bootstrapping
//
// Generates a HalfTRGSWFFT encrypting the secret key s (circular encryption).
// Uses gadget decomposition for noise-controlled multiplication in
// HomomorphicInnerProduct: noise grows as l·Bg·√n·σ instead of ||ct[0]||·σ.
// ---------------------------------------------------------------------------

template <class P>
std::unique_ptr<relinKeyFFT<P>> makeBootKeyFFT(const Key<P> &key)
{
    Polynomial<P> partkey;
    for (int i = 0; i < static_cast<int>(P::n); i++) partkey[i] = key[i];

    auto bootkeyfft = std::make_unique<relinKeyFFT<P>>();
    halftrgswSymEncrypt<P>(*bootkeyfft, partkey, key);
    return bootkeyfft;
}

// ---------------------------------------------------------------------------
// HomomorphicInnerProduct<P> — first step of BFV bootstrapping
//
// Given a BFV ciphertext ct = (a, b) where phase = b - a·s = Δ·m + noise,
// and a boot key (HalfTRGSWFFT encrypting s), computes a fresh TRLWE
// encrypting the phase of ct.
//
// TFHEpp convention: TRLWE = (ct[0], ct[1]) = (a, b), phase = b - a·s.
//
// The mask a = ct[0] is decomposed via the gadget and multiplied with the
// boot key (HalfTRGSW(s)) using relinKeySwitch. This keeps noise growth
// controlled regardless of the magnitude of a's coefficients.
// ---------------------------------------------------------------------------

template <class P>
void HomomorphicInnerProduct(TRLWE<P> &res, const TRLWE<P> &ct,
                             const relinKeyFFT<P> &bootkeyfft)
{
    // relinKeySwitch(ct[0], bootKey) → TRLWE(a·s)
    relinKeySwitch<P>(res, ct[0], bootkeyfft);
    // Subtract from trivial TRLWE(b):
    //   new_mask = -res_mask,  new_body = b - res_body
    // Phase: (b - res_body) - (-res_mask)·s = b - (res_body - res_mask·s) = b - a·s
    for (int k = 0; k < P::k; k++)
        for (uint32_t i = 0; i < P::n; i++)
            res[k][i] = -res[k][i];
    for (uint32_t i = 0; i < P::n; i++)
        res[P::k][i] = ct[P::k][i] - res[P::k][i];
}

// ---------------------------------------------------------------------------
// Baby-step/giant-step homomorphic polynomial evaluation
//
// Evaluates a cleartext polynomial f(X) = Σ coeffs[i]*X^i on an encrypted
// value x, producing Enc(f(x)).  All slots are evaluated in parallel (SIMD).
//
// Algorithm:
//   1. Choose k ≈ √degree, m = ⌈log₂(degree/k)⌉
//   2. Baby step: precompute x, x², ..., x^k  (k-1 multiplications)
//   3. Giant step: precompute x^k, x^{2k}, ..., x^{2^{m-1}·k}  (m-1 squarings)
//   4. Recursively evaluate using Horner-like decomposition at x^{k·2^j} boundaries
//
// Coefficients are unsigned integers mod t (or t²).  The constant term is
// encoded as a BFV plaintext (scaled by Δ) and added to the body.
// Non-constant terms are applied via scalar-ciphertext multiplication
// (multiply each TRLWE component by the scalar coefficient).
// ---------------------------------------------------------------------------

namespace polyeval {

// Find optimal baby-step size k and giant-step depth m for degree d.
// Minimizes total non-scalar multiplications ≈ k + m + (recursive overhead).
inline std::pair<int, int> FindBSGSParams(int d)
{
    if (d <= 0) return {1, 0};
    int best_k = 1, best_m = 0, best_cost = d + 1;
    for (int k = 1; k <= d; k++) {
        int m = 0;
        { int reach = k; while (reach < d) { reach *= 2; m++; } }
        // Cost: (k-1) baby muls + (m-1 if m>0) giant squarings + m recursive muls
        int cost = (k - 1) + (m > 0 ? m - 1 : 0) + m;
        if (cost < best_cost) {
            best_cost = cost; best_k = k; best_m = m;
        }
        if (k > 2 * best_k) break;  // prune search
    }
    return {best_k, best_m};
}

// Largest power of 2 ≤ n (for balanced binary splitting in baby-step sum).
inline int FloorPow2(int n)
{
    int p = 1;
    while (2 * p <= n) p *= 2;
    return p;
}

// Recursive evaluation of coefficient vector using precomputed powers.
// coeffs[0..len-1] are the coefficients (index 0 = lowest power).
// baby[1..k] = x^i.  giant[0..m-1] = x^{k·2^j}.
// level = current giant-step level (counts down from m to 0).
// Returns TRLWE encrypting Σ coeffs[i] * x^i.
template <class P>
void EvalRecursive(TRLWE<P> &res,
                   const uint64_t *coeffs, int len,
                   const std::vector<TRLWE<P>> &baby,
                   const std::vector<TRLWE<P>> &giant,
                   int level, int k,
                   const relinKeyFFT<P> &rlk)
{
    if (len <= 0) {
        // Zero
        for (int c = 0; c <= static_cast<int>(P::k); c++)
            for (uint32_t i = 0; i < P::n; i++)
                res[c][i] = 0;
        return;
    }

    if (level == 0) {
        // Baby-step inner product: res = Σ_{i=0}^{len-1} coeffs[i] * x^i
        // coeffs[0] is the constant term (add Δ*c to body).
        // coeffs[i>0] multiply baby[i] by the scalar.
        bool first = true;
        for (int i = 0; i < len; i++) {
            if (coeffs[i] == 0) continue;
            if (i == 0) {
                // Constant term: trivial TRLWE(Δ*c)
                if (first) {
                    for (int c = 0; c <= static_cast<int>(P::k); c++)
                        for (uint32_t j = 0; j < P::n; j++)
                            res[c][j] = 0;
                    constexpr typename P::T delta = P::delta_int;
                    constexpr uint64_t r = P::Q_mod_t;
                    constexpr uint64_t t_val = static_cast<uint64_t>(P::plain_modulus);
                    res[P::k][0] = static_cast<typename P::T>(coeffs[0]) * delta
                                 + static_cast<typename P::T>(coeffs[0] * r / t_val);
                    first = false;
                } else {
                    constexpr typename P::T delta = P::delta_int;
                    constexpr uint64_t r = P::Q_mod_t;
                    constexpr uint64_t t_val = static_cast<uint64_t>(P::plain_modulus);
                    res[P::k][0] += static_cast<typename P::T>(coeffs[0]) * delta
                                  + static_cast<typename P::T>(coeffs[0] * r / t_val);
                }
            } else {
                // Non-constant: scalar * baby[i]
                auto c_val = static_cast<typename P::T>(coeffs[i]);
                if (first) {
                    for (int c = 0; c <= static_cast<int>(P::k); c++)
                        for (uint32_t j = 0; j < P::n; j++)
                            res[c][j] = baby[i][c][j] * c_val;
                    first = false;
                } else {
                    for (int c = 0; c <= static_cast<int>(P::k); c++)
                        for (uint32_t j = 0; j < P::n; j++)
                            res[c][j] += baby[i][c][j] * c_val;
                }
            }
        }
        if (first) {
            // All coefficients were zero
            for (int c = 0; c <= static_cast<int>(P::k); c++)
                for (uint32_t j = 0; j < P::n; j++)
                    res[c][j] = 0;
        }
        return;
    }

    // Giant-step split: split at index k * 2^(level-1)
    int split = k;
    for (int i = 0; i < level - 1; i++) split *= 2;
    if (split > len) split = len;

    // Lower part: coeffs[0..split-1], recurse with level-1
    TRLWE<P> lower;
    EvalRecursive<P>(lower, coeffs, split, baby, giant, level - 1, k, rlk);

    if (split >= len) {
        // No upper part
        res = lower;
        return;
    }

    // Upper part: coeffs[split..len-1], recurse with level-1
    TRLWE<P> upper;
    EvalRecursive<P>(upper, coeffs + split, len - split, baby, giant, level - 1, k, rlk);

    // res = lower + upper * giant[level-1]
    TRLWE<P> prod;
    TRLWEMultFullDD<P>(prod, upper, giant[level - 1], rlk);
    for (int c = 0; c <= static_cast<int>(P::k); c++)
        for (uint32_t i = 0; i < P::n; i++)
            res[c][i] = lower[c][i] + prod[c][i];
}

}  // namespace polyeval

// ---------------------------------------------------------------------------
// PolyEval<P> — Evaluate cleartext polynomial on encrypted value
//
//   coeffs: f(X) = coeffs[0] + coeffs[1]*X + ... + coeffs[d]*X^d
//           Each coefficient is in [0, t) or [0, t²) (unsigned).
//   x:      TRLWE encrypting the input value (all slots evaluated in parallel)
//   rlk:    Relinearization key for ciphertext-ciphertext multiplication
//   res:    TRLWE encrypting f(x)
// ---------------------------------------------------------------------------

template <class P>
void PolyEval(TRLWE<P> &res,
              const std::vector<uint64_t> &coeffs,
              const TRLWE<P> &x,
              const relinKeyFFT<P> &rlk)
{
    int d = static_cast<int>(coeffs.size()) - 1;
    if (d <= 0) {
        // Constant or empty
        for (int c = 0; c <= static_cast<int>(P::k); c++)
            for (uint32_t i = 0; i < P::n; i++)
                res[c][i] = 0;
        if (!coeffs.empty() && coeffs[0] != 0) {
            constexpr typename P::T delta = P::delta_int;
            constexpr uint64_t r = P::Q_mod_t;
            constexpr uint64_t t_val = static_cast<uint64_t>(P::plain_modulus);
            res[P::k][0] = static_cast<typename P::T>(coeffs[0]) * delta
                         + static_cast<typename P::T>(coeffs[0] * r / t_val);
        }
        return;
    }

    auto [k, m] = polyeval::FindBSGSParams(d);

    // Baby step: baby[i] = x^i for i=1..k
    std::vector<TRLWE<P>> baby(k + 1);
    baby[1] = x;
    for (int i = 2; i <= k; i++) {
        int a = polyeval::FloorPow2(i - 1);  // balanced split for minimal depth
        int b = i - a;
        TRLWEMultFullDD<P>(baby[i], baby[a], baby[b], rlk);
    }

    // Giant step: giant[j] = x^{k*2^j} for j=0..m-1
    std::vector<TRLWE<P>> giant(m);
    if (m > 0) {
        giant[0] = baby[k];
        for (int j = 1; j < m; j++)
            TRLWEMultFullDD<P>(giant[j], giant[j - 1], giant[j - 1], rlk);
    }

    // Recursive evaluation
    polyeval::EvalRecursive<P>(res, coeffs.data(), d + 1,
                               baby, giant, m, k, rlk);
}

}  // namespace TFHEpp
