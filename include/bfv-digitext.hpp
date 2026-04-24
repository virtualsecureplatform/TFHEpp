#pragma once

#include <cassert>
#include <cstdint>
#include <vector>

#include "bfv-slots.hpp"

// bfv-digitext.hpp: Digit-extraction polynomial generators and orchestration.
//
// Mirrors the Magma proof-of-concept in
//   Bootstrapping_BGV_BFV/Digit extraction/{Arithmetic,Digit_extraction,Poly_eval}.m
//
// Polynomials operate in Z / p^e Z with centered-reduction semantics for odd p.
// For a prime power plaintext modulus t = p^e encoded in balanced digits,
// these polynomials drive the Halevi-Shoup / Chen-Han / "our" digit extraction.

namespace TFHEpp {
namespace digitext {

// ---------------------------------------------------------------------------
// Modular helpers (64-bit, to match the PoC which works with Magma integers)
// ---------------------------------------------------------------------------

inline uint64_t powmod(uint64_t base, uint64_t exp, uint64_t mod)
{
    uint64_t r = 1 % mod;
    base %= mod;
    while (exp > 0) {
        if (exp & 1) r = static_cast<uint64_t>((static_cast<unsigned __int128>(r) * base) % mod);
        base = static_cast<uint64_t>((static_cast<unsigned __int128>(base) * base) % mod);
        exp >>= 1;
    }
    return r;
}

// Modular inverse via extended Euclidean. Requires gcd(a, mod) = 1.
inline uint64_t modinv(uint64_t a, uint64_t mod)
{
    // Extended GCD on signed 128-bit to safely handle any 64-bit mod.
    __int128 old_r = static_cast<__int128>(a % mod);
    __int128 r = static_cast<__int128>(mod);
    __int128 old_s = 1, s = 0;
    while (r != 0) {
        __int128 q = old_r / r;
        __int128 tmp = old_r - q * r; old_r = r; r = tmp;
        tmp = old_s - q * s; old_s = s; s = tmp;
    }
    // old_r should be gcd; expect 1.
    assert(old_r == 1 || old_r == -1);
    if (old_r == -1) old_s = -old_s;
    old_s %= static_cast<__int128>(mod);
    if (old_s < 0) old_s += mod;
    return static_cast<uint64_t>(old_s);
}

inline uint64_t power(uint64_t base, uint64_t exp)
{
    uint64_t r = 1;
    for (uint64_t i = 0; i < exp; i++) r *= base;
    return r;
}

// Centered reduction mod p: result in (-p/2, p/2] (for odd p) or [-p/2, p/2) (for p=2).
inline int64_t centered_reduction(int64_t x, uint64_t p)
{
    int64_t r = x % static_cast<int64_t>(p);
    if (r < 0) r += p;
    if (p == 2) {
        // p=2 special case: use non-centered (PoC convention).
        return r;
    }
    const int64_t half = static_cast<int64_t>(p) / 2;
    if (r > half) r -= static_cast<int64_t>(p);
    return r;
}

// ---------------------------------------------------------------------------
// GetLowestDigitRemovalPolynomial(p, e)
//
// Returns coefficients (lowest to highest, mod p^e) of a polynomial f(x) such
// that f(x) ≡ x - (x mod p) (mod p^e) for all integers x (with balanced
// reduction of x mod p when p is odd).  Degree is (p-1)(e-1) + 1.
//
// Implementation follows the Magma PoC: build the value table f(i) for
// i = 0..degree, compute forward differences iteratively, then convert the
// Newton form into monomial coefficients.
// ---------------------------------------------------------------------------

inline std::vector<uint64_t> GetLowestDigitRemovalPolynomial(uint64_t p, uint64_t e)
{
    assert(p >= 2);
    assert(e >= 1);
    const uint64_t mod = power(p, e);
    const uint64_t degree = (p - 1) * (e - 1) + 1;

    // f(i) = i - (i mod p)  (with balanced reduction if p is odd)
    std::vector<int64_t> fd(degree + 1);
    for (uint64_t i = 0; i <= degree; i++) {
        int64_t signed_i = static_cast<int64_t>(i);
        int64_t reduced = (p == 2) ? (signed_i % 2 + (signed_i % 2 < 0 ? 2 : 0))
                                   : centered_reduction(signed_i, p);
        fd[i] = signed_i - reduced;
    }

    // Compute forward differences in place, storing Δ^k f(0), Δ^(k-1) f(1), ...
    // at positions [k, k+1, ..., degree] after k iterations (Magma's layout).
    // We follow the Magma PoC exactly: for iteration in 1..degree,
    //   for index = 0..degree-iteration: fd[degree-index] -= fd[degree-index-1]
    for (uint64_t iter = 1; iter <= degree; iter++) {
        for (uint64_t idx = 0; idx <= degree - iter; idx++) {
            int64_t &target = fd[degree - idx];
            target -= fd[degree - idx - 1];
            // Reduce mod p^e
            int64_t m = static_cast<int64_t>(mod);
            target %= m;
            if (target < 0) target += m;
        }
    }

    // Now fd[k] = (Δ^k f)(0) for k = 0..degree (by the same layout as the PoC).
    // Newton form: f(x) = Σ_{k=0..degree} Δ^k f(0) · x·(x-1)·...·(x-k+1) / k!
    // Convert to monomials mod p^e:
    //   poly = Σ_{k=0..degree} fd[k] / p_factors[k] * modinv(prod[k], p^e) * factor[k]
    //   where factor[k] = x·(x-1)·...·(x-k+1),
    //         prod[k]   = k! / (product of p^{valuation(i,p)} for i=1..k),
    //         p_factors[k] = product of p^{valuation(i,p)} for i=1..k
    //
    // The p_factors term exactly absorbs the p-valuation of k!, so prod[k] is
    // coprime to p and modinv is well-defined in Z/p^e Z.
    std::vector<int64_t> factor_poly = {1};  // factor = 1 initially
    std::vector<int64_t> result_poly(degree + 1, 0);
    uint64_t prod = 1;
    uint64_t p_factors = 1;
    for (uint64_t k = 0; k <= degree; k++) {
        // Term: fd[k] / p_factors * modinv(prod) * factor_poly
        uint64_t coeff = static_cast<uint64_t>(fd[k]) / p_factors;  // exact division
        coeff = static_cast<uint64_t>(
            (static_cast<unsigned __int128>(coeff) * modinv(prod, mod)) % mod);

        // Add coeff * factor_poly to result_poly
        for (size_t j = 0; j < factor_poly.size(); j++) {
            int64_t m = static_cast<int64_t>(mod);
            int64_t v = result_poly[j] + static_cast<int64_t>(
                (static_cast<unsigned __int128>(coeff) *
                 static_cast<unsigned __int128>(factor_poly[j] % m + (factor_poly[j] % m < 0 ? m : 0))) % mod);
            v %= m;
            if (v < 0) v += m;
            result_poly[j] = v;
        }

        // factor_poly *= (x - k)
        if (k < degree) {
            std::vector<int64_t> new_factor(factor_poly.size() + 1, 0);
            int64_t m = static_cast<int64_t>(mod);
            for (size_t j = 0; j < factor_poly.size(); j++) {
                new_factor[j] -= static_cast<int64_t>(k) * factor_poly[j];
                new_factor[j + 1] += factor_poly[j];
            }
            for (auto &v : new_factor) {
                v %= m;
                if (v < 0) v += m;
            }
            factor_poly = std::move(new_factor);

            // Update prod, p_factors for next iteration's k+1
            uint64_t next_k = k + 1;
            uint64_t val_p = 0;
            uint64_t tmp = next_k;
            while (tmp % p == 0) { val_p++; tmp /= p; }
            uint64_t p_v = power(p, val_p);
            prod = static_cast<uint64_t>(
                (static_cast<unsigned __int128>(prod) * (next_k / p_v)) % mod);
            p_factors *= p_v;
        }
    }

    std::vector<uint64_t> out(degree + 1);
    for (uint64_t i = 0; i <= degree; i++)
        out[i] = static_cast<uint64_t>(result_poly[i]);
    return out;
}

// ---------------------------------------------------------------------------
// GetLowestDigitRetainPolynomial(p, e)
//
// Returns coefficients of f(x) ≡ (x mod p) (mod p^e), i.e. f(x) = x - removal(x).
// For odd p, the resulting polynomial contains only odd-degree coefficients
// (matching the balanced reduction).  For p=2, only even-degree coefficients
// after the symmetrization step.
// ---------------------------------------------------------------------------

inline std::vector<uint64_t> GetLowestDigitRetainPolynomial(uint64_t p, uint64_t e)
{
    assert(p >= 2);
    assert(e >= 1);
    // Magma: if p=2, compute at precision e+1, then average (p=2) or antisymmetrize (p odd).
    const uint64_t e_new = (p == 2) ? (e + 1) : e;
    const uint64_t mod_new = power(p, e_new);
    const uint64_t mod = power(p, e);

    std::vector<uint64_t> removal = GetLowestDigitRemovalPolynomial(p, e_new);
    // retain(x) = x - removal(x), reduce mod p^e_new
    std::vector<int64_t> poly(std::max<size_t>(2, removal.size()), 0);
    for (size_t i = 0; i < removal.size(); i++) poly[i] = -static_cast<int64_t>(removal[i]);
    poly[1] += 1;
    for (auto &v : poly) {
        int64_t m = static_cast<int64_t>(mod_new);
        v %= m;
        if (v < 0) v += m;
    }

    // Symmetrize to keep only even (p=2) or odd (p odd) powers, then reduce mod p^e.
    std::vector<int64_t> sym(poly.size(), 0);
    if (p == 2) {
        // (poly(x) + poly(-x)) / 2 — keeps even degrees.
        for (size_t i = 0; i < poly.size(); i++) {
            int64_t sign = (i % 2 == 0) ? 1 : -1;
            sym[i] = (poly[i] + sign * poly[i]) / 2;
        }
    } else {
        // (poly(x) - poly(-x)) / 2 — keeps odd degrees.
        for (size_t i = 0; i < poly.size(); i++) {
            int64_t sign = (i % 2 == 0) ? 1 : -1;
            sym[i] = (poly[i] - sign * poly[i]) / 2;
        }
    }

    std::vector<uint64_t> out(sym.size());
    int64_t m = static_cast<int64_t>(mod);
    for (size_t i = 0; i < sym.size(); i++) {
        int64_t v = sym[i] % m;
        if (v < 0) v += m;
        out[i] = static_cast<uint64_t>(v);
    }
    // Trim trailing zeros so PolyEval sees the true degree.
    while (out.size() > 1 && out.back() == 0) out.pop_back();
    return out;
}

// ---------------------------------------------------------------------------
// Plaintext evaluation — useful for cross-checking polynomial generators.
// ---------------------------------------------------------------------------

inline uint64_t plainEvalMod(const std::vector<uint64_t> &coeffs, uint64_t x, uint64_t mod)
{
    unsigned __int128 result = 0, xi = 1;
    x %= mod;
    for (uint64_t c : coeffs) {
        result = (result + static_cast<unsigned __int128>(c) * xi) % mod;
        xi = (xi * x) % mod;
    }
    return static_cast<uint64_t>(result);
}

}  // namespace digitext
}  // namespace TFHEpp
