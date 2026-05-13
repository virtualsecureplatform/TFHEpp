#pragma once

#include <array>
#include <cstdint>
#include <memory>
#include <vector>

#include "bfv-slots.hpp"

// bfv-c2s.hpp — CoeffToSlot / SlotToCoeff linear transforms for BFV
// bootstrapping.
//
// The target slot-domain matrix is the SlotEncode matrix S of the current
// parameter set (S[i][j] = SlotEncode(e_j)[i] where e_j is the j-th standard
// basis slot vector).  S captures all of our NTT conventions — PSI pre/post
// twist, bit-reversal, and the GaloisPermutation between user slot index
// and NTT position — so deriving S by direct computation avoids having to
// re-derive those conventions from the Magma PoC.
//
// S has up to n non-zero diagonals in general.  We split S into two halves:
//
//   S_same[i][j]  = S[i][j]  when i and j are in the same Galois row,
//                  0 otherwise.
//   S_cross[i][j] = S[i][j]  when i and j are in different Galois rows,
//                  0 otherwise.
//
// Then S · x = S_same · x  +  S_cross · x.  The cross part can be written
// as M · σ(x) where σ is the conjugation automorphism (swap rows) and M is
// another block-diagonal (intra-row-only) matrix.  So
//
//   S · x  =  intra_row_LT(x; D_same)  +  intra_row_LT(σ(x); D_cross).
//
// Each intra-row LinearTransform costs ≤ 2√(n/2) ≈ 90 rotations via BSGS,
// plus one ConjugateSlots.  Total ~180 rotations per CoeffToSlot call for
// n=4096 — not optimal (FFT-factorization would give O(log n) = 12), but
// correct and a reasonable first implementation.

namespace TFHEpp {
namespace c2s {

namespace detail {

inline uint64_t modinv_u64(uint64_t a, uint64_t mod)
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
    old_s %= static_cast<__int128>(mod);
    if (old_s < 0) old_s += mod;
    return static_cast<uint64_t>(old_s);
}

template <class P>
uint64_t crt_stage_twiddle(int stage, uint32_t local_index, bool inverse_row)
{
    constexpr uint64_t t = static_cast<uint64_t>(P::plain_modulus);
    constexpr uint64_t twon = 2 * P::n;
    uint64_t exp = 1;
    for (uint32_t i = 0; i < local_index; i++) exp = (exp * 5) % twon;
    exp = (exp * (uint64_t{1} << stage)) % twon;
    uint64_t tw = powmod64(SIMDConstants<P>::PSI, exp, t);
    return inverse_row ? modinv_u64(tw, t) : tw;
}

template <class P>
void ApplyCRTRowButterfly(TRLWE<P> &res, const TRLWE<P> &ct, int stage,
                          bool inverse, const GaloisKey<P> &gk,
                          bool normalize_inverse = true)
{
    constexpr uint64_t t = static_cast<uint64_t>(P::plain_modulus);
    constexpr int half = static_cast<int>(P::n) / 2;
    const int stride = half >> (stage + 1);
    const uint64_t inv2 = modinv_u64(2, t);

    auto d0 = std::make_unique<std::array<uint64_t, P::n>>();
    auto d_plus = std::make_unique<std::array<uint64_t, P::n>>();
    auto d_minus = std::make_unique<std::array<uint64_t, P::n>>();
    d0->fill(0);
    d_plus->fill(0);
    d_minus->fill(0);

    for (int row = 0; row < 2; row++) {
        const int row_base = row * half;
        for (int block = 0; block < half; block += 2 * stride) {
            for (int i = 0; i < stride; i++) {
                const uint32_t first = row_base + block + i;
                const uint32_t second = first + stride;
                const uint64_t tw = crt_stage_twiddle<P>(
                    stage, static_cast<uint32_t>(i), row == 1);
                if (!inverse) {
                    (*d0)[first] = 1;
                    (*d0)[second] = (t - tw) % t;
                    (*d_plus)[first] = tw;
                    (*d_minus)[second] = 1;
                } else if (normalize_inverse) {
                    const uint64_t tw_inv = modinv_u64(tw, t);
                    (*d0)[first] = inv2;
                    (*d0)[second] = (t - mulmod64(tw_inv, inv2, t)) % t;
                    (*d_plus)[first] = inv2;
                    (*d_minus)[second] = mulmod64(tw_inv, inv2, t);
                } else {
                    const uint64_t tw_inv = modinv_u64(tw, t);
                    (*d0)[first] = 1;
                    (*d0)[second] = (t - tw_inv) % t;
                    (*d_plus)[first] = 1;
                    (*d_minus)[second] = tw_inv;
                }
            }
        }
    }

    std::vector<std::array<uint64_t, P::n>> diagonals;
    diagonals.reserve(3);
    diagonals.push_back(*d0);
    diagonals.push_back(*d_plus);
    diagonals.push_back(*d_minus);
    const std::vector<int> offsets{0, stride, -stride};
    LinearTransform<P>(res, ct, diagonals, offsets, gk);
}

template <class P>
void ApplyCRTCrossButterfly(TRLWE<P> &res, const TRLWE<P> &ct, bool inverse,
                            const GaloisKey<P> &gk,
                            bool normalize_inverse = true)
{
    constexpr uint64_t t = static_cast<uint64_t>(P::plain_modulus);
    constexpr int half = static_cast<int>(P::n) / 2;
    const uint64_t a = powmod64(SIMDConstants<P>::PSI, half, t);
    const uint64_t b = modinv_u64(a, t);

    auto d_id = std::make_unique<std::array<uint64_t, P::n>>();
    auto d_conj = std::make_unique<std::array<uint64_t, P::n>>();
    if (!inverse) {
        for (int i = 0; i < half; i++) {
            (*d_id)[i] = 1;
            (*d_conj)[i] = a;
            (*d_id)[half + i] = b;
            (*d_conj)[half + i] = 1;
        }
    } else if (normalize_inverse) {
        const uint64_t det = (b + t - a) % t;
        const uint64_t inv_det = modinv_u64(det, t);
        const uint64_t neg_a_over_det = mulmod64((t - a) % t, inv_det, t);
        const uint64_t neg_one_over_det = (t - inv_det) % t;
        for (int i = 0; i < half; i++) {
            (*d_id)[i] = mulmod64(b, inv_det, t);
            (*d_conj)[i] = neg_a_over_det;
            (*d_id)[half + i] = inv_det;
            (*d_conj)[half + i] = neg_one_over_det;
        }
    } else {
        for (int i = 0; i < half; i++) {
            (*d_id)[i] = b;
            (*d_conj)[i] = (t - a) % t;
            (*d_id)[half + i] = 1;
            (*d_conj)[half + i] = t - 1;
        }
    }

    auto id_part = std::make_unique<TRLWE<P>>();
    auto conj_ct = std::make_unique<TRLWE<P>>();
    auto conj_part = std::make_unique<TRLWE<P>>();
    SlotPtxtMul<P>(*id_part, ct, *d_id);
    ConjugateSlots<P>(*conj_ct, ct, gk);
    SlotPtxtMul<P>(*conj_part, *conj_ct, *d_conj);
    for (int k = 0; k <= static_cast<int>(P::k); k++)
        for (uint32_t i = 0; i < P::n; i++)
            res[k][i] = (*id_part)[k][i] + (*conj_part)[k][i];
}

template <class P>
void CRTSlotToCoeffProduct(TRLWE<P> &res, const TRLWE<P> &ct,
                           const GaloisKey<P> &gk)
{
    auto cur = std::make_unique<TRLWE<P>>(ct);
    auto next = std::make_unique<TRLWE<P>>();

    ApplyCRTCrossButterfly<P>(*next, *cur, false, gk);
    std::swap(cur, next);

    for (int stage = static_cast<int>(P::nbit) - 2; stage >= 0; stage--) {
        ApplyCRTRowButterfly<P>(*next, *cur, stage, false, gk);
        std::swap(cur, next);
    }
    res = *cur;
}

template <class P>
void CRTCoeffToSlotProduct(TRLWE<P> &res, const TRLWE<P> &ct,
                           const GaloisKey<P> &gk)
{
    constexpr uint64_t t = static_cast<uint64_t>(P::plain_modulus);
    constexpr int half = static_cast<int>(P::n) / 2;
    auto cur = std::make_unique<TRLWE<P>>(ct);
    auto next = std::make_unique<TRLWE<P>>();

    for (int stage = 0; stage <= static_cast<int>(P::nbit) - 2; stage++) {
        ApplyCRTRowButterfly<P>(*next, *cur, stage, true, gk, false);
        std::swap(cur, next);
    }
    ApplyCRTCrossButterfly<P>(*next, *cur, true, gk, false);

    const uint64_t a = powmod64(SIMDConstants<P>::PSI, half, t);
    const uint64_t b = modinv_u64(a, t);
    uint64_t scale = (b + t - a) % t;
    for (int i = 0; i < static_cast<int>(P::nbit) - 1; i++)
        scale = mulmod64(scale, 2, t);
    const uint64_t scale_inv = modinv_u64(scale, t);

    auto scale_slots = std::make_unique<std::array<uint64_t, P::n>>();
    scale_slots->fill(scale_inv);
    SlotPtxtMul<P>(res, *next, *scale_slots);
}

}  // namespace detail

// ---------------------------------------------------------------------------
// Compute the n × n SlotEncode matrix S where S[i][j] = SlotEncode(e_j)[i].
// Entries are stored mod t (uint64_t).  Roughly n · n entries = 16M uint64 =
// 128 MB for n = 4096 — caller should allocate on the heap.
// ---------------------------------------------------------------------------

template <class P>
void ComputeSlotEncodeMatrix(std::vector<std::array<uint64_t, P::n>> &S)
{
    constexpr uint64_t t = static_cast<uint64_t>(P::plain_modulus);
    S.assign(P::n, std::array<uint64_t, P::n>{});

    std::array<uint64_t, P::n> e_j{};
    Polynomial<P> poly;
    for (uint32_t j = 0; j < P::n; j++) {
        if (j > 0) e_j[j - 1] = 0;
        e_j[j] = 1;
        SlotEncode<P>(poly, e_j);
        for (uint32_t i = 0; i < P::n; i++)
            S[i][j] = static_cast<uint64_t>(poly[i]) % t;
    }
}

// ---------------------------------------------------------------------------
// Apply M ∈ Z_t^{n×n} directly to a length-n vector (reference implementation).
// ---------------------------------------------------------------------------

template <class P>
void MatVecMulMod(std::array<uint64_t, P::n> &y,
                  const std::vector<std::array<uint64_t, P::n>> &M,
                  const std::array<uint64_t, P::n> &x)
{
    constexpr uint64_t t = static_cast<uint64_t>(P::plain_modulus);
    for (uint32_t i = 0; i < P::n; i++) {
        unsigned __int128 acc = 0;
        for (uint32_t j = 0; j < P::n; j++)
            acc = (acc + static_cast<unsigned __int128>(M[i][j]) * x[j]) % t;
        y[i] = static_cast<uint64_t>(acc);
    }
}

// ---------------------------------------------------------------------------
// Diagonal decomposition of the "same-row" part of S.
//
// For a Galois row h ∈ {0, 1} and an intra-row offset r ∈ [0, n/2),
// the diagonal vector D_same[r][i] stores the coefficient that multiplies
// x[i → row h, shifted by r within the row]:
//
//   D_same[r][i] = S[i][i_same_row_shift_by_r]
//
// where i_same_row_shift_by_r = h·(n/2) + ((i_local + r) mod n/2) if
// i = h·(n/2) + i_local, so the input slot is in the SAME row as i.
//
// When we apply LinearTransform(x; D_same, {r}), each output slot i gets
//   y[i] = Σ_r D_same[r][i] · x[rotate_same_row(x, r)[i]]
//        = Σ_r S[i][rotate_same_row_j_from_i(r)] · x[j]
//   which is exactly the same-row portion of Σ_j S[i][j] · x[j].
// ---------------------------------------------------------------------------

template <class P>
void ExtractSameRowDiagonals(std::vector<std::array<uint64_t, P::n>> &D_same,
                             const std::vector<std::array<uint64_t, P::n>> &S)
{
    constexpr int half = static_cast<int>(P::n) / 2;
    D_same.assign(half, std::array<uint64_t, P::n>{});
    for (int r = 0; r < half; r++) {
        for (int i = 0; i < static_cast<int>(P::n); i++) {
            int i_local, row_base;
            if (i < half) { i_local = i; row_base = 0; }
            else          { i_local = i - half; row_base = half; }
            const int j = row_base + (i_local + r) % half;
            D_same[r][i] = S[i][j];
        }
    }
}

// ---------------------------------------------------------------------------
// Diagonal decomposition of the "cross-row" part of S, expressed as diagonals
// to be applied to σ(x) (the conjugated ciphertext).  σ swaps rows:
//   σ(x)[i] = x[i + half]  for i < half,
//   σ(x)[i] = x[i - half]  for i ≥ half.
//
// So for an output row 0 slot i, the cross-row contribution to y[i] from
// inputs in row 1 is Σ_{j ≥ half} S[i][j] · x[j] = Σ_{k < half} S[i][k+half]
// · σ(x)[k].  After conjugation, k indexes "opposite row", so an offset r
// within σ(x) corresponds to mapping position (k+r) mod half in σ(x), which
// is the INPUT row's slot (j - half + r) mod half + half — but this is still
// a same-row shift in the conjugated vector, so D_cross uses the same intra-
// row offset indexing as D_same (applied after ConjugateSlots).
// ---------------------------------------------------------------------------

template <class P>
void ExtractCrossRowDiagonals(std::vector<std::array<uint64_t, P::n>> &D_cross,
                              const std::vector<std::array<uint64_t, P::n>> &S)
{
    constexpr int half = static_cast<int>(P::n) / 2;
    D_cross.assign(half, std::array<uint64_t, P::n>{});
    for (int r = 0; r < half; r++) {
        for (int i = 0; i < static_cast<int>(P::n); i++) {
            int i_local, other_row_base;
            if (i < half) { i_local = i; other_row_base = half; }
            else          { i_local = i - half; other_row_base = 0; }
            // The slot index of σ(x) that RotateSlots_intra maps to position i
            // after a rotation by r is (i_local + r) mod half in the OPPOSITE
            // row of σ's view; after conjugation that corresponds to input slot
            // (i_local + r) mod half in the opposite row of x.
            const int j = other_row_base + (i_local + r) % half;
            D_cross[r][i] = S[i][j];
        }
    }
}

// ---------------------------------------------------------------------------
// Plaintext application of the decomposition:
//   y = ApplySameRow(x; D_same) + ApplyCrossRow(Conj(x); D_cross)
//
// ApplySameRow computes Σ_r D_same[r][i] · RotateIntra(x, r)[i] per output
// slot i.  This is exactly what LinearTransform(x; D_same, {0..half-1})
// would compute homomorphically.  Likewise for ApplyCrossRow, but the input
// to ApplyCrossRow is σ(x) (row-swapped).
//
// Combined correctness: ApplySameRow + ApplyCrossRow ∘ σ = S · x.
// ---------------------------------------------------------------------------

template <class P>
static inline void RotateIntraSlot(std::array<uint64_t, P::n> &out,
                                   const std::array<uint64_t, P::n> &in,
                                   int r)
{
    constexpr int half = static_cast<int>(P::n) / 2;
    r = ((r % half) + half) % half;
    for (int i = 0; i < half; i++) {
        out[i]        = in[(i + r) % half];
        out[half + i] = in[half + (i + r) % half];
    }
}

template <class P>
static inline void ConjSlot(std::array<uint64_t, P::n> &out,
                            const std::array<uint64_t, P::n> &in)
{
    constexpr int half = static_cast<int>(P::n) / 2;
    for (int i = 0; i < half; i++) {
        out[i]        = in[half + i];
        out[half + i] = in[i];
    }
}

// Plaintext "LinearTransform": Σ_r D[r] pointwise RotateIntra(x, r).
// Mirrors what LinearTransform<P> would compute homomorphically on a
// ciphertext whose slots equal x.
template <class P>
void PlainLinearTransform(std::array<uint64_t, P::n> &y,
                          const std::array<uint64_t, P::n> &x,
                          const std::vector<std::array<uint64_t, P::n>> &D)
{
    constexpr uint64_t t = static_cast<uint64_t>(P::plain_modulus);
    std::array<uint64_t, P::n> rotated{};
    y.fill(0);
    for (size_t r = 0; r < D.size(); r++) {
        RotateIntraSlot<P>(rotated, x, static_cast<int>(r));
        for (uint32_t i = 0; i < P::n; i++) {
            unsigned __int128 v = static_cast<unsigned __int128>(D[r][i]) * rotated[i];
            y[i] = static_cast<uint64_t>((y[i] + v) % t);
        }
    }
}

// Plaintext evaluation of the full decomposition:
//   y = PlainLinearTransform(x; D_same) + PlainLinearTransform(Conj(x); D_cross)
template <class P>
void PlainApplyS(std::array<uint64_t, P::n> &y,
                 const std::array<uint64_t, P::n> &x,
                 const std::vector<std::array<uint64_t, P::n>> &D_same,
                 const std::vector<std::array<uint64_t, P::n>> &D_cross)
{
    constexpr uint64_t t = static_cast<uint64_t>(P::plain_modulus);
    std::array<uint64_t, P::n> y_same{}, y_cross{}, conj_x{};
    PlainLinearTransform<P>(y_same, x, D_same);
    ConjSlot<P>(conj_x, x);
    PlainLinearTransform<P>(y_cross, conj_x, D_cross);
    for (uint32_t i = 0; i < P::n; i++)
        y[i] = (y_same[i] + y_cross[i]) % t;
}

// ---------------------------------------------------------------------------
// Build both diagonal decompositions for CoeffToSlot from the precomputed S.
// Convenience wrapper; allocates two heap vectors of n diagonals.
// ---------------------------------------------------------------------------

template <class P>
void BuildCoeffToSlotKeys(std::vector<std::array<uint64_t, P::n>> &D_same,
                          std::vector<std::array<uint64_t, P::n>> &D_cross)
{
    auto S = std::make_unique<std::vector<std::array<uint64_t, P::n>>>();
    ComputeSlotEncodeMatrix<P>(*S);
    ExtractSameRowDiagonals<P>(D_same, *S);
    ExtractCrossRowDiagonals<P>(D_cross, *S);
}

// ---------------------------------------------------------------------------
// Compute the n × n SlotDecode matrix S^{-1} where
//   S^{-1}[i][j] = SlotDecode(e_j)[i]
// with e_j the standard basis polynomial (coefficient 1 at position j).
// This is the inverse of the SlotEncode matrix and the target for
// SlotToCoeff.
// ---------------------------------------------------------------------------

template <class P>
void ComputeSlotDecodeMatrix(std::vector<std::array<uint64_t, P::n>> &Sinv)
{
    constexpr uint64_t t = static_cast<uint64_t>(P::plain_modulus);
    Sinv.assign(P::n, std::array<uint64_t, P::n>{});

    Polynomial<P> e_j{};
    std::array<uint64_t, P::n> slots{};
    for (uint32_t j = 0; j < P::n; j++) {
        if (j > 0) e_j[j - 1] = 0;
        e_j[j] = 1;
        SlotDecode<P>(slots, e_j);
        for (uint32_t i = 0; i < P::n; i++)
            Sinv[i][j] = slots[i] % t;
    }
}

// ---------------------------------------------------------------------------
// Build both diagonal decompositions for SlotToCoeff from S^{-1}.
// ---------------------------------------------------------------------------

template <class P>
void BuildSlotToCoeffKeys(std::vector<std::array<uint64_t, P::n>> &D_same,
                          std::vector<std::array<uint64_t, P::n>> &D_cross)
{
    auto Sinv = std::make_unique<std::vector<std::array<uint64_t, P::n>>>();
    ComputeSlotDecodeMatrix<P>(*Sinv);
    ExtractSameRowDiagonals<P>(D_same, *Sinv);
    ExtractCrossRowDiagonals<P>(D_cross, *Sinv);
}

// ---------------------------------------------------------------------------
// Homomorphic CoeffToSlot: applies S (the SlotEncode matrix) in slot domain
// to a TRLWE ciphertext.
//
// Given TRLWE ct whose polynomial has BFV-plaintext coefficients (p_0, ...,
// p_{n-1}) (i.e. phase ≈ Δ · P), the resulting TRLWE encrypts a polynomial
// whose user slots are exactly (p_0, ..., p_{n-1}).
//
// Implementation:
//   result = LinearTransformBSGS(ct,         D_same,  offsets=[0..n/2), k)
//          + LinearTransformBSGS(Conj(ct),   D_cross, offsets=[0..n/2), k)
//
// k is the baby-step factor; k ≈ √(n/2) is near-optimal.  Caller must provide
// the galois key; the conjugation key at index P::nbit is used for Conj.
// ---------------------------------------------------------------------------

template <class P>
void CoeffToSlot(TRLWE<P> &res, const TRLWE<P> &ct,
                 const std::vector<std::array<uint64_t, P::n>> &D_same,
                 const std::vector<std::array<uint64_t, P::n>> &D_cross,
                 int k_step,
                 const GaloisKey<P> &gk)
{
    constexpr int half = static_cast<int>(P::n) / 2;
    assert(static_cast<int>(D_same.size())  == half);
    assert(static_cast<int>(D_cross.size()) == half);

    std::vector<int> offsets(half);
    for (int r = 0; r < half; r++) offsets[r] = r;

    auto res_same = std::make_unique<TRLWE<P>>();
    auto res_cross = std::make_unique<TRLWE<P>>();
    auto ct_conj = std::make_unique<TRLWE<P>>();

    LinearTransformBSGS<P>(*res_same, ct, D_same, offsets, k_step, gk);

    ConjugateSlots<P>(*ct_conj, ct, gk);
    LinearTransformBSGS<P>(*res_cross, *ct_conj, D_cross, offsets, k_step, gk);
    ct_conj.reset();

    for (int k = 0; k <= static_cast<int>(P::k); k++)
        for (uint32_t i = 0; i < P::n; i++)
            res[k][i] = (*res_same)[k][i] + (*res_cross)[k][i];
}

// ---------------------------------------------------------------------------
// SlotToCoeff: applies S^{-1} (the SlotDecode matrix) in slot domain.
//
// Given a TRLWE ct whose user slots are (q_0, ..., q_{n-1}), produces a
// TRLWE whose polynomial has coefficients (q_0, ..., q_{n-1}) directly —
// i.e. the inverse of CoeffToSlot. The implementation is identical to
// CoeffToSlot; only the precomputed diagonals differ (built from S^{-1}
// via BuildSlotToCoeffKeys instead of BuildCoeffToSlotKeys).
// ---------------------------------------------------------------------------

template <class P>
void SlotToCoeff(TRLWE<P> &res, const TRLWE<P> &ct,
                 const std::vector<std::array<uint64_t, P::n>> &D_same,
                 const std::vector<std::array<uint64_t, P::n>> &D_cross,
                 int k_step,
                 const GaloisKey<P> &gk)
{
    CoeffToSlot<P>(res, ct, D_same, D_cross, k_step, gk);
}

// ---------------------------------------------------------------------------
// FFT-style CRT-factor transforms used by BFV bootstrapping.
//
// These implement the sparse power-of-two CRT factors from
// Bootstrapping_BGV_BFV.  The product is the SlotDecode/SlotEncode map up to a
// bit-reversal of coefficient positions inside each Galois row.  BFV
// bootstrapping uses the forward and inverse product as a matched pair around
// NoisyDecrypt, so that permutation cancels without paying for an additional
// sparse permutation stage.
// ---------------------------------------------------------------------------

template <class P>
void SlotToCoeffCRT(TRLWE<P> &res, const TRLWE<P> &ct, const GaloisKey<P> &gk)
{
    detail::CRTSlotToCoeffProduct<P>(res, ct, gk);
}

template <class P>
void CoeffToSlotCRT(TRLWE<P> &res, const TRLWE<P> &ct, const GaloisKey<P> &gk)
{
    detail::CRTCoeffToSlotProduct<P>(res, ct, gk);
}

}  // namespace c2s
}  // namespace TFHEpp
