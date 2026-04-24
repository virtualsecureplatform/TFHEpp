#pragma once

#include <array>
#include <cstdint>
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

    TRLWE<P> res_same, res_cross, ct_conj;

    LinearTransformBSGS<P>(res_same, ct, D_same, offsets, k_step, gk);

    ConjugateSlots<P>(ct_conj, ct, gk);
    LinearTransformBSGS<P>(res_cross, ct_conj, D_cross, offsets, k_step, gk);

    for (int k = 0; k <= static_cast<int>(P::k); k++)
        for (uint32_t i = 0; i < P::n; i++)
            res[k][i] = res_same[k][i] + res_cross[k][i];
}

}  // namespace c2s
}  // namespace TFHEpp
