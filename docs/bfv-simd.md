# BFV SIMD Slot Operations

TFHEpp supports BFV-style SIMD (Single Instruction Multiple Data) operations over the 128-bit torus using the `lvl3simdparam` parameter set. This enables packing up to 4096 independent plaintext values into a single ciphertext and operating on all of them in parallel.

## Parameter Set

| Parameter | Value | Description |
|-----------|-------|-------------|
| n | 4096 | Ring dimension (polynomial degree) |
| Q | 2^128 | Ciphertext modulus (implicit, `__uint128_t`) |
| t | 114689 = 7·2^14+1 | Plaintext modulus (prime) |
| Δ | floor(Q/t) ≈ 2^111.19 | Scaling factor (NOT a power of 2) |
| α | 2^{-125} | Noise parameter (σ ≈ 3.2) |
| l, Bgbit | 4, 21 | Standard decomposition |
| l̅, B̅gbit | 8, 16 | Double Decomposition (DD) |

The plaintext modulus t = 114689 is prime with t ≡ 1 (mod 2n = 8192), which causes x^n + 1 to split completely into n linear factors over Z_t. This gives the CRT isomorphism:

    Z_t[x]/(x^n+1) ≅ Z_t × Z_t × ... × Z_t   (n copies)

Each copy is an independent "slot" that can hold a value in Z_t.

## Slot Structure

The 4096 NTT evaluation points split into two **rows** of 2048 under the Galois group (Z/2n)* ≅ Z_2 × Z_{n/2}:

- **Row 0**: slots 0..2047 — indexed by powers of 5 (generator of Z_{n/2})
- **Row 1**: slots 2048..4095 — indexed by negated powers of 5

Slot rotation (via Galois automorphism x → x^5) cyclically shifts slots **within each row independently**. Cross-row operation uses conjugation (x → x^{-1}).

## API

All functions are in namespace `TFHEpp`, templated on `P = lvl3simdparam`.

### Encoding/Decoding

```cpp
// slot values → polynomial coefficients (with Galois permutation)
void SlotEncode<P>(Polynomial<P>& poly, const std::array<uint64_t, P::n>& slots);

// polynomial coefficients → slot values (with inverse Galois permutation)
void SlotDecode<P>(std::array<uint64_t, P::n>& slots, const Polynomial<P>& poly);
```

The encode/decode apply the Galois permutation so that user slot indices align with the Galois group orbit. This makes d=5 automorphism correspond to rotation by 1.

### Encrypt/Decrypt

```cpp
// Encrypt slot vector as TRLWE
void trlweSlotEncrypt<P>(TRLWE<P>& ct,
                         const std::array<uint64_t, P::n>& slots,
                         const Key<P>& key);

// Decrypt TRLWE to slot vector
void trlweSlotDecrypt<P>(std::array<uint64_t, P::n>& slots,
                         const TRLWE<P>& ct,
                         const Key<P>& key);
```

- **Encrypt** uses `floor(m·Q/t)` per coefficient, giving at most 1 unit of rounding error independent of the message value m. This is more precise than `m · floor(Q/t)` which accumulates error proportional to m.
- **Decrypt** uses `round(phase·t/Q)` — the standard BFV decoding that works with the non-power-of-2 Δ.

### Arithmetic

```cpp
// Addition: coefficient-wise TRLWE addition (free, no key needed)
for (int k = 0; k <= P::k; k++)
    for (int i = 0; i < P::n; i++)
        ct_sum[k][i] = ct_a[k][i] + ct_b[k][i];

// Multiplication: Full DD tensor product + relinearization
void TRLWEMultFullDD<P>(TRLWE<P>& res,
                        const TRLWE<P>& a, const TRLWE<P>& b,
                        const relinKeyFFT<P>& relinkeyfft);
```

Multiplication uses the Double Decomposition (DD) approach:
1. Decompose both ciphertexts into l̅=8 base-B̅g digits
2. Transform the digit polynomials once and compute the l̅² pairwise digit products with FFT
3. Accumulate into a 384-bit wide accumulator (`Wide384`) to hold the full 256-bit raw product
4. Divide by Δ using 384-bit / 128-bit long division
5. Relinearize using DD ExternalProduct

### Key Generation

```cpp
// Relinearization key (for multiplication)
auto relinkeyfft = makeRelinKeyFFT<P>(key);

// Galois keys (for rotation)
GaloisKey<P> gk;
GaloisKeyGen<P>(gk, key);
```

### Rotation

```cpp
// Rotate slots by k positions within each row of n/2
void RotateSlots<P>(TRLWE<P>& res, const TRLWE<P>& ct,
                    int steps, const GaloisKey<P>& gk);

// Conjugation (swap rows)
void ConjugateSlots<P>(TRLWE<P>& res, const TRLWE<P>& ct,
                       const GaloisKey<P>& gk);
```

Rotation by k applies automorphism x → x^{5^k} via binary decomposition (at most nbit-1 = 11 sequential EvalAuto calls).

## Implementation Details

### BFV Scaling (Δ = floor(Q/t))

Standard BFV requires Δ = floor(Q/t), NOT a power-of-2 approximation. With power-of-2 Δ, the BFV rescaling wraps at Q/Δ ≠ t, corrupting multiplication when polynomial product coefficients exceed Q/(2Δ). The non-power-of-2 Δ ensures wrapping by multiples of ≈t, which is transparent to the mod-t reduction.

### 384-bit Accumulator for Multiplication

The raw product of two 128-bit torus polynomials can reach 256 bits. The FullDD multiplication accumulates digit products shifted by positional weights (up to 2^224) into a 384-bit `Wide384` accumulator, then divides by Δ (≈2^111) to get the 128-bit torus result. The accumulator handles signed digit products correctly via separate add/subtract paths.

### FFT Fixes for `__uint128_t`

Two critical fixes enable the 64-bit FFT (spqlios) to work correctly for 128-bit torus DD operations:

1. **Sign extension**: TwistFFT converts 64-bit FFT output to 128-bit via `int64_t → __int128_t → __uint128_t`. Without this, negative values (e.g., -5 stored as uint64 = 2^64-5) get zero-extended to 2^64-5 instead of sign-extended to 2^128-5, corrupting any subsequent shift by < 64 bits in RecombineTRLWEFromDD.

2. **Rounding conversion**: TwistFFT uses `execute_direct_torus64_rescale(D=1.0)` which rounds to nearest integer (via `std::llround`), instead of `f64_to_i64` which truncates toward zero. The ±0.5 truncation error per DD level gets amplified by positional shifts (up to 2^112), causing catastrophic error in ExternalProduct.

### Galois Permutation

The `SlotEncode`/`SlotDecode` functions apply a permutation between user-facing slot indices and NTT evaluation points. User slot i corresponds to NTT position `(5^i - 1)/2 mod n` (row 0) or `(2n - 5^i - 1)/2 mod n` (row 1). This makes the d=5 Galois automorphism correspond to cyclic rotation by 1 within each row.

## Example

```cpp
#include <tfhe++.hpp>

using P = TFHEpp::lvl3simdparam;

// Key generation
TFHEpp::Key<P> key;
// ... initialize key with random ±1 values ...
auto relinkeyfft = TFHEpp::makeRelinKeyFFT<P>(key);

// Encrypt two slot vectors
std::array<uint64_t, P::n> slots_a, slots_b;
// ... fill with values in [0, 114689) ...
TFHEpp::TRLWE<P> ct_a, ct_b;
TFHEpp::trlweSlotEncrypt<P>(ct_a, slots_a, key);
TFHEpp::trlweSlotEncrypt<P>(ct_b, slots_b, key);

// Multiply (slot-wise)
TFHEpp::TRLWE<P> ct_mul;
TFHEpp::TRLWEMultFullDD<P>(ct_mul, ct_a, ct_b, *relinkeyfft);

// Decrypt
std::array<uint64_t, P::n> result;
TFHEpp::trlweSlotDecrypt<P>(result, ct_mul, key);
// result[i] == (slots_a[i] * slots_b[i]) % 114689
```

## Files

| File | Contents |
|------|----------|
| `include/bfv/bfv-slots.hpp` | SlotEncode/Decode, encrypt/decrypt, GaloisKey, rotation |
| `include/bfv/bfv++.hpp` | TRLWEMultFullDD, Wide384, relinearization |
| `include/params/128bit.hpp` | `lvl3simdparam` definition |
| `include/mulfft.hpp` | TwistFFT/TwistIFFT with sign-extension + rounding |
| `test/bfv/simd_encode.cpp` | Encode/decode round-trip test |
| `test/bfv/simd_ops.cpp` | Addition + multiplication test |
| `test/bfv/simd_rotation.cpp` | Slot rotation test |

## Limitations

- **Multiplication depth**: One level only (the noise budget after one multiplication is nearly exhausted). Bootstrapping or modulus switching would be needed for deeper circuits.
- **Performance**: The FullDD multiplication uses FFT-backed digit products when the digit-product bound fits safely below the `double` mantissa. The 384-bit accumulation and Δ-division are still scalar and remain the next major optimization target.
- **AVX512**: The current build requires `USE_AVX512=OFF` due to a missing s=1 inverse twiddle table entry in the AVX512 Stockham FFT path. The AVX2 path works correctly.
