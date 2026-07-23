# Gentry–Lee Matrix FHE

TFHEpp includes an experimental implementation of the complex Gentry–Lee
(GL) scheme in `include/gl/gl.hpp`, plus the low-depth SHIP bootstrap from
*Low-Depth Bootstrapping for Matrix-Native FHE* in
`include/gl/gl-bootstrap.hpp`. It follows the matrix encoding and
relinearization from *Fully Homomorphic Encryption for Matrix Arithmetic* and
uses TFHEpp's CKKS/BFV Double Decomposition (DD) arithmetic for products, key
switching, masked selection, and HMux.

## Representation

The implementation uses

```text
R' = Z_q[I,X,Y,W] / (I^2+1, X^n-I, Y^n-I, Phi_p(W))
R  = Z_q[I,X,W]   / (I^2+1, X^n-I, Phi_p(W)).
```

The ternary secret is in `R`. A ciphertext in `R'` is stored as `n`
independent RLWE slices, one per `Y` coefficient, and decrypts with the paper's
`body + mask * secret` convention. `GLMatrixBatch<GLP>` contains `phi(p)`
complex `n` by `n` matrices.

`GLParameter<BaseP, N, P, AuxiliaryLogQ>` checks the implementation's current
constraints:

- `N` is a power of two;
- `P` is an odd prime;
- `BaseP::n == 2 * N * (P - 1)` and `BaseP::k == 1`;
- the torus is `__uint128_t` or TFHEpp's multi-limb type; and
- every switch satisfies `ciphertext LogQ + AuxiliaryLogQ <= torus width`.

The paper-profile storage aliases are:

| Alias | Slice degree | Torus storage | Paper `log Q` | Paper `log P` | Paper `log(PQ)` |
|---|---:|---:|---:|---:|---:|
| `GL256p17Parameter` | 8192 | 256 bits | 180 | 34 | 214 |
| `GL512p17Parameter` | 16384 | 448 bits | 338 | 92 | 430 |
| `GL1024p17Parameter` | 32768 | 896 bits | 641 | 220 | 861 |

The former degree-8192 alias used `lvl4param` and only 128 bits of storage;
that is insufficient for the paper's 214-bit `P*Q` requirement. The new
four-limb degree-8192 base fixes the representation limit. The compile-time
`GLSHIPPaperParameterProfile` specializations record Table 1 and statically
check that each profile fits its torus storage.

## Double Decomposition

DD is used in both places where a product would otherwise exceed the torus
word:

1. Each operand of a matrix or Hadamard product is decomposed at its active
   modulus in base `Bbar`. Small digit products are accumulated in a 384-bit
   accumulator for 128-bit parameters, or a wide signed limb accumulator for
   multi-limb parameters, and then rounded/rescaled.
2. A switch input is decomposed in the primary base. At key generation, each
   unseeded evaluation-key RLWE row is independently decomposed in `Bbar` and
   stored as signed native digits, retaining both RLWE components. Evaluation
   expands one stored row at a time, multiplies it by the primary digit, and
   recombines the result; it does not decompose an evaluation-key ciphertext
   online.

Evaluation keys are generated at `q * q0`. Their plaintext is
`q0 * gadget * source_secret`; after DD recomposition, switching performs a
signed rounded division by `q0`. This is the auxiliary-modulus step from the
GL construction and prevents evaluation-key error from being amplified by an
unscaled same-modulus switch.

The SHIP path preserves the paper's hybrid-RNS operation boundaries while
using DD for the underlying decomposition:

- masked-column ciphertext/plaintext products are accumulated modulo `P*Q`
  and Equation (15)'s single ModDown removes `P` after the sum;
- each HMux stage accumulates the raw DD body and mask switches for every
  radix branch modulo `P*Q`, then performs one ModDown for the stage; and
- a product-tree node multiplies and relinearizes modulo its input `Q`, then
  rescales the resulting two-component ciphertext once.

The last point avoids the former component-wise tensor floor
`r00 + (r01+r10)s + r11*s^2`; relinearize-before-rescale has the standard
two-component floor `r0+r1*s`. Product relinearization keys are consequently
generated at each tree node's input modulus, not its post-rescale modulus.

Fresh ciphertexts and evaluation keys default to coefficient-domain Gaussian
noise with standard deviation 3.2, matching the paper's reference choice.
`GLNoiseAtLevel<LogQ>(sigma)` converts that absolute standard deviation into
the modular value expected by `CKKSNoise`.

## Exact NTT Multiplication

The `p=17` base-ring multiplier uses the scheme-neutral arithmetic in
`include/modular_ntt.hpp`.  That header provides radix-2 cyclic and negacyclic
NTTs, Rader prime-length NTTs, prime-cyclotomic evaluation/interpolation, and
centered two-prime CRT reconstruction.  The GL adapter identifies `I=X^n`,
applies a length-`2n` negacyclic transform on the combined `(I,X)` axis, and
uses Rader's length-17 transform on the `W` axis.

Two 62-bit transform primes reconstruct the production `n512p17` 85-by-16-bit
DD product exactly.  A 16-by-16-bit digit product needs one prime.  For the
wide-torus-by-small products used by encryption and key generation, the wide
operand is split into unsigned power-of-two chunks, each chunk is multiplied
exactly with two primes, and the results are recombined modulo the native
power-of-two torus.  Products that exceed these proved bounds retain the
coefficient-domain reference fallback.  NTT/CRT is therefore an exact
arithmetic backend and adds no approximation noise to the estimator.

The stored ciphertexts and packed DD evaluation keys remain in coefficient
form.  Transform primes are temporary multiplication machinery; no RNS limbs
are added to the persistent key representation.

Small-DD switching reuses those temporary transforms within one switch.  Each
packed key row is transformed once, each primary input digit is transformed
once per `Y` slice, and all primary-row products are accumulated before one
inverse transform.  The `n512p17` path uses a 54 MiB transient key-spectrum
cache; it does not enlarge the serialized evaluation key.  Slice-at-a-time
decomposition and fused Bbar recomposition also remove about 25.4 GiB of old
full-polynomial scratch at `n512p17`.  Beyond the caller-owned 896 MiB output
ciphertext, the optimized switch needs roughly 60 MiB of working storage.  The
full switch's per-prime transform counts fall from 221,184 forward and 110,592
inverse transforms to 2,264 forward and 27,648 inverse transforms.  The
regression's dense raw base switch takes about one second on the development
host, and its measured fixed/eight-slice costs project the full 512-slice
switch at roughly two minutes.  This is a component projection, not a measured
end-to-end bootstrap runtime.

## Basic Use

```cpp
#include <tfhe++.hpp>

using GLP = TFHEpp::GL256p17Parameter;
constexpr std::uint32_t input_log_q = 100;
constexpr std::uint32_t log_delta = 40;

using Input = TFHEpp::GLCiphertext<GLP, input_log_q, log_delta>;
using Product = TFHEpp::GLMatrixMultResult<
    GLP, input_log_q, log_delta, input_log_q, log_delta>;

TFHEpp::Key<typename GLP::baseP> secret;
TFHEpp::keyGen<typename GLP::baseP>(secret);

TFHEpp::GLMatrixBatch<GLP> a;
TFHEpp::GLMatrixBatch<GLP> b;
// Fill a(batch, row, column) and b(batch, row, column).

Input encrypted_a;
Input encrypted_b;
TFHEpp::GLEncrypt(encrypted_a, a, secret);
TFHEpp::GLEncrypt(encrypted_b, b, secret);

TFHEpp::GLMatrixRelinKey<GLP, Product::log_q> matrix_key;
TFHEpp::GLMatrixRelinKeyGen(matrix_key, secret);

Product encrypted_product;
TFHEpp::GLMatrixMultiply(encrypted_product, encrypted_a, encrypted_b,
                         matrix_key);

TFHEpp::GLMatrixBatch<GLP> decoded;
TFHEpp::GLDecrypt(decoded, encrypted_product, secret);
```

The native trace product decodes to `A * B^* / n`. The default
`GLMatrixMultiply` and `GLPlaintextMatrixMultiply` cancel the `1/n` factor, so
their result is `A * B^*`. Use `GLMatrixMultiply<false>` to retain the paper's
unscaled trace result. To compute an ordinary `A * B`, encode the adjoint of
`B` as the right operand.

After multiplication, the output active modulus is
`min(lhs LogQ, rhs LogQ) - max(lhs LogDelta, rhs LogDelta)`, and its scale is
`min(lhs LogDelta, rhs LogDelta)`. Generate the relinearization key for that
output modulus, as in the example.

## Operations

The public API currently provides:

| Operation | API | Evaluation key |
|---|---|---|
| Encode/decode | `GLEncode`, `GLDecode` | none |
| Encrypt/decrypt | `GLEncrypt`, `GLDecrypt` | none |
| Addition | `GLAdd` | none |
| Ciphertext matrix product | `GLMatrixMultiply` | two big DD switches in `GLMatrixRelinKey` |
| Plaintext/ciphertext matrix product | `GLPlaintextMatrixMultiply` | none |
| Hadamard product | `GLHadamardMultiply` | one small DD switch |
| Complex conjugation | `GLConjugate` | one small DD switch |
| Row rotation | `GLRotateRows` | one small DD switch, bound to the rotation key |
| Column rotation | `GLRotateColumns` | none |
| Batch (`W`) rotation | `GLRotateBatches` | one small DD switch, bound to the rotation key |
| Transpose | `GLTranspose` | one big DD switch |
| Conjugate transpose | `GLConjugateTranspose` | one big DD switch |
| Grouped slots-to-coefficients | `GLSHIPSlotsToCoefficients` | one big conjugate-transpose key plus BSGS W-rotation keys |
| SHIP half-bootstrap | `GLSHIPHalfBootstrap` | dense-to-sparse, encrypted-mask, HMux, relin, and conjugation keys |
| Full SHIP bootstrap | `GLSHIPBootstrap` | StC and half-bootstrap keys |

The conjugate-transpose switch source is the same
`s(-I,Y^-1,W^-1)` source used by `GLMatrixRelinKey::conjugate_key`, so that
member can be reused at the same active modulus.

## Low-Depth SHIP Bootstrap

`GLSHIPBootstrapSchedule` describes a complete reference schedule:

- a natural-order grouped GL StC, with matrix-native X inversion and a
  diagonal BSGS W inversion followed by one delayed rescale;
- DD dense-to-sparse switching at the level-zero modulus;
- two Gaussian channels implementing Equation (1) of the bootstrap paper;
- encrypted row-wise masked-column selectors at modulus `P*Q`;
- radix-B, X-only hidden HMux stages implemented as DD external key switches;
- a balanced `h+1` product tree with one relinearization/rescale per layer;
- conjugate-and-add sine extraction, multiplication of the second channel by
  formal `I`, and direct Y-slice reassembly.

The main interfaces are:

```cpp
TFHEpp::GLSHIPBootstrapKeyGen(key, dense_secret, sparse_secret, intervals);
TFHEpp::GLSHIPSlotsToCoefficients<Schedule>(coeff, input, key.stc_key);
TFHEpp::GLSHIPHalfBootstrap<Schedule>(output, coeff, key);
TFHEpp::GLSHIPBootstrap<Schedule>(output, input, key);
```

`sparse_secret` must have the schedule's Hamming weight and coefficients in
`{0,+1,-1}`. A coefficient's flat W-major index determines its `X`, `W`, and
Gaussian `I` component. Each public `GLSHIPSupportInterval` must contain the
corresponding nonzero coefficient. As in SHIP, those intervals leak bounded
support information; the exact position, sign, Gaussian phase, and HMux
digits remain encrypted.

The regression schedule uses `n=2`, `p=5`, sparse weight 3, and nonzero
Gaussian noise. It is deliberately insecure and exists to exercise every
algebraic branch quickly. Production callers must select a schedule consistent
with the chosen Table-1 profile and perform a fresh estimator/noise analysis.

The companion estimator in `../../Parameter-Selection/python/GLnoise.py` keeps a
`legacy` comparison mode for the old operation boundaries. The fused-DD model
first reaches the paper's measured precision for `n512p17` and `n1024p17` at
47- and 50-bit uniform tree scales without increasing `Q` or `P`. Because the
47-bit n512 result has only about 0.09 bit of modeled margin, the selected
n512 alias uses 49 bits, four 85-bit primary rows (covering its 338-bit `Q`),
and pre-decomposed signed 16-bit evaluation-key rows. The reconstructed
`n256p17` schedule remains below target at its largest feasible 26-bit tree
scale, so that profile is not yet certified by the available schedule data.

The four-row DD layout is a TFHEpp storage/arithmetic choice that mirrors the
paper's `dnum=4` partition count; the paper's `dnum` denotes RNS partitions and
does not prescribe an 85-bit gadget base. Likewise, the packed unseeded key
format is an implementation representation, not a seeded or online
evaluation-key decomposition technique claimed by either paper.

The two estimator-screened reference aliases are:

| Alias | q0 | X+W StC | Half-bootstrap Q | Tree scale | Output Q | Estimated precision |
|---|---:|---:|---:|---:|---:|---:|
| `GLSHIP512p17FusedDDSchedule` | 48 | 18+19 | 338 | 49 | 93 | 16.27 bits |
| `GLSHIP1024p17FusedDDSchedule` | 50 | 19+20 | 641 | 50 | 391 | 16.43 bits |

Using the paper's aggregate count of 1,504 masked columns, the n512 alias has a
compile-time coefficient-payload estimate of 8,458,338,304 bytes (7.88 GiB).
This is the unseeded representation used by the implementation: the figure
includes both components of every DD key row and native-width masked-column
ciphertexts. `GLSHIPBootstrapKeyPackedPayloadBytes(key)` reports the exact
coefficient payload after concrete support intervals are known. Archive and
allocator metadata add a small overhead; saving does not require a second
in-memory key, although the atomic writer temporarily needs room for the
destination file.

## Current Boundary

This remains a correctness-oriented implementation. Its bootstrap follows the
paper's algebra, and base-ring multiplication now has an exact Rader/NTT
backend.  Small-DD switches now reuse transforms across primary rows and `Y`
slices, but big polynomial switches do not yet have a corresponding full-ring
transform, and the implementation does not fuse Rader butterflies as
aggressively as the paper or provide its optimized modular matrix-multiplication
kernels.  Consequently, the component benchmarks must not be interpreted as a
full-bootstrap runtime.

Storage sufficiency is also not a security or precision proof. The paper's
`n256p17` and `n512p17` values use the full 214- and 430-bit 128-bit-security
limits, leaving no modulus margin; `n1024p17` uses 861 of 868 permitted bits.
TFHEpp's power-of-two, multi-limb backend has enough bits to represent those
profiles. The fused-DD estimator validates the modulus budget for the latter
two profiles under its average-case model, but full-size bootstrap noise has
not yet been measured empirically. The selected `n512p17` tree scale estimates
16.27 bits of precision, 1.33 bits above the paper's reported 14.94 bits.
Do not label these aliases production-secure until that validation and
constant-time optimized kernels are complete.

Build and run the regression with:

```bash
cmake -S . -B build -DENABLE_TEST=ON
cmake --build build --target modular_ntt gl_ntt gl_scheme gl_bootstrap
./build/test/common/modular_ntt
./build/test/gl/gl_ntt
./build/test/gl/gl_scheme
./build/test/gl/gl_bootstrap
```

The regression uses `n=2`, `p=5` and nonzero Gaussian noise. It covers
encoding, the paper trace identity, encryption, both DD switch sizes, matrix
and Hadamard products, all listed automorphisms, and auxiliary-modulus
switching.

The bootstrap regression additionally covers slice encoding, natural-order
grouped StC, dense-to-sparse switching, masked candidates, nontrivial X HMux,
both Gaussian channels, the balanced product tree, half-bootstrap, and the
full bootstrap wrapper.
