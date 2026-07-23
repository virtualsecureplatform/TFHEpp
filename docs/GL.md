# Gentry–Lee Matrix FHE

TFHEpp includes an experimental, leveled implementation of the complex
Gentry–Lee (GL) scheme in `include/gl/gl.hpp`. It follows the matrix encoding
and relinearization from *Fully Homomorphic Encryption for Matrix Arithmetic*
and uses TFHEpp's CKKS/BFV Double Decomposition (DD) arithmetic for products
and key switching.

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

`GL256p17Parameter` supplies the paper's `n=256`, `p=17`, degree-8192 ring
shape on `lvl4param`, with a 32-bit auxiliary switch modulus. It is a baseline
type, not a claim that the existing `lvl4param` modulus/noise schedule matches
every parameter in the paper.

## Double Decomposition

DD is used in both places where a product would otherwise exceed the torus
word:

1. Each operand of a matrix or Hadamard product is decomposed at its active
   modulus in base `Bbar`. Small digit products are accumulated in a 384-bit
   accumulator for 128-bit parameters, or a wide signed limb accumulator for
   multi-limb parameters, and then rounded/rescaled.
2. A switch input is decomposed in the primary base. Each evaluation-key RLWE
   row is independently decomposed in `Bbar`, multiplied by the primary digit,
   and recombined.

Evaluation keys are generated at `q * q0`. Their plaintext is
`q0 * gadget * source_secret`; after DD recomposition, switching performs a
signed rounded division by `q0`. This is the auxiliary-modulus step from the
GL construction and prevents evaluation-key error from being amplified by an
unscaled same-modulus switch.

Fresh ciphertexts and evaluation keys default to coefficient-domain Gaussian
noise with standard deviation 3.2, matching the paper's reference choice.
`GLNoiseAtLevel<LogQ>(sigma)` converts that absolute standard deviation into
the modular value expected by `CKKSNoise`.

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

The conjugate-transpose switch source is the same
`s(-I,Y^-1,W^-1)` source used by `GLMatrixRelinKey::conjugate_key`, so that
member can be reused at the same active modulus.

## Current Boundary

This is a correctness-oriented coefficient-domain implementation. The
degree-8192 type is intentionally exposed, but the current multiplication and
trace kernels are not yet replaced by the paper's Rader/NTT and optimized
modular matrix-multiplication path. The alias currently serves compile-time
integration and parameter-shape work; full-size encoding and evaluation are
not practical with these reference kernels.

The SHIP-based low-depth bootstrap from *Low-Depth Bootstrapping for
Matrix-Native FHE* is not implemented here. That construction additionally
requires GL StC/CtS transforms, masked-column candidate generation, and hidden
HMux selection over the non-power-of-two `(X,W)` algebra; the current TFHEpp
bootstrap machinery only supplies the power-of-two CKKS path. The GL code is
therefore leveled and should not be presented as a bootstrapped GL instance.

Build and run the regression with:

```bash
cmake -S . -B build -DENABLE_TEST=ON
cmake --build build --target gl_scheme
./build/test/gl/gl_scheme
```

The regression uses `n=2`, `p=5` and nonzero Gaussian noise. It covers
encoding, the paper trace identity, encryption, both DD switch sizes, matrix
and Hadamard products, all listed automorphisms, and auxiliary-modulus
switching.
