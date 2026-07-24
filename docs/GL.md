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

One, two, or three 62-bit transform primes are selected from a proved
coefficient bound. Two primes reconstruct the production `n512p17`
85-by-16-bit DD product exactly; the full-polynomial big switch uses the third
prime where its accumulation bound requires it. A 16-by-16-bit digit product
needs one prime. For wide-torus-by-small products used by encryption and key
generation, the wide operand is split into unsigned power-of-two chunks, each
chunk is multiplied exactly, and the results are recombined modulo the native
power-of-two torus. Products outside these proved bounds retain the
coefficient-domain reference fallback. NTT/CRT is therefore an exact
arithmetic backend and adds no approximation noise to the estimator.

Wide-by-wide multiplication includes every chunk pair accumulated on an
output diagonal in its CRT bound. It uses two primes and chooses the largest
safe chunk width for the active modulus. At n512/logQ 338 this is seven chunks
per operand instead of the former sixteen one-prime chunks; diagonals whose
shift is already zero modulo `2^LogQ` are not evaluated.

Ciphertext products share those chunk transforms across all four input
components. The three tensor outputs `(a0*b0, a0*b1+a1*b0, a1*b1)` are formed
directly, with the cross term accumulated before inversion and its extra
factor included in the CRT bound. Across the five n512 product-tree levels,
128-call batches are 24--31% faster than four separate base products.

The common radix-2 implementation precomputes stage twiddles and bit-reversal
swaps, and omits the modular product for every unity twiddle. Its
pseudo-Mersenne reduction exploits `p = 2^62-c`, avoiding the compiler's
128-by-64 division helper. Inverse negacyclic and Rader transforms fold their
normalization constants into already-required coefficient or kernel
multipliers. These optimizations are in the scheme-neutral `modular_ntt.hpp`,
not in a GL-specific FFT wrapper. The fixed 16-point convolution inside the
`p=17` Rader transform additionally uses flat tables and compile-time-unrolled
stages instead of the dynamic radix-2 plan. Fixed twiddles, negacyclic twists,
Rader kernels, and cyclotomic interpolation weights use precomputed Shoup
quotients, replacing pseudo-Mersenne reduction in constant products. Larger
radix-2 stages retain residues in `[0,2p)` and canonicalize once at the
transform boundary; the 16-point Rader convolution remains canonical because
lazy boundary work is slower at that size. On the development host, these
changes reduce a single-thread n512 full-ring forward/inverse pair from about
0.87/0.92 seconds to about 0.40/0.40 seconds.

Configuring with `-DUSE_HEXL=ON` selects Intel HEXL for the large power-of-two
transforms while retaining the exact fixed Rader implementation for `p=17`:

```sh
cmake -S . -B build-hexl -DENABLE_TEST=ON -DUSE_HEXL=ON
cmake --build build-hexl -j
```

The custom transform root is the same root used by the native plan. GL keeps
HEXL's bit-reversed spectra through pointwise operations and translates the X
automorphism maps accordingly; the public common-arithmetic transform still
returns natural order. A cyclic GL Y transform is obtained exactly from
HEXL's negacyclic transform by the standard `psi^-i`/`psi^i` coefficient
twists. Direct spectrum-evaluation tests cover both cyclic and negacyclic
ordering, in addition to round-trip and convolution tests. The native path
remains the default and has no HEXL dependency.

The stored ciphertexts and packed DD evaluation keys remain in coefficient
form.  Transform primes are temporary multiplication machinery; no RNS limbs
are added to the persistent key representation.

Small-DD switching reuses those temporary transforms across calls. Each
packed key row is transformed once, each primary input digit is transformed
once per `Y` slice, and all primary-row products are accumulated before one
inverse transform.  The `n512p17` path uses a 54 MiB transient key-spectrum
cache per 338-bit switch; it does not enlarge the serialized evaluation key.
HMux also sums all eight body/mask branch switches in the NTT domain. The
native path uses batched 128-bit pointwise accumulators; the optional HEXL
path uses AVX-512 vector products and modular adds. Slice-at-a-time
decomposition and fused Bbar recomposition also remove about 25.4 GiB of old
full-polynomial scratch at `n512p17`.  Beyond the caller-owned 896 MiB output
ciphertext, the optimized switch needs roughly 60 MiB of working storage.  The
full switch's per-prime transform counts fall from 221,184 forward and 110,592
inverse transforms to 2,264 forward and 27,648 inverse transforms.  The
regression's current dense raw base-switch measurement is about 0.05 second,
or roughly 27--30 seconds when multiplied across 512 slices. This is a
component projection, not a measured full switch.

Within each radix-4 HMux stage, symmetric DD decomposition is hoisted across
the four X automorphisms. The two source components are decomposed and
transformed once; exact signed coefficient permutations become NTT-spectrum
permutations. This reduces the stage's input forward transforms from 64 to 16
without changing its two-prime CRT bound. On the development host, 128 warm
complete stages take about 2.1 seconds instead of 2.34 seconds with unhoisted
transforms on the native path. With HEXL, 128 complete stages take about
1.8--1.9 seconds with 32 SMT workers; 16 physical-core workers sustain about
89--90 stages per second for this bandwidth-bound kernel.

Masked-column candidate construction likewise hoists the expensive torus
phase roots out of the candidate loop. Roots for the two Gaussian channels
and the extended `W` basis are evaluated once per input mask; candidate
differences, signs, and Gaussian phases are then exact complex products and
conjugations. Candidates with the same `W` index and Gaussian phase are
encoded once at fine-X zero; the remaining fine-X shifts use exact base-ring
automorphisms. Those shifts are applied directly as NTT-spectrum permutations,
so each group also needs only one pair of forward transforms. All three
rewrites remain bit-for-bit equal to direct candidate encoding and
multiplication. On the development host, 128 warm 48-candidate columns take
about 1.47 seconds on the native path, down from 1.88 seconds before spectrum
reuse, 2.12 seconds with phase-root hoisting alone, and about 2.50 seconds
before these optimizations. HEXL vector products reduce the same batch to
about 1.23 seconds.

HMux X automorphisms act directly on the `(I,X,W)` base polynomial. The former
generic route lifted one base polynomial into two 448 MiB full GL temporaries
per automorphism at n512; the direct coefficient permutation is checked
exactly against that generic map and avoids this allocation entirely.

The big DD switch similarly caches full-ring key spectra. At n512 the cache is
about 4.5 GiB; StC prepares it once, reuses it for both conjugate transposes,
and releases it before the W transform. Representative 32-thread measurements
on the development host are 3.6 seconds to prepare that cache, 6.2 seconds per
warm big switch, and 9.3 seconds cold. The specialized X trace avoids a
448 MiB plaintext and takes about 1.1 seconds. Four W diagonals are accumulated
before inverse transforms, making the complete W stage about 1.5 seconds.
With HEXL's cyclic and negacyclic transforms, representative cache-prepare,
warm-switch, and cold-switch measurements are 1.08, 3.37, and 4.67 seconds,
respectively. The complete n512 StC path measures 30.7 seconds with 32 workers.

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

The default keeps MaskedColumn and HMux in one OpenMP loop. On an SMT host
where HMux saturates memory bandwidth before MaskedColumn, callers can split
bounded factor tiles and give HMux a smaller team:

```cpp
TFHEpp::GLSHIPBootstrapExecutionOptions execution{
    .hmux_threads = 16,
    .factor_tile_size = 256,
    .batch_hmux_products = true,
    .block_hmux_key_spectra = true,
};
TFHEpp::GLSHIPBootstrap<Schedule>(output, input, key, execution);
```

The worker count is intentionally explicit and should be benchmarked per
machine. A zero `hmux_threads` retains the portable default. At n512 a
256-factor staging tile adds 448 MiB of ciphertext storage while it is live.
With `USE_HEXL`, `batch_hmux_products` uses the common AVX-512DQ
pseudo-Mersenne MAC to reduce four exact pointwise products together; without
AVX-512DQ the regular HEXL multiply/add path is retained. With the fast MAC,
`block_hmux_key_spectra` stores each transient HMux key spectrum in
4096-coefficient blocks with its four primary rows adjacent. It does not
change the serialized key size or the 54 MiB transient-cache size per n512
switch; only one cache layout is materialized by the production path.

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

NTT spectra and workspaces are derived state and are omitted from archives.
Retaining every n512 HMux cache simultaneously would consume about 39.2 GiB,
and retaining all masked-column spectra would add roughly another 6.6 GiB.
`GLSHIPHalfBootstrap` therefore evaluates sparse terms outside the slice loop:
it prepares one term's approximately 1.3 GiB of HMux spectra and roughly
0.2 GiB masked-column cache, reuses them for all 1,024 Y/channel slices, then
releases both. Its online binary product tree preserves exactly the original
balanced factor pairing. A batch of 1,024 base ciphertexts is 1.75 GiB; the
largest carry transient is about seven such batches (12.25 GiB), and it does
not overlap the per-term NTT caches. Worker-local NTT scratch is reused within
each term and then freed.

### Performance status

The native production n512 StC path has been measured end-to-end at 47.1
seconds with 32 OpenMP threads. A combined sequential component run peaked at
about 17.1 GiB RSS. Native warm throughput for the dominant half-bootstrap
kernels is about 65 complete HMux stages per second and 94 48-candidate masked
columns per second on the same Ryzen 9 7950X3D host. The n512 schedule needs
95,232 HMux stages and 31,744 masked columns, projecting about 24.4 and 5.6
minutes. Native level-by-level throughput projects the complete 31,744-node
product tree at about 3.9 minutes. Including StC and the smaller
dense-to-sparse, cache-preparation, and final-conjugation phases, roughly 35--36
minutes remains the native component projection.

With `USE_HEXL=ON` and 32 workers, measured component rates are about 68--70
complete HMux stages and 104 masked columns per second; the product tree
projects to about 165--168 seconds and StC measures 30.7 seconds. This gives a
roughly 32--33 minute component projection. HMux is memory-bandwidth-bound on
this 16-core/32-thread CPU: running with
`OMP_NUM_THREADS=16 OMP_PROC_BIND=spread OMP_PLACES=cores` raises HMux
throughput to about 89--90 stages per second, while masked-column throughput
is about 98 per second, the product tree projects to about 191 seconds, and
StC measures 32.7 seconds. That host-specific profile gives a roughly 27--28
minute component projection. Other CPUs should benchmark both physical-core
and SMT worker counts; no topology-dependent limit is hard-coded.

With `OMP_NUM_THREADS=32 OMP_PROC_BIND=spread OMP_PLACES=cores`, a
production-shape batch of 256 MaskedColumn plus one-stage HMux factors took
4.71--4.81 seconds with the execution options above. In paired runs, the same
staged path with row-major key spectra took 5.06--5.31 seconds, and the fused
32-worker loop took 6.49--6.51 seconds. Thus blocking saved 4.9--11.3% over
the otherwise identical staged path and the complete tuning saved
25.9--27.7% over the fused loop. The benefit was not consistent for smaller
tiles, so 256 remains the recommended n512 tile. The common 16K, 32-product
modular MAC microbenchmark took 0.125--0.140 seconds versus 0.162--0.184
seconds for separate HEXL multiply/add calls (23--24% faster). The factor
benchmark uses a synthetic one-stage key; the separately measured three-stage
rates still suggest about a 22% reduction in the complete sparse-factor phase,
not a measured full-bootstrap runtime. Run the tuning benchmark with
`TFHEPP_GL_N512_MASKED_BENCH=1` and `TFHEPP_GL_N512_FACTOR_BENCH=1`; run the
arithmetic microbenchmark with `TFHEPP_MODULAR_MAC_BENCH=1`.

These are component projections, not measured end-to-end runtimes. The full
7.88 GiB production key has deliberately not been generated merely to obtain
a timing number.

## Current Boundary

This remains a correctness-oriented implementation. Its bootstrap follows the
paper's algebra, and base-ring and full-polynomial DD operations now use exact
Rader/NTT backends with reusable spectra. The specialized X trace, grouped W
transform, big transpose switches, masked columns, fused HMux branches, and
online product tree all compile for the n512 schedule. The toy regression runs
the complete wrapper, but a full n512 bootstrap with a real 7.88 GiB key has
not yet been run. Component benchmarks must therefore not be reported as a
measured full-bootstrap runtime.

The remaining performance gap is chiefly the 95,232 exact HMux stages and
31,744 masked columns required by the algorithm, followed by the product tree.
The optional HEXL backend supplies SIMD power-of-two NTTs and pointwise
products, but the code does not yet provide the paper implementation's full
set of fused kernels, NUMA-aware key streaming, or constant-time production
hardening.

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
