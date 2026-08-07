# CKKS

TFHEpp includes experimental CKKS support in `include/ckks/ckks.hpp`. The
current target path is a leveled CKKS multiplication followed by dense CKKS
bootstrapping. The current storage-practical full-size correctness target is
seeded hybrid giant streamed bootstrapping with compact seeded
double-decomposition (DD) EvalMod and product relinearization. This 2^15-ring
target is a 128-bit-secure parameter set under the local estimator when keys
are generated with the current lvl6 noise parameter.

The most practical full-size configuration is the tuned lvl6 seeded hybrid
giant streamed bootstrap path:

- schedule: `TFHEpp::lvl6CKKSDenseBootstrapTunedSchedule`
- ring degree: `n = 32768`
- torus storage: `896` bits
- boot modulus: `boot_log_q = 896`
- input level: `input_log_q = 58`
- output level: `output_log_q = 116`
- scale bits: `log_delta = 52`
- lvl6 CKKS noise: `α = 2^-872`
- low-level CKKS noise floor: `σ >= 3.2`
- bounded bootstrap secret weight: `16`
- DD relin bases: `Bg = 2^4`, `Bbar = 2^16`

CKKS APIs and schedules are compile-time types. Changing a schedule, level,
scale, sparse key weight, CKKS noise parameter, or key layout requires
regenerating the relevant bootstrap and product evaluation keys.

## Core Types

| Type | Purpose |
|------|---------|
| `CKKSSlotVector<P>` | Plain complex slot vector for tests and encoding helpers. |
| `CKKSCiphertext<P, LogQ, LogDelta>` | Two-component CKKS ciphertext at compile-time level and scale. |
| `CKKSMultResult<...>` | Three-component multiplication result before relinearization. |
| `CKKSRelinKey<P, LogQ>` | Relinearization key for CKKS multiplication at level `LogQ`. |
| `CKKSDDRelinKey<P, LogQ>` | Double-decomposition CKKS relinearization key at level `LogQ`. |
| `CKKSDenseBootstrapSchedule<...>` | Compile-time dense bootstrap parameter schedule. |
| `Schedule::InputCiphertext` | Ciphertext type accepted by a dense bootstrap schedule. |
| `Schedule::OutputCiphertext` | Ciphertext type produced by a dense bootstrap schedule. |

Useful aliases are defined near the end of `include/ckks/ckks.hpp`, including:

- `lvl6CKKSDenseBootstrapInput`
- `lvl6CKKSDenseBootstrapOutput`
- `lvl6CKKSDenseBootstrapTunedSeededHybridGiantStreamedFilesystemKeyProvider`
- `lvl6CKKSDenseBootstrapTunedSeededHybridGiantStreamedDDEvalModFilesystemKeyProvider`

## Dense Bootstrapping Pipeline

Dense CKKS bootstrapping follows this structure:

1. Normalize the input ciphertext if it comes from a multiplication.
2. Optionally switch from an external/product key to the bootstrap key.
3. Raise the ciphertext level to `Schedule::boot_log_q`.
4. Apply coefficient-to-slot linear transforms.
5. Split real and imaginary components.
6. Run bounded cosine EvalMod, including inverse and double-angle stages when
   configured by the schedule.
7. Apply slot-to-coefficient linear transforms.
8. Optionally switch back from the bootstrap key to the external/product key.

The practical product path uses prebuilt product evaluation keys. Runtime
bootstrap/product calls should load these keys from disk rather than generating
them on demand.

## Functional EvalLUT (IACR 2024/1623)

`include/ckks/functional_bootstrapping.hpp` implements the first-order
trigonometric-Hermite LUT construction from Theorem 1 and Corollary 1 of
*General Functional Bootstrapping using CKKS*. Given a table
`f(0), ..., f(p-1)`, `CKKSBuildFunctionalBootstrapLUT` constructs the complex
coefficients of Equation (4). Its real part has both required properties:

- `R(k/p) = f(k)` for every plaintext representative;
- `R'(k/p) = 0`, which gives quadratic input-noise reduction.

The reference evaluator is periodic, so modulus-raising overflow terms that are
integer multiples of its period do not change the result.

TFHEpp provides two encrypted integrations. The direct prototype uses a real
Chebyshev representation of `R` on a caller-selected bounded interval:

```cpp
auto lut = TFHEpp::CKKSBuildFunctionalBootstrapLUT({0.0, 1.0});
auto polynomial =
    TFHEpp::CKKSBuildFunctionalBootstrapChebyshevPolynomial(lut, 24);
```

`CKKSFunctionalBootstrapEvalLUTNormalized` expects encrypted inputs normalized
to `[-1, 1]`. For a polynomial built with `input_bound = B`, pass an encryption
of `x/B`. The result type and required relinearization-key chain are available
through `CKKSFunctionalBootstrapEvalLUTResult` and
`CKKSFunctionalBootstrapEvalLUTRelinKey`; generate the latter with
`CKKSFunctionalBootstrapEvalLUTKeyGen`.

`CKKSDenseFunctionalBootstrap` and
`CKKSDenseFunctionalBootstrapWithKeyProvider` wire this EvalLUT into the full
Algorithm 1 path:

```text
ModRaise -> CtS -> split -> EvalLUT(real, imag) -> StC
```

The LUT builder for this path is
`CKKSBuildDenseFunctionalBootstrapChebyshevPolynomial<Schedule>`. It accounts
for the C2S input normalization and the existing StC output scaling. The timed
overloads fill `CKKSDenseBootstrapTimings` in the same way as ordinary dense
bootstrapping.

The paper's FHE-friendly path is also implemented. It evaluates a complex
Chebyshev approximation of `exp(2*pi*i*x/2^r)`, applies `r` double-angle
squarings, evaluates Equation (4) as a power polynomial in the encrypted
complex exponential, and extracts the real part with conjugation. The complete
dense API is:

```cpp
using Schedule = TFHEpp::lvl6CKKSDenseFunctionalBootstrapP8Schedule;
auto lut = TFHEpp::CKKSBuildFunctionalBootstrapLUT(values); // 2 <= p <= 8
auto exponential =
    TFHEpp::CKKSBuildFunctionalBootstrapComplexExponentialPolynomial(
        Schedule::exponential_degree, Schedule::functional_double_angle,
        Schedule::functional_input_bound);

using FunctionalKey =
    TFHEpp::CKKSDenseFHEFriendlyFunctionalBootstrapHybridGiantSeededKey<
        Schedule>;
FunctionalKey key;
TFHEpp::CKKSDenseFHEFriendlyFunctionalBootstrapHybridGiantSeededKeyGen<
    Schedule>(key, sk);

TFHEpp::CKKSDenseBootstrapTimings timings;
TFHEpp::CKKSDenseFHEFriendlyFunctionalBootstrapConsumeTimed<Schedule>(
    output, input, lut, exponential, key, timings);
```

The hybrid giant-step/seeded key is the practical in-memory lvl6 variant. The
`ConsumeTimed` call releases each phase's evaluation keys as soon as they are
no longer needed, so that key object is single-use. Use the non-consuming
`CKKSDenseFHEFriendlyFunctionalBootstrap` entry point when the key must remain
reusable and enough memory is available.

The lvl6 p<=8 schedule keeps the already checked 896-bit largest modulus. Its
periodic LUT permits reducing the message ratio from `2^6` to `2`, which leaves
enough budget for a 51-bit EvalLUT scale in TFHEpp's torus/FFT backend. It uses
exponential degree 58, three double-angle squarings, 34-bit plaintext
coefficients, and a degree-7 LUT power polynomial, then restores the 52-bit
scale before StC. Its level budget is:

```text
component q=764
  scale adjustment: 1 bit
  exponential:      6*51 + 34 = 340 bits
  double angle:     3*51      = 153 bits
  LUT powers:       3*51 + 34 = 187 bits
EvalLUT output q=83; two 14-bit StC levels; output q=55
```

The output remains above the 53-bit input level. Since the maximum modulus is
unchanged, the relevant `Parameter-Selection` check remains the existing
`n=32768`, `q=2^896`, sparse-H16 estimate. The functional evaluation key is
different from the ordinary bootstrap key and must be generated or loaded
separately.

### Runtime check

The functional-bootstrap test has an opt-in matched-throughput benchmark:

```sh
./build/test/ckks/ckks_functional_bootstrapping --benchmark
```

The full lvl6 benchmark optionally caches its portable evaluation key:

```sh
./build/test/ckks/ckks_functional_bootstrapping \
  --benchmark-lvl6 /path/to/lvl6-functional-key.bin
```

If the file does not exist, the benchmark generates and saves it. Later runs
load the same key before invoking the single-use consuming bootstrap.

On the current test host (release build), it measured:

| Path | Parameters | Packed values | Total | Amortized |
|------|------------|--------------:|------:|----------:|
| Direct degree-63 functional bootstrap | toy, not security-sized | 16 | 56.89 ms | 3.56 ms/value |
| FHE-friendly complex/double-angle bootstrap | toy, not security-sized | 16 | 22.03 ms | 1.38 ms/value |
| TFHE gate bootstrap | TFHEpp default | 16 | 694.88 ms | 43.43 ms/value |

The real lvl6 run used an Intel i9-9900X (10 cores/20 threads), 62.48 GiB RAM,
and 8 GiB swap. Its arbitrary p=8 LUT passed over all 32768 coefficients:

| Stage | Runtime |
|-------|--------:|
| C2S | 1114.590 s |
| Component split | 25.362 s |
| Real EvalLUT | 1171.920 s |
| Imaginary EvalLUT | 1174.400 s |
| StC | 365.519 s |
| **Functional bootstrap total** | **3851.790 s** |

The resulting amortized time is **117.547 ms/value**, with maximum error
`0.0030527`. The matched-host TFHEpp default gate bootstrap measured
`44.0333 ms/gate`, so the current lvl6 implementation has `0.3746x` TFHE's
amortized throughput: it is **2.67x slower**, not faster. The recorded ordinary
lvl6 CKKS bootstrap takes about `100.55 ms/value`, making this functional path
about 16.9% slower than ordinary bootstrapping as well.

The portable seeded evaluation-key cache is 57,597,521,697 bytes (53.64 GiB).
Generation plus atomic serialization took 47.69 minutes. The consuming run
peaked at 61,899,136 KiB RSS (59.03 GiB), then fell to roughly 11--13 GiB after
C2S keys were released. Thus the 62.48 GiB host needs no additional RAM for
this path; a reusable non-consuming key would require more memory.

On the toy encrypted input, the paper's FHE-friendly EvalLUT is 2.58x faster
than the direct polynomial path. The toy result is useful as a regression but
does not predict the full-size comparison; the measured lvl6 result above is the
relevant throughput result.

## Bounded Sparse Bootstrap Key

The tuned lvl6 bootstrap is validated for a bounded sparse secret key. For the
current tuned schedule:

- `evalmod_k = 18`
- `modraise_mask_bound = 0`
- `bounded_sparse_secret_key_weight = 16`

The bounded sparse key keeps the modulus-raising error inside the EvalMod
interval assumed by the schedule. A dense bootstrap secret can exceed that bound
and invalidate the practical error analysis.

The current lvl6 noise is tuned for the 896-bit bootstrap level. With
`n = 32768`, `q = 2^896`, and `α = 2^-872`, the top-level integer noise
standard deviation is `σ = 2^24`. Lower active levels use the `3.2` floor. This
matches the small-integer Gaussian noise shape used by RNS CKKS libraries and
avoids the large absolute EvalMod noise from the old relative-noise setting.

This schedule is close to Lattigo's local dense `N15QP880H16384H32`
bootstrapping reference in ring size and modulus size: both use `LogN = 15`,
total bootstrap modulus around 880 bits, and dense ternary main secrets.
TFHEpp's current tuned schedule uses inverse-correction bounded-cos EvalMod
(`degree=63`, `double_angle=2`, `inv_degree=7`). OpenFHE's local CKKS
bootstrapping benchmark uses larger `ringDim = 2^16`/`2^17` sets with 50 to 59
bit RNS primes, so the TFHEpp retuned lvl6 path is in the smaller full-size
bootstrapping neighborhood rather than the earlier high-noise 1108-bit
experiment.

The tuned schedule spends 780 bits across C2S, component split, EvalMod, and
STC, leaving `output_log_q = 116` and a post-bootstrap product level of 64 bits.
That is 6 bits above the 58-bit input level. The lvl6 parameter also defines a
`3.2` minimum CKKS integer noise standard deviation so direct low-active-level
encryption does not round the default noise to zero. Reproduce the security
check from the workspace root with:

```sh
cd Parameter-Selection/python
sage -python estimates/CKKS_lvl6.py
```

Validation commands record the requested sparse weight and CKKS noise parameter
in the key directory and reject a mismatch. This prevents accidentally reusing a
bootstrap key generated for a different sparse-weight or security/noise
assumption.

## Key Directories

Full-size CKKS bootstrap keys are too large for convenient in-memory workflows.
Use filesystem key directories for practical runs.

Recommended practical keygen:

```bash
DIR=/tmp/tfhepp_ckks_lvl6_tuned_seeded_streamed
./build-ckks/test/ckks/ckks_bootstrap_validation \
  --lvl6-tuned-seeded-hybrid-streamed-keygen "$DIR"
```

Recommended DD EvalMod keygen:

```bash
DIR=/tmp/tfhepp_ckks_lvl6_tuned_seeded_streamed_dd_evalmod
./build-ckks/test/ckks/ckks_bootstrap_validation \
  --lvl6-tuned-seeded-hybrid-streamed-dd-evalmod-keygen "$DIR"
```

The generated directory contains a manifest and the key material needed by the
streamed filesystem key provider:

- linear transform plan
- rotation-key usage metadata
- EvalMod polynomial metadata
- streamed seeded coefficient-to-slot keys
- packed conjugation key
- seeded EvalMod relinearization keys
- streamed seeded slot-to-coefficient keys

The DD EvalMod directory uses the same linear-transform file layout, but its
EvalMod relinearization files contain compact seeded DD relin keys and its
manifest format is intentionally different. Standard streamed providers reject
DD EvalMod directories, and DD EvalMod providers reject standard streamed
directories.

The library checks manifests before loading or appending to a directory. If a
schedule changes, regenerate the directory.

Relevant public helpers:

```cpp
using Schedule = TFHEpp::lvl6CKKSDenseBootstrapTunedSchedule;
using P = typename Schedule::Param;

TFHEpp::CKKSDenseBootstrapKeyDirectoryOptions options;
options.overwrite_existing = true;

TFHEpp::CKKSDenseBootstrapSeededHybridGiantStreamedKeyGenToDirectory<
    Schedule>(dir, secret_key, {P::α, 0}, options);

TFHEpp::CKKSDenseBootstrapSeededHybridGiantStreamedDDEvalModKeyGenToDirectory<
    Schedule>(dd_dir, secret_key, {P::α, 0}, options);

bool complete =
    TFHEpp::CKKSDenseBootstrapSeededHybridGiantStreamedKeyDirectoryComplete<
        Schedule>(dir);

bool dd_complete =
    TFHEpp::
        CKKSDenseBootstrapSeededHybridGiantStreamedDDEvalModKeyDirectoryComplete<
            Schedule>(dd_dir);

bool manifest_ok =
    TFHEpp::
        CKKSDenseBootstrapSeededHybridGiantStreamedKeyDirectoryManifestMatches<
            Schedule>(dir);

bool dd_manifest_ok =
    TFHEpp::
        CKKSDenseBootstrapSeededHybridGiantStreamedDDEvalModKeyDirectoryManifestMatches<
            Schedule>(dd_dir);
```

## Product Evaluation Keys

The practical multiplication/bootstrap path needs an external product evaluation
key directory. For the seeded streamed path, generate it after the bootstrap key
directory:

```bash
DIR=/tmp/tfhepp_ckks_lvl6_tuned_seeded_streamed
./build-ckks/test/ckks/ckks_bootstrap_validation \
  --lvl6-tuned-seeded-hybrid-streamed-evalkeygen "$DIR"
```

The seeded product eval-key directory contains:

- `seeded_encapsulation_key.bin`
- `seeded_product_relin_key.bin`
- `seeded_post_bootstrap_product_relin_key.bin`
- an eval-key manifest

For the full-DD path, generate a seeded-DD product eval-key directory instead:

```bash
DIR=/tmp/tfhepp_ckks_lvl6_tuned_seeded_streamed_full_dd
./build-ckks/test/ckks/ckks_bootstrap_validation \
  --lvl6-tuned-seeded-hybrid-streamed-full-dd-evalkeygen "$DIR"
```

That directory uses the same seeded encapsulation key file, but stores
double-decomposition product relin keys and a different manifest format:

- `seeded_encapsulation_key.bin`
- `seeded_dd_product_relin_key.bin`
- `seeded_dd_post_bootstrap_product_relin_key.bin`
- an eval-key manifest with the seeded-DD product format

Public helpers:

```cpp
using Schedule = TFHEpp::lvl6CKKSDenseBootstrapTunedSchedule;
using P = typename Schedule::Param;

const bool include_post_bootstrap_product = true;

TFHEpp::CKKSDenseBootstrapSeededProductEvalKeyGenToDirectory<Schedule>(
    eval_key_dir, external_key, bootstrap_key, include_post_bootstrap_product,
    {P::α, 0}, options);

TFHEpp::CKKSDenseBootstrapProductWithSeededHybridGiantStreamedFilesystemSeededEvalKeyDirectoryTimed<
    Schedule>(out, lhs, rhs, key_dir, eval_key_dir, timings);

TFHEpp::CKKSDenseBootstrapPostBootstrapProductWithSeededHybridGiantStreamedFilesystemSeededEvalKeyDirectoryTimed<
    Schedule>(out, lhs, rhs, key_dir, eval_key_dir, timings);

TFHEpp::CKKSDenseBootstrapSeededDDProductEvalKeyGenToDirectory<Schedule>(
    dd_eval_key_dir, external_key, bootstrap_key,
    include_post_bootstrap_product, {P::α, 0}, options);

TFHEpp::CKKSDenseBootstrapProductWithSeededHybridGiantStreamedDDEvalModFilesystemSeededDDEvalKeyDirectoryTimed<
    Schedule>(out, lhs, rhs, dd_key_dir, dd_eval_key_dir, timings);

TFHEpp::CKKSDenseBootstrapPostBootstrapProductWithSeededHybridGiantStreamedDDEvalModFilesystemSeededDDEvalKeyDirectoryTimed<
    Schedule>(out, lhs, rhs, dd_key_dir, dd_eval_key_dir, timings);
```

Use the `PostBootstrapProduct...` function when both inputs are already at the
post-bootstrap product level and the operation must use the
post-bootstrap-product relinearization key.

## DD Relinearization

CKKS multiplication supports both the original full-`Bbar` relinearization key
and DD relinearization keys:

- `CKKSDDRelinKey<P, LogQ>`
- `CKKSSeededDDRelinKey<P, LogQ>`
- `CKKSDDRelinKeyChain<P, StartLogQ, LogDelta, Depth>`
- `CKKSEvalModBoundedCosDDRelinKeys<...>`
- `CKKSDenseBootstrapDDRelinKey<Schedule>`

The generic `CKKSRelinearization` overload dispatches through the key type, so a
caller can use standard relin or DD relin without changing the multiply result
type. The dense in-memory bootstrap provider also accepts
`CKKSDenseBootstrapDDRelinKey<Schedule>`, which makes the EvalMod polynomial,
double-angle, and inverse-correction stages use DD relin.

For filesystem bootstrapping, use the compact seeded DD EvalMod path:

```cpp
using Schedule = TFHEpp::lvl6CKKSDenseBootstrapTunedSchedule;
using P = typename Schedule::Param;

TFHEpp::CKKSDenseBootstrapSeededHybridGiantStreamedDDEvalModKeyGenToDirectory<
    Schedule>(dir, bootstrap_key, {P::α, 0}, options);

TFHEpp::CKKSDenseBootstrapWithSeededHybridGiantStreamedDDEvalModFilesystemKey<
    Schedule>(out, in, dir);
```

`CKKSSeededDDRelinKey` stores one seeded primary encryption per `Bg` digit. At
evaluation time, the provider deterministically expands that primary row,
decomposes it by `Bbar`, and applies the active `Bg` digit. This is compact
seeded storage, not runtime key generation. It exposes DD tuning while avoiding
the full pre-decomposed DD EvalMod key, which is tens of GB at lvl6.

For multiplication followed by bootstrapping, the full-DD filesystem path pairs
the DD EvalMod bootstrap key directory with a seeded-DD product eval-key
directory. Standard seeded product eval-key manifests and seeded-DD product
eval-key manifests are intentionally incompatible so the wrong relin format is
rejected before evaluation.

## Validation Commands

Build the validation binaries:

```bash
cmake --build build-ckks --target ckks_bootstrap_workflow -j2
cmake --build build-ckks --target ckks_bootstrap_validation -j2
```

Run quick tests:

```bash
./build-ckks/test/ckks/ckks_bootstrap_workflow
./build-ckks/test/ckks/ckks_bootstrap_validation --toy-inverse
```

Check the current tuned lvl6 readiness report:

```bash
./build-ckks/test/ckks/ckks_bootstrap_validation --lvl6-tuned-readiness
```

Run the focused in-memory EvalMod diagnostic before regenerating a large key
directory:

```bash
./build-ckks/test/ckks/ckks_bootstrap_validation \
  --lvl6-tuned-inmemory-debug-evalmod
```

The zero-noise variant isolates arithmetic/approximation errors from
key-switching noise:

```bash
./build-ckks/test/ckks/ckks_bootstrap_validation \
  --lvl6-tuned-inmemory-debug-evalmod-zero-noise
```

The readiness report prints:

- schedule levels and EvalMod parameters
- output margin and product slack
- bounded sparse weight
- estimated key/eval-key artifact sizes
- recommended practical commands

Run the target path as separate explicit phases:

```bash
DIR=/tmp/tfhepp_ckks_lvl6_tuned_seeded_streamed

./build-ckks/test/ckks/ckks_bootstrap_validation \
  --lvl6-tuned-seeded-hybrid-streamed-keygen "$DIR"

./build-ckks/test/ckks/ckks_bootstrap_validation \
  --lvl6-tuned-seeded-hybrid-streamed-evalkeygen "$DIR"

./build-ckks/test/ckks/ckks_bootstrap_validation \
  --lvl6-tuned-seeded-hybrid-streamed-run-chained-product-encap "$DIR"
```

Run the DD EvalMod target path as separate explicit phases:

```bash
DIR=/tmp/tfhepp_ckks_lvl6_tuned_seeded_streamed_dd_evalmod

./build-ckks/test/ckks/ckks_bootstrap_validation \
  --lvl6-tuned-seeded-hybrid-streamed-dd-evalmod-keygen "$DIR"

./build-ckks/test/ckks/ckks_bootstrap_validation \
  --lvl6-tuned-seeded-hybrid-streamed-dd-evalmod-evalkeygen "$DIR"

./build-ckks/test/ckks/ckks_bootstrap_validation \
  --lvl6-tuned-seeded-hybrid-streamed-dd-evalmod-run-chained-product-encap "$DIR"
```

Or run the same DD EvalMod target path as one direct command:

```bash
DIR=/tmp/tfhepp_ckks_lvl6_tuned_seeded_streamed_dd_evalmod
./build-ckks/test/ckks/ckks_bootstrap_validation \
  --lvl6-tuned-seeded-hybrid-streamed-dd-evalmod-all "$DIR"
```

Run the full-DD target path as separate explicit phases:

```bash
DIR=/tmp/tfhepp_ckks_lvl6_tuned_seeded_streamed_full_dd

./build-ckks/test/ckks/ckks_bootstrap_validation \
  --lvl6-tuned-seeded-hybrid-streamed-full-dd-keygen "$DIR"

./build-ckks/test/ckks/ckks_bootstrap_validation \
  --lvl6-tuned-seeded-hybrid-streamed-full-dd-evalkeygen "$DIR"

./build-ckks/test/ckks/ckks_bootstrap_validation \
  --lvl6-tuned-seeded-hybrid-streamed-full-dd-run-chained-product-encap "$DIR"
```

Run only the first full-DD product bootstrap when tuning parameters. This uses
the same key directory and product eval-key directory, but avoids the second
post-bootstrap product:

```bash
DIR=/tmp/tfhepp_ckks_lvl6_tuned_seeded_streamed_full_dd
./build-ckks/test/ckks/ckks_bootstrap_validation \
  --lvl6-tuned-seeded-hybrid-streamed-full-dd-run-product-encap "$DIR"
```

Or run the same full-DD target path as one direct command:

```bash
DIR=/tmp/tfhepp_ckks_lvl6_tuned_seeded_streamed_full_dd
./build-ckks/test/ckks/ckks_bootstrap_validation \
  --lvl6-tuned-seeded-hybrid-streamed-full-dd-all "$DIR"
```

Or run the same target path as one direct command:

```bash
DIR=/tmp/tfhepp_ckks_lvl6_tuned_seeded_streamed
./build-ckks/test/ckks/ckks_bootstrap_validation \
  --lvl6-tuned-seeded-hybrid-streamed-all "$DIR"
```

The seeded-streamed path does not require a polling loop. It generates
the complete key directory, generates the product eval-key directory, then runs
the chained product bootstrap.

## Current Full-Size Status

The 896-bit tuned path replaces the earlier 1108-bit experiment. Existing
1108-bit key directories are intentionally incompatible and must be regenerated.
Use `--lvl6-tuned-readiness` for current artifact estimates; the exact key size
depends on the schedule, sparse weight, and seeded/non-seeded layout.

The current readiness report is green for static level budgeting, plaintext
EvalMod approximation, disk advisory checks, and compact seeded DD EvalMod
artifact accounting. Small encrypted regression tests pass for both in-memory DD
EvalMod relin, the compact seeded DD EvalMod filesystem provider, and a tiny
full-DD filesystem chained product bootstrap using seeded-DD product relin
files. The full-size full-DD chained product bootstrap also passes locally with
prebuilt seeded key directories.

The current `Parameter-Selection/python/estimates/CKKS_lvl6.py` estimate for
this exact `n=32768`, `q=2^896`, `α=2^-872` setting passes the 128-bit target:

```text
dense-ternary weakest = 129.183 bits
sparse-H16 weakest = 128.254 bits
```

The full-DD path should therefore be treated as the current practical
full-size, 128-bit-secure dense-bootstrapping baseline, provided bootstrap and
product evaluation keys are regenerated after any parameter change. Local
estimator sweeps still show that `n=65536` gives a larger security margin, but
the active `n=32768` lvl6 path now has both estimator and full-size correctness
coverage.

The measured lvl6 full-DD candidate is:

- `log_delta = 52`, `ratio = 6`, `boot_log_q = 896`
- `c2s_plain_log_delta = 44`, `split_plain_log_delta = 44`
- `slot_to_coeff_plain_log_delta = 14`
- `evalmod_log_scale = 24`, `degree = 63`, `double_angle = 2`
- `inv_degree = 7`, `dd_relin_bgbit = 4`, `dd_relin_bbarbit = 16`

With a prebuilt full-DD seeded directory, the full-size chained product
bootstrap passes the current tolerance:

```text
first_max_error = 0.0952908
chained_max_error = 0.0773982
first_bootstrap_ms = 3.29487e+06
chained_bootstrap_ms = 3.36438e+06
wall_time = 2:16:26
max_rss = 20652888 KiB
```

That run used a regenerated α=2^-872 seeded full-DD key directory of about
35 GB of key material on the local test host. It validates both the first
product-bootstrap step and the chained post-bootstrap product step. The full-DD
chained command is intentionally not part of CI because the key directory is
much larger than the repository artifact limit.

Earlier tuning points are useful negative controls. A `c2s=40`/`split=40`
schedule failed the same streamed full-DD product path with roughly `0.33`
maximum error. An older `log_delta=50`/`c2s=48` schedule produced about `0.064`
single-product error locally, but had less current output-level margin and was
not retained as the active schedule.

The old 1108-bit schedule generated full seeded-streamed keys locally, but the
end-to-end chained product validation failed in encrypted EvalMod because
`α = 2^-850` produced about `2^102` integer noise at the active EvalMod level.
That run should not be used as a correctness reference for the retuned schedule.
Exact timings depend heavily on CPU, OpenMP settings, filesystem, and build
configuration.

## Practical Notes

- Keep CKKS bootstrap and product eval keys outside the repository; `/tmp` or a
  dedicated artifact volume is preferable.
- Regenerate keys when changing schedules, sparse weight, level/scale choices,
  or seeded/non-seeded layout.
- Do not mix bootstrap key directories and product eval-key directories from
  different schedules; manifests are intended to catch this.
- The tuned lvl6 path is the recommended full-size tuning target. Other schedule
  variants are useful for diagnostics, tuning, or comparing key-size/runtime
  tradeoffs.
- The implementation is still experimental. Treat validation commands as part of
  the workflow when changing CKKS code.
