# CKKS

TFHEpp includes experimental CKKS support in `include/ckks/ckks.hpp`. The current
target path is a leveled CKKS multiplication followed by dense CKKS
bootstrapping. The current storage-practical full-size target is seeded hybrid
giant streamed bootstrapping with compact seeded double-decomposition (DD)
EvalMod relinearization.

The most practical full-size configuration is the tuned lvl6 seeded hybrid
giant streamed bootstrap path:

- schedule: `TFHEpp::lvl6CKKSDenseBootstrapTunedSchedule`
- ring degree: `n = 32768`
- torus storage: `896` bits
- boot modulus: `boot_log_q = 888`
- input level: `input_log_q = 48`
- output level: `output_log_q = 108`
- scale bits: `log_delta = 40`
- lvl6 CKKS noise: `α = 2^-886`
- low-level CKKS noise floor: `σ >= 3.2`
- bounded bootstrap secret weight: `16`

CKKS APIs and schedules are compile-time types. Changing a schedule, level,
scale, sparse key weight, or key layout requires regenerating the relevant
bootstrap and product evaluation keys.

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

## Bounded Sparse Bootstrap Key

The tuned lvl6 bootstrap is validated for a bounded sparse secret key. For the
current tuned schedule:

- `evalmod_k = 18`
- `modraise_mask_bound = 0`
- `bounded_sparse_secret_key_weight = 16`

The bounded sparse key keeps the modulus-raising error inside the EvalMod
interval assumed by the schedule. A dense bootstrap secret can exceed that bound
and invalidate the practical error analysis.

The current lvl6 noise is tuned for the 888-bit bootstrap level. With
`n = 32768`, `q = 2^888`, and `α = 2^-886`, the top-level integer noise
standard deviation is `σ = 2^2`. Lower active levels use the `3.2` floor. This
matches the small-integer Gaussian noise shape used by RNS CKKS libraries and
avoids the large absolute EvalMod noise from the old relative-noise setting.

This schedule is close to Lattigo's local dense `N15QP880H16384H32`
bootstrapping reference in ring size and modulus size: both use `LogN = 15`,
total bootstrap modulus around 880 bits, and dense ternary main secrets.
TFHEpp's current tuned schedule uses inverse-correction bounded-cos EvalMod
(`degree=34`, `double_angle=4`, `inv_degree=5`). OpenFHE's local CKKS
bootstrapping benchmark uses larger `ringDim = 2^16`/`2^17` sets with 50 to 59
bit RNS primes, so the TFHEpp retuned lvl6 path is in the smaller full-size
bootstrapping neighborhood rather than the earlier high-noise 1108-bit
experiment.

The tuned schedule spends 780 bits across C2S, component split, EvalMod, and
STC, leaving `output_log_q = 108` and a post-bootstrap product level of 68 bits.
That is 20 bits above the 48-bit input level. The lvl6 parameter also defines a
`3.2` minimum CKKS integer noise standard deviation so direct low-active-level
encryption does not round the default noise to zero. Reproduce the security
check from the workspace root with:

```sh
cd Parameter-Selection/python
sage -python estimates/CKKS_lvl6.py
```

Validation commands record the requested sparse weight in the key directory and
reject a mismatch. This prevents accidentally reusing a bootstrap key generated
for a different sparse-weight assumption.

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

The DD EvalMod filesystem path changes only EvalMod relinearization. Product
evaluation keys still use the existing seeded relin format unless a future
product-DD path is added.

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

The 888-bit tuned path replaces the earlier 1108-bit experiment. Existing
1108-bit key directories are intentionally incompatible and must be regenerated.
Use `--lvl6-tuned-readiness` for current artifact estimates; the exact key size
depends on the schedule, sparse weight, and seeded/non-seeded layout.

The current readiness report is green for static level budgeting, plaintext
EvalMod approximation, disk advisory checks, and compact seeded DD EvalMod
artifact accounting. Small encrypted regression tests pass for both in-memory DD
EvalMod relin and the compact seeded DD EvalMod filesystem provider.

The remaining full-size validation step is generating the lvl6 tuned DD EvalMod
directory, generating the seeded product eval-key directory, and running the
chained product bootstrap command. That run is intentionally not part of CI
because the key directory is much larger than the repository artifact limit.

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
