# CKKS

TFHEpp includes experimental CKKS support in `include/ckks/ckks.hpp`. The current
practical path is a leveled CKKS multiplication followed by dense CKKS
bootstrapping, using Double Decomposition (DD) through the underlying TFHEpp
external-product parameter path.

The most practical full-size configuration is the tuned lvl6 seeded hybrid
giant streamed bootstrap path:

- schedule: `TFHEpp::lvl6CKKSDenseBootstrapTunedSchedule`
- ring degree: `n = 32768`
- boot modulus: `boot_log_q = 1152`
- input level: `input_log_q = 60`
- output level: `output_log_q = 156`
- scale bits: `log_delta = 52`
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
| `CKKSDenseBootstrapSchedule<...>` | Compile-time dense bootstrap parameter schedule. |
| `Schedule::InputCiphertext` | Ciphertext type accepted by a dense bootstrap schedule. |
| `Schedule::OutputCiphertext` | Ciphertext type produced by a dense bootstrap schedule. |

Useful aliases are defined near the end of `include/ckks/ckks.hpp`, including:

- `lvl6CKKSDenseBootstrapInput`
- `lvl6CKKSDenseBootstrapOutput`
- `lvl6CKKSDenseBootstrapTunedSeededHybridGiantStreamedFilesystemKeyProvider`

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

Validation commands record the requested sparse weight in the key directory and
reject a mismatch. This prevents accidentally reusing a bootstrap key generated
for a different sparse-weight assumption.

## Key Directories

Full-size CKKS bootstrap keys are too large for convenient in-memory workflows.
Use filesystem key directories for practical runs.

Recommended practical keygen:

```bash
DIR=/tmp/tfhepp_ckks_lvl6_tuned_seeded_streamed
./build-ckks/test/ckks_bootstrap_validation \
  --lvl6-tuned-seeded-hybrid-streamed-keygen "$DIR"
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

bool complete =
    TFHEpp::CKKSDenseBootstrapSeededHybridGiantStreamedKeyDirectoryComplete<
        Schedule>(dir);

bool manifest_ok =
    TFHEpp::
        CKKSDenseBootstrapSeededHybridGiantStreamedKeyDirectoryManifestMatches<
            Schedule>(dir);
```

## Product Evaluation Keys

The practical multiplication/bootstrap path needs an external product evaluation
key directory. For the seeded streamed path, generate it after the bootstrap key
directory:

```bash
DIR=/tmp/tfhepp_ckks_lvl6_tuned_seeded_streamed
./build-ckks/test/ckks_bootstrap_validation \
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

## Validation Commands

Build the validation binaries:

```bash
cmake --build build-ckks --target ckks_bootstrap_workflow -j2
cmake --build build-ckks --target ckks_bootstrap_validation -j2
```

Run quick tests:

```bash
./build-ckks/test/ckks_bootstrap_workflow
./build-ckks/test/ckks_bootstrap_validation --toy-inverse
```

Check the current tuned lvl6 readiness report:

```bash
./build-ckks/test/ckks_bootstrap_validation --lvl6-tuned-readiness
```

The readiness report prints:

- schedule levels and EvalMod parameters
- output margin and product slack
- bounded sparse weight
- estimated key/eval-key artifact sizes
- recommended practical commands

Run the practical path as separate explicit phases:

```bash
DIR=/tmp/tfhepp_ckks_lvl6_tuned_seeded_streamed

./build-ckks/test/ckks_bootstrap_validation \
  --lvl6-tuned-seeded-hybrid-streamed-keygen "$DIR"

./build-ckks/test/ckks_bootstrap_validation \
  --lvl6-tuned-seeded-hybrid-streamed-evalkeygen "$DIR"

./build-ckks/test/ckks_bootstrap_validation \
  --lvl6-tuned-seeded-hybrid-streamed-run-chained-product-encap "$DIR"
```

Or run the same practical path as one direct command:

```bash
DIR=/tmp/tfhepp_ckks_lvl6_tuned_seeded_streamed
./build-ckks/test/ckks_bootstrap_validation \
  --lvl6-tuned-seeded-hybrid-streamed-all "$DIR"
```

The practical seeded-streamed path does not require a polling loop. It generates
the complete key directory, generates the product eval-key directory, then runs
the chained product bootstrap.

## Current Practical Reference

On the tuned lvl6 seeded-streamed path, the readiness report currently estimates
about 33.5 GB of seeded bootstrap/product artifacts with the product eval keys.
The non-seeded hybrid path is roughly twice that size.

A full local validation of
`--lvl6-tuned-seeded-hybrid-streamed-all /tmp/tfhepp_ckks_lvl6_tuned_seeded_streamed`
completed with:

- bootstrap key directory complete: `159/159` files
- total files including product eval keys: `165`
- first product bootstrap max error: about `0.00133`
- chained product bootstrap max error: about `0.00140`
- first bootstrap time: about `1.42e6 ms`
- chained bootstrap time: about `1.42e6 ms`

Exact timings depend heavily on CPU, OpenMP settings, filesystem, and build
configuration.

## Practical Notes

- Keep CKKS bootstrap and product eval keys outside the repository; `/tmp` or a
  dedicated artifact volume is preferable.
- Regenerate keys when changing schedules, sparse weight, level/scale choices,
  or seeded/non-seeded layout.
- Do not mix bootstrap key directories and product eval-key directories from
  different schedules; manifests are intended to catch this.
- The tuned lvl6 path is the recommended practical target. Other schedule
  variants are useful for diagnostics, tuning, or comparing key-size/runtime
  tradeoffs.
- The implementation is still experimental. Treat validation commands as part of
  the workflow when changing CKKS code.
