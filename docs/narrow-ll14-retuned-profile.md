# Retuned CPU LL14 profile for the 128-bit classical target

The retuned N=16384 profile passes the tested 128-bit classical
Parameter-Selection screens under BDGL16 and MATZOV. It also passes the original
accuracy gates for real GeLU, Tanh and Sigmoid and their following multiplication,
confirmed by independent integer phase decoding. CUDA is excluded.

| Exposed dense family | log2(q) | Continuous sigma | Integer RMS | BDGL16 bits | MATZOV bits |
|---|---:|---:|---:|---:|---:|
| Direct selectors | 440 | 32768 | 32767.6011 | 131.62 | 129.30 |
| Native automorphism keys | 435 | 4096 | 4095.6011 | 132.37 | 130.08 |
| DD master | 433 | 4096 | 4095.6011 | 133.07 | 130.65 |

The native automorphism primary base is now 32 bits (14 rows). Smaller digits
absorb the increased noise; retaining the historical 44-bit digits failed the
accuracy gate. The DD master retains primary/secondary bases 20/20; evaluation
can publicly regroup primary digits to 40 bits, as in the qualified command
below. Selector scale 57, internal/output scale 44, the original approximation
degrees 12/40/12, and all accuracy gates are preserved.

For the unchanged q56 sparse switch and its public one-signed-position-per-block
secret, the uniform-sparse proxy gives 143.89 bits under BDGL16 and 143.81 under
MATZOV. The separate structured rotation checks use the actual block geometry.
Posterior joint-decoding grids over weights 9–30 reach 137.05 bits. Broader exact
rotation-count grids over multiple guessing dimensions and layouts reach
134.28 bits with rotation-labelled dictionaries, or 129.09 under an optimistic
dictionary-union size assumption. The latter is a sensitivity calculation and
does not demonstrate a feasible attack. The grids are bounded searches.

These are attack-cost estimates under the pinned classical models, scalar GSA
and dictionary/decoder assumptions. The sampler truncates a continuous normal
toward zero; estimates use variance-matched discrete Gaussian proxies. The
paper's structured-RLWE and circular/KDM assumptions remain. This qualification
is not a cryptographic proof or a global lower bound on attack costs. The JSON
field `security_validated: false` retains that distinction; the parameter-screen
and numerical qualification fields record their passing results.

Historical low-noise evaluation-key archives remain historical fixtures. New
archive tags and native noise/base metadata prevent the retuned runner from
accepting them or falling back to historical packed keys. Changing parameters
does not repair weak public material already disclosed for the same secret.
The diagnostic tests reuse retained owner secrets; deployment needs a fresh
secret if weak evaluation-key material for that secret was exposed.

## Reuse saved setup

From the TFHEpp repository root:

```sh
cmake -S . -B build -DENABLE_TEST=ON
cmake --build build --target ckks_narrow_ll14_retuned_setup ckks_narrow_ll14_retuned_native_setup ckks_narrow_ll14_retuned_fbs ckks_narrow_ll14_retuned_archives -j2
build/test/ckks/ckks_narrow_ll14_retuned_archives
```

The setup/replay tools below use the retained LL14 diagnostic checkpoint and
its matching saved switch key and encrypted input. The checkpoint contains
owner secrets and is not a public evaluation-key artifact. Generate replacement
evaluation keys once in a new directory, then keep using those saved files:

```sh
ckks_checkpoint=/path/to/retained/checkpoint-directory
ckks_switch56=/path/to/matching/s2c-switch56.bin
ckks_input=/path/to/matching/candidate-input.bin
ckks_keys=/path/to/new-retuned-key-directory
ckks_output=/path/to/new-replay-directory
mkdir "$ckks_keys" "$ckks_output"
cp "$ckks_input" "$ckks_output/candidate-input.bin"
OMP_NUM_THREADS=8 build/test/ckks/ckks_narrow_ll14_retuned_setup "$ckks_checkpoint" "$ckks_keys/selectors" full
OMP_NUM_THREADS=8 build/test/ckks/ckks_narrow_ll14_retuned_native_setup "$ckks_checkpoint" "$ckks_keys/native"
```

Setup generated no private keys or input ciphertexts in qualification. The
selector/DD setup took 286.96s and the native setup 3.28s on the qualification
machine. A second native setup loaded all eight saved sets and generated zero
new rows. Full selector setup refuses an existing directory. Evaluation never
generates keys or inputs.

The qualified retained-cache evaluation command is:

```sh
OMP_NUM_THREADS=8 build/test/ckks/ckks_narrow_ll14_retuned_fbs \
  "$ckks_checkpoint" "$ckks_switch56" all "$ckks_output" \
  "$ckks_keys/selectors" "$ckks_keys/selectors/dd433-primary20.bin" unused \
  --cache-gib 45 --s2c-cache-gib 4 --direct-rns 1 --s2c-rns 1 \
  --tensor-rns 1 --auto-rns 1 --signed-basis 1 --dd-wide-rns 1 \
  --matrix-tiled 1 --parallel-conjugation 1 --parallel-tail-alignment 1 \
  --parallel-root-columns 1 --native-auto-directory "$ckks_keys/native" \
  --batch-product-trees 1 --dd-primary-bits 40 --dd-level-chunks 1 \
  --matrix-retained-rows 1 --matrix-joint 1 --s2c-baby-workers 8 \
  --parallel-table-ntt 1 --s2c-group-workers 8 --basis-accumulate 1 \
  --basis-group-tensors 1 --native-auto-batch 1 --dd-windowed 1
```

Observed hot evaluations used eight pinned P cores of an i9-12900K:

| Function | Maximum error | Following-square error | Evaluation |
|---|---:|---:|---:|
| GeLU | 4.74e-8 | 3.21e-7 | 5.98s |
| Tanh | 3.61e-8 | 6.41e-8 | 10.96s |
| Sigmoid | 2.07e-8 | 3.74e-8 | 5.83s |

All 8192 slots pass the original 2^-23 function gate and amplitude-aware
post-multiplication gate. These observations exclude saved-key/cache setup;
background estimator jobs overlapped some measurements. Hardware and parameters
differ from the paper. The earlier historical tenfold-latency qualification
remains a measurement of its older low-noise profile.

## Evidence

[The qualification record](narrow-ll14-retuned-profile-validation.json) binds
source revisions, actual saved-key hashes, input/output hashes, commands,
attack estimates, public archive rejection checks, and independent noise and
phase checks. Sampled key residual RMS values were 32930.41 (selectors), 4087.20
(native) and 4110.49 (DD), each checked over 16384 coefficients. This verifies
noise scale on the saved samples, not an exact PMF or security proof.

The full Parameter-Selection assessment is
`python/estimates/narrow_ckks_ll14_retuned_dense_security.json` in that repository.
The qualification record pins its commit and SHA-256. The lattice estimator is
pinned to `3e48ef421ec256afddb3e7d2249a77eab6e9ba12`.
Its [attack documentation](https://lattice-estimator.readthedocs.io/en/latest/algorithms/lwe-primal.html)
describes the underlying primal attack models.
