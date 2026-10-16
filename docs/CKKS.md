# CKKS

## CPU narrow bootstrapping from 2026/2119

The CPU implementation uses Double Decomposition for the product tree and
implements the paper's head/tail selectors, NTT permutations, frequency batching,
packed-key reconstruction, and functional bootstrapping. CUDA is excluded.
The native-modulus port uses different parameters and hardware from the paper;
overall 128-bit security and matching paper latency are not claimed.

The retained N16384 mask42 pipeline independently passes real and componentwise
complex GeLU, Tanh and Sigmoid, including a subsequent multiplication. Both FFT
and the optional combined NTT automorphism/DD backend pass. The full default CPU
regression also passes against the qualified headers now installed. The NTT backend is now
available in production headers; its cache lifecycle preserves the original
serialized key format. Evidence: `narrow-dd-combined-ntt-production-integration-validation.json`,
`narrow-dd-combined-default-regression-validation.json`, and the two
`narrow-dd-ll14-packed-mask42-combined-ntt-*-validation.json` reports.

| Function | Real NTT maximum error | Real evaluation | Complex NTT maximum error |
|---|---:|---:|---:|
| GeLU | 3.48e-8 | 167.7 s | 3.95e-8 |
| Tanh | 3.77e-8 | 406.2 s | 4.76e-8 |
| Sigmoid | 1.71e-8 | 169.8 s | 1.87e-08 |

These are observed timings with eight OpenMP threads on an Intel i9-12900K;
some runs overlap other work. They exclude loading keys and preparing caches.
Paper Table 7 reports 2.71/4.51/2.79 seconds on 48 CPU cores of a dual EPYC 7K62
server. This is numerical qualification, not a matched performance benchmark.
See `narrow-dd-mask42-cpu-platform-observation.json` and
`narrow-dd-mask42-security-scope-review.json` for the platform and security scope.

Build the saved-fixture runners with:

```sh
cmake -S TFHEpp -B TFHEpp/build -DENABLE_TEST=ON
cmake --build TFHEpp/build --target ckks_narrow_ll14_ntt_fbs ckks_narrow_ll14_ntt_complex_fbs -j4
```

Both accept `checkpoint-directory dd-chain-file switch56-file function
output-directory reconstructed-selector-directory saved-dd433 saved-packed435`.
The output directory must contain the retained `candidate-input.bin`.
These are fixed LL14 validation fixtures: selectors379/57 with mask42,
factor324/44, DD324/280/236/192/148, output101/44 and square57/44.
The FFT reference runners `ckks_narrow_ll14_packed_fbs` and
`ckks_narrow_ll14_packed_complex_fbs` accept the first six arguments.
No runner generates keys, encrypts a new input, or reconstructs selectors.
The direct-selector reference is `ckks_narrow_ll14_retained_fbs` with its
separately qualified43-bit chain and direct-selector directory.

Generate setup once using `CKKSNarrowRefreshSetupKeyGen`, then generate/save the
packed public key with `CKKSNarrowPackedKeyGen` and `CKKSNarrowSavePackedKey`.
Retain the reconstructed selector files. For the saved mask42 LL14 archive,
reconstruction uses `ckks_narrow_packed_archive --reconstruct-mask42 ARCHIVE
NEW_SELECTOR_DIRECTORY 14 NEW_SCRATCH_DIRECTORY`. The measured one-time
reconstruction takes95.7min and retains34.95GB of selectors. Subsequent refreshes
reuse them. Earlier mask34 failures remain archived as diagnostic evidence.

For the retained mask42 schedule, reattach its deferred-rescale/public-mask
factor source after loading. Include `ckks/narrow_dd_ntt.hpp` and
`ckks/narrow_auto_ntt.hpp`, then select the cached backend once:

```cpp
CKKSNarrowUseStagedFFT(key.narrow, 40);
CKKSNarrowUseAuxiliaryNTTDD(key.narrow, master);
CKKSNarrowUsePublicNTTAutomorphisms<337>(key.narrow, packed);
```

To evaluate independent DD product-tree pairs concurrently, opt in with
`CKKSNarrowUseAuxiliaryNTTDD<100, 32, true>(key.narrow, master)`.
The real saved-fixture runner accepts `--parallel-pairs` after its eight
arguments. Use `OMP_NUM_THREADS=8 OMP_MAX_ACTIVE_LEVELS=1` for the recorded
comparison: lower tree levels distribute independent multiplications across
the team, while a single remaining pair uses its internal parallelism.
The prepared keys are shared read-only and failures propagate after joining
the workers. The existing sequential pair schedule remains the default.
See `narrow-dd-cpu-gap-optimization-validation.json` for the measured scope.
In a saved-fixture GeLU comparison, the tree stage drops from 33.0 s to 16.8 s
and the full refresh from 134.1 s to 112.0 s. Both full output archives are
byte-identical, including the subsequent multiplication. These are single
observations with host variability, not a matched paper benchmark or a 10x
speedup. Selector matrix work and transforms remain the largest costs.

The real runner also accepts `--small-prime-factors`. This opt-in backend uses
four primes below 2^30, 52-bit selector digits, native AVX2 Montgomery products,
and exact signed CRT. It preserves the ciphertext modulus, scale, saved keys,
selectors and approximation. Derived spectrum payload is unchanged. Four-column
remainders stay vectorized, and raw reads overlap the current chunk's evaluation.
Include `ckks/narrow_factor_ntt32.hpp` and use
`CKKSNarrowDeferredRescaleFactorSource32<TargetSchedule, TableSchedule, SelectorSchedule>`
in place of the corresponding deferred factor source; retain the public-mask
normalization wrapper.

Saved-fixture real GeLU/Tanh/Sigmoid refreshes with parallel DD pairs take
98.0/249.2/99.8 seconds in the recorded runs and preserve all output bytes,
including the subsequent multiplication. GeLU's preceding wide-prime paired
run took 112.0 s. These are unpinned observations; the 10x gap remains open.
The first synchronous integration was slower and remains recorded separately.
See `narrow-dd-small-prime-arithmetic-validation.json`,
`narrow-dd-small-prime-factor-validation.json`, and
`narrow-dd-small-prime-factor-synchronous-observation.json`.

Here `master` is the matching saved public DD433 key and `packed` is the matching
saved packed435 key. Keep `key` alive across evaluations to reuse all caches.
Call `CKKSNarrowPrepareCoefficients` once per input when several functions share
it, then call `CKKSNarrowFunctionalRefreshPrepared` for each function.
These evaluation entry points do not generate keys or reconstruct selectors.
Serialized refresh keys omit external selector/factor sources and derived caches;
reattach them after loading. Explicit FFT selection clears the corresponding NTT
backend, while automorphism selection preserves the independent DD backend.
The generic `CKKSNarrowRefreshKeyGen` instead generates expanded direct selectors.

The current prepared DD relinearization stores inverse-FFT products in native
64-bit rows and sign-extends them during torus recombination. This reduces
product workspace by nine times for the 576-bit profile without changing the
FFT order, rounding, decomposition bases, or mantissa guard. All seven N32768 retained-key DD levels match both frozen FFT and coefficient-
domain switching (1,376,256 coefficients). Signed recombination passes an
independent integer oracle. N512 all-function replay preserves all three output
archives and all 36 retained files. Paired component trials take 25–45% of the
previous switching time; these are component timings, not full-refresh speedups.
See `narrow-dd-native-products-full-levels-validation.json`,
`narrow-dd-native-products-recombine-validation.json`, and
`narrow-dd-native-products-small-refresh-validation.json`. A preceding N32768 replay with these storage changes passes all three output archive
comparisons, independent phase checks and all 37 retained-file hashes
(`narrow-dd-native-products-full-replay-validation.json`).

Prepared automorphism switching now also retains inverse-FFT products in
native64 rows until signed torus recombination. Retained N32768 Q523 prototype
controls preserve all outputs and take about 72% of the previous switching time.
Current-source direct switches and N32768 head/tail factors pass exact archive
comparisons and 720,896 slot checks at both selector widths. N512 all-function
refresh preserves all three output archives and all 36 retained files. Full
checkpoint and selected factor files are unchanged. See
`narrow-dd-auto-native-products-integration-validation.json`. Neither this
component timing nor the DD timing establishes a full-refresh speedup; complete
N32768 retained replay passes for these source changes. Shared coefficient
preparation and 22-root generation take 2012.36 s (33.54 min), versus 2181.40 s
for the preceding run. This observed 7.75% reduction includes host/compiler/runtime
variation. Reconstruction adds 2.54 s for Sigmoid, 2.45 s for GeLU and 8.74 s for
Tanh; these additions are not standalone functional-bootstrapping latencies.
Paper latency and security remain unqualified. See
`narrow-dd-native-products-full-runtime-profile.json`.

The public frequency-basis API accepts `complex_input=false` to evaluate the
real component and return a real-valued result. Current N512 retained-input
checks call this API separately for Sigmoid, GeLU and Tanh, using 4, 4 and 9
positive-frequency columns respectively. Independent integer phase decoding,
post-refresh multiplication and normalization checks pass, with all 36 retained
files unchanged and no new keys. These are small-ring API checks; complete
N32768 real-component API and paper latency remain unqualified. The completed
full diagnostic instead shares 22 columns across all three functions on both
complex components, so its total is not a Table 7 single-function latency.
See `narrow-dd-current-real-basis-public-api-validation.json`.

Narrow S2C now uses an evaluation-owned lazy prepared FFT Galois cache when
its mantissa guard permits. The existing BSGS boundaries release derived keys;
serialized keys and inputs can be reused. A saved N32768 adapter comparison
preserves all dense and sparse-secret coefficients and takes 198.12 s versus
409.12 s for the raw path, including cache construction (one trial per path).
Integrated N512 all-function replay preserves all three archives, passes independent
phase/product checks and leaves all 36 retained files unchanged. The combined integrated
N32768 replay now passes exact output, independent accuracy and all 37 retained-file
checks. Paper latency remains unqualified. See
`narrow-dd-s2c-prepared-integration-validation.json`.

The current S2C/tensor executable also reuses the previously created N512
selector56 disk cache: all three archives match the current raw-path replay,
all 36 checkpoint files and 35 cache files remain unchanged, and selector NTT
time is zero. A full N32768 selector56 cache contains about 399 GB of spectrum
payload and still requires raw-source reads for validation. The cache exceeds
host RAM; its additional I/O can offset saved transforms. A complete full-size
cache has not been created or timed. Alternating current N32768 head30/tail4
40-record trials give warm prepared-load medians of 0.366 s versus raw-load plus
NTT medians of 0.322/0.312 s, including raw/payload verification. Exact spectra
match. These concurrent component observations support retaining the raw path
on this host (`narrow-dd-selector56-full-chunk-cost-validation.json`). See
`narrow-dd-s2c-tensor-existing-cache-replay-validation.json` and
`narrow-dd-selector56-current-cache-reuse-cost-audit.json`.

Tensor-product FFT evaluation now also accumulates native signed 64-bit
inverse-FFT products directly, avoiding expansion into multi-limb torus values
before reading them back as signed integers. This reduces each worker's product
buffer ninefold for the 576-bit profile. FFT batching, mantissa guards, signed
rounding and the coefficient fallback are preserved. Independent integer/NTT
controls and the frozen tensor implementation agree for N512 and N32768,
including all seven tree levels, asymmetric moduli, negative inputs, rounding
boundaries and dense large digits. Alternating N32768 component trials take
91.4–95.6% of the old time. Combined N512 all-function replay passes exact archive,
independent accuracy and all 36 retained-file checks without generating keys.
Combined N32768 replay with this change and the prepared S2C cache passes all
three exact archive comparisons, independent phase/product checks and all 37
retained-file hashes. Keys and input are reused. Shared S2C and 22-root preparation
take 1910.13 s (31.84 min); total evaluator wall time is 32:12.67. The observed
5.08% reduction from the preceding run includes host and concurrent-work variation.
This shared all-function complex workload is not a Table 7 single-function latency;
paper performance and security remain unqualified. See
`narrow-dd-s2c-tensor-full-replay-validation.json` and
`narrow-dd-s2c-tensor-full-runtime-profile.json`.
Evidence: `narrow-dd-tensor-native-products-integration-validation.json`,
`narrow-dd-tensor-native-products-timing-validation.json` and
`narrow-dd-tensor-native-products-small-refresh-validation.json`.

The frozen N16384 packed-key creation now exits 0 and passes independent
full-ring refresh and post-refresh multiplication checks. Reconstruction took
97118.6 s and refresh took 900.43 s; maximum refresh error was 4.6877e-6 against
this fixture's 2^-16 gate. Both creation and replay exited 0, passed independent phase checks and produced
identical refresh/square archives. All 38 retained files are unchanged
(`narrow-dd-full-packed-retained-completion-validation.json`). This is a
standard refresh check, not the high-precision functional result in Table 7.

Including `tfhe++.hpp` exposes the narrow refresh API and packed-key save/load
functions. With BLAKE3 enabled, it also exposes the optional prepared selector-file
cache. Consumer compilation passes with and without BLAKE3
(`narrow-dd-public-header-validation.json`). Packed keys can be saved once and
reloaded for later evaluations.

`CKKSReduceAutoKey<P, TargetQ>` and `CKKSReduceGaloisKey<P, TargetQ>`
(`ckks/key_projection.hpp`, also exposed by `tfhe++.hpp`) reuse retained public
keys at a lower modulus. They select the matching gadget rows and reduce their
coefficients modulo the destination modulus; they do not rescale or generate
keys. The ring, secret and primary radix must match, and the destination modulus
cannot exceed the source modulus. This also applies to the automorphism-key
alias used for secret switching. Changing an evaluation scale still requires
an accuracy check. Saved N512 controls at scale 47 pass independent integer
gadget-phase checks and exact raw/prepared switches for 523→345, 197→101,
197→100 and 197→197. A scale-41 trial with the retained radix-8 key fails the
23-bit conjugation gate; a smaller refresh profile is not thereby qualified.
See `narrow-dd-auto-key-modulus-reuse-validation.json`.

A retained N16384 radix-4 key also passes public reduction from Q345 to Q101
and Q109. At scale 41, direct conjugation has 1.71e-6 error and fails the
23-bit gate; an eight-bit public auxiliary lift using the reduced Q109 key has
9.47e-9 error and passes. Raw/prepared ciphertexts match exactly; independent
integer phases confirm all 8192 slots, and the retained checkpoint is unchanged
(`narrow-dd-auto-key-modulus-full-validation.json`). These are conjugation
controls, not a complete functional-refresh result.

`CKKSNarrowRescaledSelectorSource<Target, Source>` in
`ckks/narrow_selector_projection.hpp` wraps retained selectors and performs public
rounded rescaling for each requested raw chunk. Ring and selector geometry must
match; the source modulus/scale must exceed the destination by the same number
of bits. Original prepared spectra cannot be reused after conversion. The
wrapper preserves raw prefetch support and keeps the source alive. N16384
379/49→371/41 head30/tail4 edge samples pass 131,072 independent signed-rounding
checks and independent integer phase decoding. Their maximum slot perturbation
is 9.60e-9. Both sampled source files and the consumed checkpoint secret prefix
remain unchanged. This qualifies four records, not every selector or a complete
refresh (`narrow-dd-selector-rescale-reuse-full-validation.json`). The smaller
functional profile and its remaining unimplemented/unchecked gates are recorded
in `narrow-dd-ll14-retained-functional-candidate.json`.

`CKKSNarrowRefreshSchedule` accepts an optional final S2C auxiliary-scale argument;
its default remains the narrow schedule's auxiliary scale. Existing N512 default
all-function output archives and all 36 retained files remain identical, with
independent phase checks passing. Byte-for-byte source restoration confirms that
existing evaluator bodies are unchanged; this preserves the preceding full
qualification scope without claiming another full default computation
(`narrow-dd-s2c-independent-scale-default-source-projection.json`).

`CKKSRebindCompactLinearInputScale` reuses encoded compact plaintext plans at a
different input scale, and also at a different modulus when both levels have
at least 64 bits. Plaintext precision stays fixed. A fresh public input over
[-4,4] under retained N16384 private keys rejects S2C scale44/auxiliary0 at
3.11e-7 error. Auxiliary8 passes independent full-vector decoding at 2.63e-9
error, with S2C taking 16.63 s in this component trial. Its Q93/Q56 Galois keys
are reductions of the retained Q345 dense-key Galois set. The Q56 structured
switch key could not be raised from the old Q48 key: exactly one new evaluation
key was generated under retained private keys and saved for reuse. No packed
selectors were reconstructed. See
`narrow-dd-s2c-independent-scale-aux8-full-validation.json`; this component
result does not qualify the full functional profile or paper latency.

`makeCKKSCompactDDRelinKeyAtLevel<TargetQ>` publicly reduces retained compact DD
keys. It selects matching primary gadget rows, reconstructs signed secondary
digits, reduces encrypted coefficients and decomposes them at the destination
level. N512 retained Q244→Q197/Q132/Q101/Q100/Q244 controls pass 46,080 integer
key-coefficient checks, archive roundtrips, exact raw/prepared before-rescale
squares and independent phase decoding. Private and evaluation keys are reused.
Full-ring qualification follows below
(`narrow-dd-dd-key-modulus-reuse-validation.json`).

The proposed N16384 DD chain is now component-qualified. Retained Q296 keys
reduce to Q296/Q255/Q214/Q173/Q101 with 8,585,216 independent integer encrypted-
coefficient comparisons. Exactly one missing Q337 evaluation key was generated
under the retained dense private secret and saved. All six before-rescale
multiplication checks have identical raw/prepared ciphertexts and independently
pass the 23-bit gate. The complete production 32-leaf, five-level tree also
matches exactly and independently has 4.87e-8 error. All 114,688 output phase
coefficients are checked; the original checkpoint remains unchanged. The saved
chain and post-product key are reusable
(`narrow-dd-dd-retained-tree-full-validation.json`).

`ckks_narrow_ll14_retained_fbs` is an experimental real-input replay runner. It
loads those saved evaluation keys, the saved Q56 S2C switch key and the retained
packed-checkpoint Galois keys; it rescales raw selectors in 40-record chunks and
uses 60-bit matrix scratch digits under the existing exact CRT bound. It saves
one public encrypted input for reuse across functions and runs the direct
weighted-frequency series with before-rescale DD. The output is publicly reduced
to Q101 at scale41 and squared using the retained Q101 key, leaving Q60/scale41.
The sequential full functional validation is pending, beginning with Sigmoid;
this runner does not yet establish paper parameter, security or latency matching.
See `narrow-dd-ll14-retained-fbs-validation.json`.

`CKKSNarrowRefreshPrepared` and `CKKSNarrowFunctionalRefreshPrepared` evaluate
an existing `RefreshSchedule::Coefficients` ciphertext produced by
`CKKSNarrowPrepareCoefficients` with the same refresh key. They preserve that
ciphertext, allowing diagnostics and multiple functions to reuse S2C. The
standard input APIs still perform preparation automatically. The diagnostic
benchmark now uses prepared evaluation and reports full refresh time as
`coefficient_preparation_seconds + prepared_evaluation_seconds`; diagnostic
decryption and first-factor checks are excluded. Exact standard/complex/real
controls pass on direct and packed compact configurations. Existing N512
checkpoint replays pass independent checks, preserve all 37 files and match
the earlier executable's refresh and square archives exactly, without new keys
(`narrow-dd-prepared-refresh-validation.json`). The standard prepared-refresh
API also passes N16384 retained replay with
identical earlier archives. Full-size prepared functional refresh remains
unqualified.

Packed selector diagnostics now check candidates in parallel with a private
slot buffer per worker and finite per-candidate maxima. Worker validation
exceptions propagate after the parallel region. All candidate slots are still
checked, and each completed group prints progress; reports include
`selector_validation_seconds` within the enclosing setup time. N512 retained
serial/four-worker controls pass 152320 selector-slot comparisons in 35 groups,
produce identical selector maxima and match the older refresh/square archives.
All 37 files are unchanged. Validation took 0.6695 s serially and 0.2307 s with four
workers under concurrent load. N16384 replay now passes 156008448 comparisons
in 36 groups and the unchanged output/integrity gates. Its selector diagnostics
took 218.07 s. Evidence:
`narrow-dd-parallel-selector-validation.json` and
`narrow-dd-current-qualified-packed-full-replay-validation.json`.

`Parameter-Selection/python/estimates/narrow_ckks_packed_qualified_replay.py`
runs one evaluation against a completed retained-checkpoint baseline. It pins
the small-qualified executable, verifies saved baseline output hashes, checks
all selector progress groups and inclusive refresh timing, independently decodes
refresh and multiplication, and hashes every retained file after evaluation.
An optional `--wait-runner PID:start_ticks` waits for the existing baseline
controller to exit before starting, so large evaluations do not overlap. The
strict N512 workflow and current N16384 retained replay both pass without
generating keys. The N16384 controller/evaluator exited 0, all 38 files are
unchanged and both output archives match the frozen baseline. The observed
evaluator wall time was 434.05 s versus 2099.85 s for the older executable; refresh
took 209.00 s including 46.77 s preparation and 162.23 s prepared evaluation. Earlier
core optimizations, compiler and concurrent I/O conditions prevent attribution
to only the latest benchmark changes or a controlled speedup claim. This
standard refresh does not qualify Table 7 functional latency. See
`narrow-dd-current-qualified-packed-full-replay-validation.json`.

Packed substitution preparation also accepts a budget for FFT spectrum data:
`CKKSNarrowPreparePackedSubstitutionsWithinBudget<Schedule>(packed, bytes)`
selects as many deepest extraction levels as fit. Zero clears the derived
cache. `CKKSNarrowPackedSubstitutionCacheBytesPerDepth` and
`CKKSNarrowPackedSubstitutionCacheDepthForBudget` provide allocation-free size
and depth calculations. The budget covers spectra; raw keys and reconstruction
buffers need their own memory allocation. Existing explicit-depth preparation
and benchmark defaults remain unchanged.

Retained N512 controls for zero, two, four and nine cached levels preserve all
2048 ordered terms and 2097152 ciphertext coefficients, with unchanged packed
key/reference files. Extraction took 7.49 s with two levels and 4.36 s with all
nine; these are component observations. The N16384 profile needs10602676224
spectrum bytes for all14 levels (about9.9GiB), within a16GiB spectrum budget.
Full-size reconstruction and its speedup remain unqualified. Evidence:
`narrow-dd-packed-cache-budget-validation.json`.

The existing FFT precision guard permits a 12-bit DD secondary base at the
N32768 Q244/Q197 levels. Retained public-key controls pass 393216 output
coefficient comparisons against the coefficient-domain reference exactly.
Higher levels retain secondary11. The two lower-level controls do not qualify
a whole-refresh speedup. A separate secondary11/12 control now passes the
same exact reference comparisons and saves 20507232 cache bytes across these
two levels. Three-switch observations were 0.1614/0.1374 s at Q244 and
0.1154/0.1089 s at Q197; preparation timings have different scopes. The
refresh policy remains unchanged. Evidence:
`narrow-dd-secondary12-low-levels-validation.json` and
`narrow-dd-secondary11-v12-low-levels-validation.json`.

Factor evaluation reuses its two immutable public NTT plans across selector
groups and frequency tiles. Current-source retained N512 selector56 head/tail
controls preserve both earlier factor archives exactly and pass 11264 slot
checks; all retained files remain unchanged. This change is qualified at the
component scope; full-size refresh timing remains unmeasured. Evidence:
`narrow-dd-factor-plan-reuse-validation.json`.

The scalar NTT backend now uses AVX2 lazy butterflies when compiled with
AVX2 and the prime is below 2^61. Wider primes and other CPUs retain the scalar
path. Direct evaluation, schoolbook convolution, scalar-build regressions and
36175320 frozen-reference coefficient comparisons pass across all 23 selector
primes and degrees through 65536. Four alternating N32768 trials measured
about 2.2x faster isolated transforms. N32768 selector56 factor checks now pass exact earlier archive identity,
720896 accurate slots and unchanged complete selected files. Observed
raw-prefetch factor time fell from about 34.3 s to 27.5 s, including selector
NTT from about 12.4 s to 6.1 s. N512 all-function refresh also preserves all
three earlier outputs and passes independent phase checks with unchanged
retained files. The full 48/56-bit factor comparison and N32768 all-function
retained replay now pass. Shared coefficient preparation and 22 roots took
2181.40 s (36.4 min), versus 2873.07 s (47.9 min) in the earlier four-worker run. All three independent
refresh/product/normalization checks pass, all three earlier archives match
exactly and all 37 retained files are unchanged. Evaluator and controller exited 0,
with no new keys or inputs. This is an observational shared-root comparison;
standalone paper latency and security remain unqualified. Evidence:
`narrow-dd-ntt-avx2-full-replay-validation.json`,
`narrow-dd-ntt-avx2-full-runtime-profile.json`,
`narrow-dd-ntt-avx2-validation.json`,
`narrow-dd-ntt-avx2-full-factors-validation.json`,
`narrow-dd-ntt-avx2-small-refresh-validation.json`.
Separately, the pre-AVX2 plan-reuse executable passes N32768 head/tail checks:
720896 slots, both earlier archives identical and all selected retained files
unchanged (`narrow-dd-factor-plan-reuse-full-validation.json`).

A geometry-only screen compares the public disjoint-block support with SHIP
Section 5.1's `(4w)^(0.4h)` scalar-product heuristic. The N8192/h31 schedule
falls below SHIP's stated `2^115` scalar-product threshold under this analogy;
N16384 and N32768 exceed it. Transferring that exponent to disjoint ragged
blocks is an assumption, so these results neither prove an attack nor qualify
128-bit security. See
`Parameter-Selection/python/estimates/narrow_ckks_ship_structured_screen.py`
and its JSON result. Moduli, noise, ring/hybrid correlations and circular/KDM
assumptions still require separate analysis.

New basis benchmark checkpoints sample dense ternary secrets with TFHEpp's
randomly seeded key generator. The fixed benchmark seed is used only for public
test inputs. Existing checkpoints can still be replayed without generating keys;
older fixed-seed fixtures provide correctness evidence, not key secrecy evidence.
The creation change passes N512 all-function creation/replay, independent output
checks, identical ciphertext archives and retained-file integrity checks; see
`narrow-dd-random-key-checkpoint-validation.json`.

The basis runner also reports `root_tile_stage_seconds`: exclusive timings for
input/factor-cache preparation, leaf factors and inline alignment, bounded tail
alignment, DD cache preparation, and DD trees with beta encoding. The library's
optional `root_stage_observer` receives one timing record per tile and is omitted
from key archives. These scopes complement the nested factor timings; do not add
the factor timings to the root-tile timings. Shared S2C preparation and per-function
basis reconstruction remain separate report fields. Allocation teardown outside
the tile scopes can appear in the enclosing root-generation residual.
N512 shared and standalone replay preserve the previous ciphertexts and all
retained files; see `narrow-dd-root-stage-timing-validation.json`. The completed
frozen parallel secondary11 full replay below also records these scopes.

Loading a narrow-key archive clears both caller-owned timing observers along
with the derived runtime policies. The file-selector regression checks this
reset; see `narrow-dd-root-observer-load-reset-validation.json`.

Cached DD relinearization switches now distribute independent secondary rows
across CPU workers while retaining the primary-row FMA order for each output.
Each worker keeps one spectrum accumulator. Active-digit recomposition accumulates
each coefficient locally and distributes independent coefficients across workers.
Small rings use the serial path. Retained N32768/Q479 secondary8 and secondary11
switches match the coefficient-domain reference at one and four threads; a separate
Boost integer oracle checks 101,376 recomposed coefficients per thread setting.
N512 all-function replay remains ciphertext-identical with all retained files
unchanged. See `narrow-dd-parallel-switch-validation.json`. These checks do not
establish a full-size refresh speedup; already running full tests retain their
frozen executables.

Prepared automorphism switching also parallelizes independent primary forward
FFTs and coefficient-local secondary recomposition for N >= 8192. Each
coefficient retains its secondary accumulation order; small rings keep the
serial path. The normal-build `ckks_narrow_bootstrapping --prepared-auto-parallel`
control matches coefficient-domain references at N8192/Q85 and Q379 with one
and four threads, including default/secondary16 prepared and original/adaptive
streamed variants (`narrow-dd-prepared-auto-parallel-validation.json`). The
retained N32768 head30/tail4 factor check also preserves exact ciphertext
archives, passes 720,896 decoded slots and leaves complete selected files
unchanged (`narrow-dd-prepared-auto-parallel-full-factors-validation.json`).
Retained full keys are reused. These controls do not establish bounded
tail-alignment or complete-refresh speedup.

The parallel DD component checks additionally pass all seven retained N32768
moduli (479, 432, 385, 338, 291, 244 and 197 bits), with unchanged public-key
files; see `narrow-dd-parallel-switch-full-levels-validation.json`. The optional
secondary11 small replay also preserves the independently validated default
ciphertexts (`narrow-dd-parallel-secondary11-small-validation.json`). Build
`ckks_dd_recombine_exact` through CMake and run it with `OMP_NUM_THREADS=4` for the
101,376-coefficient integer-oracle regression. Its normal-build evidence is
`narrow-dd-recombine-cmake-validation.json`.

The frozen parallel secondary11 full-size evaluator completed with exit status 0,
and all three independent refresh, square and normalization checks pass. All
37 retained files are unchanged, and the controller completed with exit status 0; see
`narrow-dd-parallel-secondary11-full-resumed-validation.json`. It reused verified
frozen artifacts and the completed raw-chunk proof's exact payload baseline after
the original queue stopped before evaluation because baseline metadata was absent.
No full-size keys were generated. These results cover the earlier frozen runner;
the combined current-source replay described below is now fully validated.
The current combined selector56/secondary11 runner is now under full-size
retained-key validation (`narrow-dd-current-combined-full-replay-validation.json`),
which has completed successfully.
Its current-source small ALL replay passes all three independent decoders and
archive identities, and rebuilt retained DD controls pass all seven levels
(`narrow-dd-current-combined-small-replay-validation.json`,
`narrow-dd-current-combined-dd-levels-validation.json`). The controller rejects
wrong small-policy or explicit-width-only factor evidence before evaluation
(`narrow-dd-selector56-replay-controller-controls-validation.json`). It freezes
qualified artifacts, reuses the completed retained-file baseline, and requires
all 37 file hashes after evaluation. This run is terminal with controller exit 0, all three independent phase checks
passed, exact archive identities, and all 37 retained files unchanged. This
qualifies the retained-key functional replay; security and paper performance
remain unvalidated.

The current combined evaluator has now exited 0. Its shared coefficient
preparation took 418.79 s and shared root generation took 2454.28 s. Exclusive
root stages measured 1579.33 s for leaf factors and inline alignment, 395.02 s
for bounded tail alignment, 10.59 s for DD cache preparation and 466.49 s for
DD trees and beta. Nested matrix work took 879.30 s and selector transforms
465.97 s; these are already included in the leaf stage. The observed root time
is 29.1% lower than the preceding frozen run, without a controlled speedup
claim. All three independent decoders pass refresh, square and normalization
checks, and all three ciphertext archives match the earlier validated results.
The post-evaluation checks confirm all 37 retained files are unchanged.
See `narrow-dd-current-combined-full-runtime-profile.json`.
The replay controller accepts `--omp-threads` and records the selected count.
The default remains four. Zero, negative and noninteger worker counts are
rejected before any output directory or evaluator is created. Eight-worker
full-refresh validation has also passed; its end-to-end timing is reported below.
The eight-worker N32768 retained factor control now passes all four exact
48/56-bit trials, 720896 slot checks, both archive identities and the selected
checkpoint/head30/tail4 file hashes. Its configured 56-bit head and tail trials
took 33.57 s and 34.56 s, compared with 40.46 s and 40.81 s in the earlier
four-worker control. This is an observational component comparison; it does not
qualify eight-worker full refresh or establish controlled speedup. Evidence:
`narrow-dd-cpu-eight-worker-full-factors-validation.json`.
Eight-worker DD controls also pass all seven retained modulus levels for
secondary widths 8 and 11 against the coefficient-domain reference. All seven
public key files remain unchanged; qualified executables were reused without
rebuilding or generating keys. Evidence:
`narrow-dd-cpu-eight-worker-dd-levels-validation.json`.
An eight-worker full ALL replay has completed with the same qualified
executable, keys and input. The completed four-worker checkpoint baseline is
reused; all37 files were checked after evaluation and are unchanged. The
controller exited0 and all three independent decoders/archive identities pass
(`narrow-dd-cpu-eight-worker-full-replay-validation.json`). Its evaluator has
exited 0: coefficient preparation took497.34s and shared root generation2420.29s.
The combined2917.63s exceeds the four-worker2873.07s; the four-worker default
therefore remains unchanged. Matrix work improved to744.34s and selector
transforms to296.42s, while selector reads rose to423.88s. Different concurrent
I/O/cache conditions prevent a controlled speedup claim. All three evaluator
accuracy gates pass. All three independent decoders and ciphertext archive
identities also pass;37-file verification confirms unchanged retained data.
See `narrow-dd-cpu-eight-worker-full-runtime-profile.json`.
The factor control exposes `--small-union-wide-prefetch` and
`--full-union-wide-prefetch` to compare synchronous and one-chunk read-ahead
with the configured56-bit policy held fixed. Small retained controls pass all
exact head/tail comparisons and preserve all36 checkpoint files
(`narrow-dd-selector56-prefetch-small-validation.json`). Full-size component
validation passes all four exact head/tail comparisons,720896 slot checks,
archive identities and selected checkpoint/head30/tail4 hashes
(`narrow-dd-selector56-prefetch-full-factors-validation.json`). Synchronous
versus read-ahead head times were24.67/23.64s and tail times31.32/28.83s;
exposed read waits fell to about0.07s. These component observations do not
establish full-refresh speedup. The raw-file refresh path already enables
read-ahead through `supportsRawChunkPrefetch()`; default behavior is unchanged.

The shared 22-root batch took 3462.14 s, with coefficient preparation taking
417.01 s. Exclusive root stages measured 2495.84 s for leaf factors and inline
alignment, 471.50 s for tail alignment, 12.13 s for DD cache preparation and
479.70 s for DD trees and beta; input/cache setup took 2.88 s. Nested factor
measurements include 901.14 s of selector reads, 962.06 s of matrix work and
481.34 s of selector transforms. Factor times are included in the leaf stage and
must not be added again. Only 0.087 s of root time lies outside the root-stage
timers. See `narrow-dd-parallel-secondary11-full-runtime-profile.json`. Later
function reconstruction depends on these shared roots; concurrent runs do not
establish controlled speedup or reproduce standalone paper latency.

Raw file-selector streaming now reads at most one successor chunk on a worker
while matrix work consumes the current prepared chunk. Sources must explicitly
support worker-thread reads; generic providers remain synchronous. The current
raw buffer is released before its successor starts, and read failures propagate
through a joined future. Exact N512 file/chunk regressions and truncated-file
rejection pass (`narrow-dd-raw-selector-prefetch-validation.json`). A retained
N512 checkpoint additionally preserves 3,072 output coefficients exactly across
three inputs, bounds all 105 worker reads to eight selectors, forbids full-group
loads and rejects short reads in both modes. All 36 retained files are unchanged
(`narrow-dd-raw-prefetch-retained-validation.json`); no keys were generated. The exclusive
selector-load timer counts uncovered consumer waits; background read time already
overlapping matrix work is not added again. The retained N32768 head30/tail4 component control also passes with 22 columns
and production chunk size 40: 720,896 decoded slots pass the 23-bit gate, both
archives match earlier qualified results, and the complete checkpoint file plus
both selected selector files are unchanged
(`narrow-dd-raw-prefetch-full-factors-validation.json`). Synchronous/read-ahead
head times were 48.442/50.448 s and tail times were 48.773/47.999 s. These mixed
warm-read observations under concurrent packed reconstruction do not establish
speedup. Full-refresh speedup remains unmeasured, and completed frozen replays do
not include this change.

The factor control also supports `--full-union-wide`, comparing the existing
48-bit selector scratch decomposition with explicit 56-bit limbs on retained
keys. The unchanged compile-time CRT guard permits 56 bits at N32768 and rejects
58 bits. The 56-bit width uses ten limbs instead of eleven, reducing decomposition
rows by 1/11 without changing the ciphertext modulus or key parameters. Both
N512 and N32768 head/tail controls preserve exact ciphertexts; the full control
checks 720,896 decoded slots and complete checkpoint/head30/tail4 file hashes.
See `narrow-dd-wide56-small-factors-validation.json` and
`narrow-dd-wide56-full-factors-validation.json`. Full head times were
48.756/42.290 s (default/56-bit) and tail times were 46.782/41.206 s under concurrent
packed reconstruction. These are component observations, not a controlled or
end-to-end refresh speedup. The default 48/32/16-bit policy remains unchanged.

Selector preparation extracts scratch digits directly from one or two native
limbs, masking the last digit at the ciphertext modulus. This avoids copying and
shifting the complete torus for every digit and prime. The normal-build
`ckks_narrow_bootstrapping --selector-digits` regression checks 149,504 digits
against the previous reduction/shift path, including limb boundaries and storage
bits above the ciphertext modulus
(`narrow-dd-native-selector-digit-validation.json`). The retained N32768 head/tail control also passes: 720,896 decoded slots meet
the accuracy target, ciphertext archives match exactly, and complete selected
checkpoint/head/tail hashes are unchanged
(`narrow-dd-native-selector-digit-full-factors-validation.json`). Observed
timings are close to the previous build; a material speedup is not demonstrated.
Complete-refresh speedup is unmeasured.

A derived narrow schedule may set `selector_matrix_digit_bits = 56` to select
these scratch limbs throughout preparation and factor evaluation. The unchanged
CRT fit guard validates the width, and explicit per-call widths take precedence.
The optional replay command is
`ckks_narrow_basis_replay --small-lean-basis-fbs-replay-secondary11-selector56 all CHECKPOINT_DIRECTORY REPORT.json`
(omit `small-` for N32768). It reuses serialized coefficient keys and derives the
secondary11 DD cache. Complete N512 all-function replay preserves all ciphertext
archives, passes independent refresh/square/normalization and timing checks, and
leaves all 36 retained files unchanged
(`narrow-dd-selector56-secondary11-small-replay-validation.json`). No keys were
generated. Full-size end-to-end replay of this combined policy remains unverified;
the default selector width policy is unchanged.
Build a cache for the configured replay with
`ckks_narrow_selector_ntt_file --selector56 CHECKPOINT_DIRECTORY NEW_CACHE_DIRECTORY`.
The same option works with `--sample-chunks` to create separate partial
sample archives in a new directory. Retained N512 creation/reopening checks
12,185,600 exact spectrum words, and sample controls check 450,560 words
and 13 range/partial-cache rejections. The complete cached selector56/secondary11
replay matches all three independently validated raw-path archives and leaves
all 36 retained files unchanged. Wrong-width caches and duplicate flags reject
(`narrow-dd-selector56-cache-integration-validation.json`). The default cache
command remains supported. These results do not qualify a full-size prepared
cache or establish speedup.

The integrated schedule-configured path also passes the retained N32768
head30/tail4 factor control: all 720,896 decoded slots meet the accuracy gate,
outputs match the default and previous archives exactly, and complete
checkpoint/head/tail hashes are unchanged
(`narrow-dd-selector56-schedule-full-factors-validation.json`). Default/configured
head times were 49.615/41.485 s and tail times were 47.060/41.082 s. These
observations under concurrent packed reconstruction qualify component correctness,
not complete-refresh speedup.

In shared-input benchmark mode, the runner now discards its S2C plaintext plan
after preparing both coefficient halves. The retained checkpoint keeps the
reusable plan. Default and secondary11 N512 replay preserve all saved ciphertexts
and retained files (`narrow-dd-shared-s2c-release-validation.json`). The compiled
N32768 term payload is about 8 GiB; a full-size RSS reduction has not been measured,
and the frozen full tests do not include this benchmark-only release.

The current N32768 all-function root batch performs 2,416,258,842,624 modular
products across 22 columns and rebuilds 1,675,872 selector transforms when streaming
raw selectors. A complete prepared selector cache would require 409.15 GiB against
125.58 GiB of host memory. These are compiled algorithm/storage counts, not timing
predictions; see `narrow-dd-matrix-workload-audit.json`. The first function in
all-function mode includes root generation for every function.

The frozen N32768 raw-chunk/shared-union replay completed with exit status 0.
All three output archives match the earlier full run and pass independent
refresh, multiplication and normalization checks, with no new keys. Its shared
22-root batch took 4886.50 s, including 2215.78 s of selector reads, 974.64 s of
matrix work and 478.79 s of selector NTTs. Shared coefficient preparation took
413.81 s. Function reconstruction then took 3.70 s, 4.00 s and 14.88 s for
Sigmoid, GeLU and Tanh. These reconstruction times depend on the shared root
batch and are not standalone functional-bootstrapping latencies. The enclosing
root time outside the factor-stage timers remains unattributed. Final validation
also verified that all 37 retained files are unchanged; see
`narrow-dd-raw-chunk-streaming-full-output-validation.json` and the live parent
`narrow-dd-raw-chunk-streaming-full-validation.json`. Concurrent observations
do not establish a controlled speedup or reproduce the paper's performance.

The completed N32768 run used an earlier frozen executable and one fresh input
fixture. Independent decoding verified each refresh and a multiplication after
refresh. Its measurements were:

| Function | Maximum refresh error | Refresh time | Paper Table 7 CPU time |
|----------|----------------------:|-------------:|-----------------------:|
| Sigmoid | 5.2631e-9 | 4717.94 s | 2.79 s |
| GeLU | 1.5541e-8 | 5945.51 s | 2.71 s |
| Tanh | 2.0761e-8 | 6774.21 s | 4.51 s |

These are different parameter schedules and hardware, so the timings do not
constitute a controlled paper comparison. All three errors meet the run's
2^-23 refresh target; the actual structured-key security remains unqualified.
See [the initial full-run evidence](narrow-dd-basis-full-initial-all-independent-validation.json).
The Sigmoid retained-key replay completed in 4231.09 s with identical saved
ciphertexts and passing independent phase checks; see
[the Sigmoid replay evidence](narrow-dd-basis-full-sigmoid-retained-replay-validation.json).
The GeLU retained-key replay also passed independent phase checks and produced
identical saved ciphertexts in 4803.12 s; see
[the GeLU replay evidence](narrow-dd-basis-full-gelu-retained-replay-validation.json).
The Tanh retained-key replay also passed independent phase checks and produced
identical saved ciphertexts in 6245.22 s; see
[the Tanh replay evidence](narrow-dd-basis-full-tanh-retained-replay-validation.json).
The complete original retained-key replay passed: all 37 retained files are
unchanged and all three initial/replay archives are identical; see
[the complete replay evidence](narrow-dd-basis-full-retained-replay-validation.json).
The later qualified CPU full-refresh replays are recorded above. In this earlier
replay, loading the same keys and input took 72.24 s, compared with 2937.08 s for initial setup; this reduces setup cost but
does not establish a refresh speedup.

Build the retained-checkpoint runner with the default parameter family:

```sh
cmake -S . -B build -DUSE_BLOCK_BINARY=OFF -DENABLE_TEST=ON -DUSE_BLAKE3=ON
cmake --build build --target ckks_narrow_basis_replay --parallel 1
```

Run these commands from the `TFHEpp` directory. Reuse a compatible checkpoint
with a new report directory for each evaluation:

```sh
mkdir replay-results
OMP_NUM_THREADS=4 build/test/ckks/ckks_narrow_basis_replay \
  --small-lean-basis-fbs-replay all /path/to/checkpoint replay-results/report.json
```

The small profile uses N512. For the retained N32768 profile, use
`--lean-basis-fbs-replay`. Select `sigmoid`, `gelu`, or `tanh` instead of `all`
to evaluate one function. Replay loads the stored keys and encrypted input;
it does not generate new keys. It checks the public parameter profile before
allocating evaluation keys and rejects existing output files. The current normal
CMake build passes exact arithmetic tests, raw/cached N512 replay identity, and
N32768 factor accuracy; see
[the raw chunk validation](narrow-dd-raw-chunk-streaming-validation.json) and
[the full factor validation](narrow-dd-raw-chunk-streaming-full-factors-validation.json)
for those source snapshots. Reports include
accuracy, refresh timing, root reuse, and factor-stage timings. In an `all` run,
the first function generates the combined frequency roots for all three functions
in one selector pass; the later functions reconstruct their outputs from those
roots. These per-function timings are not standalone latency. The
`refresh_timing_scope` field identifies the work included.

The same runner can create a checkpoint once:

```sh
OMP_NUM_THREADS=4 build/test/ckks/ckks_narrow_basis_replay \
  --small-lean-basis-fbs-checkpoint all new-checkpoint initial-report.json
```

This creates keys and an encrypted benchmark input, saves them privately, and
evaluates the selected functions. Use `--lean-basis-fbs-checkpoint` for N32768.
The checkpoint directory and report must be new. An existing checkpoint is
rejected with a reminder to use replay, so repeats cannot silently regenerate
its keys. The report goes outside the new checkpoint directory. The legacy
`ckks_narrow_benchmark` commands remain available.
The full profile retains about 167 GiB of raw selectors and is substantially
slower than the paper on this implementation.

Prepared-selector disk caching is optional: append a compatible cache directory
to the replay command. It removes selector NTT transforms but adds cache reads
and integrity checks. Full-degree component trials on shared storage found cold
cache reads slower than fresh transforms, so this is not a demonstrated speedup.

An optional `--[small-]lean-basis-fbs-replay-secondary11` mode uses the same
checkpoint with 11-bit secondary digits in the derived DD product-tree FFT
cache. Serialized key digits remain unchanged. It accepts the same function,
checkpoint, report, and optional selector-cache arguments as normal replay.
N512 all-function archives match default replay exactly and pass independent
refresh, multiplication, and normalization checks; see
[the replay validation](narrow-dd-basis-secondary11-small-validation.json).
The retained N32768 public DD rows also pass component checks at all seven
consumed levels; see
[the retained-key checks](narrow-dd-secondary11-retained-key-validation.json).
Only the product-tree cache policy changes in this mode. Complete N32768
refresh validation and a measured refresh speedup remain pending; normal replay
retains its existing cache policy.

The auxiliary retained-fixture diagnostic runner also accepts
`--replay-combined-primary8 STATE SELECTORS NEW_DIRECTORY`. This enables both
projected primary8 conjugation and projected primary8 DD relinearization using
the existing coefficient keys. N512 default replay remains ciphertext-identical;
the combined mode passes independent refresh and post-refresh multiplication
checks with all 36 retained files unchanged. See
[the combined-mode validation](narrow-dd-combined-primary8-small-refresh-validation.json).
This optional mode has not been qualified at N16384 or N32768 and does not
change the normal basis runner or the frozen full-size validation runs.

The sections below also describe the separate dense CKKS implementation.

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

## CKKS / TFHE Scheme Switching

For the **4-bit TFHE → CKKS multiplication → TFHE** objective, use the
specialized [integer multiplication path](Integer4Multiplication.md). Its
components pass the 128-bit classical screen and its modeled failure bound
is below 2^-40. It uses constant CKKS plaintexts and exact NTT arithmetic,
so it needs neither the dense slot transforms nor EvalMod described below.

**General-path security calibration:** the measured full roundtrip below fails the
128-bit target (111.87-bit CKKS screen and 61.59-bit compact-key screen).
See [SchemeSwitchingCalibration.md](SchemeSwitchingCalibration.md) for the
Parameter-Selection results, component targets, and the remaining precision
and switching-key changes. Its memory and correctness results describe the
functional prototype only.

`include/ckks/scheme_switching.hpp` provides CKKS-to-TFHE slot decoding and
extraction, and TFHE-to-CKKS encrypted partial decryption followed by periodic
EvalMod. Evaluation uses public evaluation keys; secret keys are required only
for setup and result verification.

The available path is:

```text
CKKS coefficient ciphertext -> sample extraction -> modulus lift
    -> identity key switch -> TFHE TLWE
```

Generate the cross-scheme key with `CKKSToTFHEKeyGen`, then invoke
`CKKSCoeffToTFHE` for a decoded coefficient.  `InputLogQ` bits must be covered
by the identity-key-switch decomposition (and therefore fit in the TFHE torus;
reduce to an OpenFHE-style intermediate CKKS modulus first when necessary),
and the output has
`CKKSToTFHEOutputScaleBits` fractional bits.  The API statically rejects a
target TFHE torus that would leave no message/noise gap. The current
cross-key implementation accepts a ternary CKKS secret and a conventional
32- or 64-bit TFHE torus; it also supports multi-limb CKKS ciphertexts.

For the compact form used by OpenFHE, generate
`CKKSToTFHERingSwitchKey` with `CKKSToTFHERingSwitchKeyGen`. It switches once
to a CKKS ring secret containing the TFHE secret followed by zeros, then uses
sample extraction to produce TFHE TLWEs. `CKKSCoeffToTFHEViaRingSwitch` is the
single-coefficient entry point; `CKKSCoeffsToTFHEViaRingSwitch` amortizes the
ring switch over a consecutive coefficient batch. The compact key is seeded
per CKKS gadget row, avoiding the direct cross-key's TLWE ciphertext for every
source-secret coefficient. This compact route requires a multi-limb CKKS
torus. Its scheme-switch key rows and products use an exact chunked
negacyclic NTT rather than the ordinary floating CKKS FFT, preserving the
low-modulus bits needed by TFHE extraction. The exact compact implementation
currently supports power-of-two rings through degree `2^15` and Bbar digits
through 20 bits. `ckksEncrypt` uses the same exact product for these
multi-limb parameters, so freshly encrypted coefficient ciphertexts retain
the bits consumed by the compact switch.
For non-consecutive or sparse layouts, use
`CKKSCoeffIndicesToTFHEViaRingSwitch`; its output order follows the supplied
coefficient indices. `CKKSSlotsToTFHEViaRingSwitchAtIndices` adds packed-slot
decoding before the same one-switch batched extraction.
`CKKSOpenFHEExtractionIndices` derives OpenFHE's standard
`RingDim / (2 * numSlots)` stride, while
`CKKSSlotsToTFHEOpenFHEPacked` and
`CKKSSlotsToTFHEViaRingSwitchOpenFHEPacked` apply it directly for the direct
and compact paths.

For directly packed CKKS inputs, `CKKSSlotsToCoeffs` applies the same
SlotToCoeff linear map and `CKKSSlotsToTFHE` composes that map with extraction.
Both require the input- and output-level Galois keys for the chosen BSGS
layout. `CKKSSlotsToCoeffs` now generates and evaluates diagonals on demand;
it does not construct `CKKSRealLinearTransformPlan`. The binary and hybrid
rotation-key collectors also derive their schedules without encoding any
plaintexts. `CKKSSlotsToCoeffsRotationKeyIndices` exposes precisely the required
binary rotation indices, and `CKKSSlotsToCoeffsSparseGaloisKeyGen` creates the
matching sparse input/output key pair. It accepts either unseeded or
`CKKSSeededSparseGaloisKey` outputs; seeded rows regenerate masks on use.
`CKKSSlotsToCoeffsHybridRotationKeyIndices` and
`CKKSSlotsToCoeffsHybridSparseGaloisKeyGen` additionally choose direct versus
binary keys independently for the BSGS baby and giant rotations. The materialized full-Galois transform builders remain available as reference
implementations. The streamed evaluator generates one giant-step group at a time while
reusing the materialized evaluator's fused FFT products and baby-step cache.
`CKKSSlotsToTFHEViaRingSwitch` composes this same SlotToCoeff transform with
the compact ring switch; generate its ring key at the decoded output level
`InputLogQ - PlainLogDelta`.
To match OpenFHE's distinct extraction modulus `Q'`, use
`CKKSSlotsToTFHEAtLevel` or `CKKSSlotsToTFHEViaRingSwitchAtLevel` and select
`ExtractionLogQ <= InputLogQ - PlainLogDelta`. The output TFHE scale changes
accordingly through `CKKSToTFHEOutputScaleBits`; the compact ring key must be
generated at `ExtractionLogQ`. For OpenFHE's usual combination of `Q'` and a
sparse packed-slot stride, use `CKKSSlotsToTFHEAtLevelOpenFHEPacked` or
`CKKSSlotsToTFHEViaRingSwitchAtLevelOpenFHEPacked`.

The reverse key is `TFHEToCKKSSwitchKey`, generated by
`TFHEToCKKSKeyGen`. `TFHEToCKKSPhase` produces the coefficient-zero CKKS
encryption of the unwrapped normalized LWE phase `b/Q - a·s/Q`.

For a batch of TLWEs, `TFHEToCKKSPackedSwitchKey`,
`TFHEToCKKSPackedKeyGen`, and `TFHEToCKKSPackedPhase` instead encrypt the
TFHE secret in real CKKS slots and apply the public-mask linear transform.
The key repeats the secret in real-slot blocks; phase `i` is emitted at slot
`i * lwe_dimension`, allowing `floor(N / (2 * lwe_dimension))` TLWEs per
CKKS ciphertext. The packed routine is a correctness/reference path: it uses
the direct-diagonal (`N/2` BSGS-step) layout, so callers should provision the
small rotations used by their masks or replace it with a fused sparse schedule
for production parameters. As in OpenFHE, either reverse output is only
meaningful modulo one until the following periodic modular-reduction stage.
`TFHEToCKKSPackedPhaseRotationKeyIndices` and
`TFHEToCKKSPackedPhaseSparseGaloisKeyGen` provide the equivalent sparse-key
provisioning for this transform; the replicated-secret layout needs only the
binary rotations that compose offsets below `lwe_dimension`.

For OpenFHE-compatible dense output packing, use
`TFHEToCKKSSlotPackedSwitchKey`, `TFHEToCKKSSlotPackedKeyGen`, and
`TFHEToCKKSSlotPackedPhase`. This stores the encrypted TFHE secret once and
uses a streamed rectangular BSGS transform to place up to `N/2` unwrapped
phases in consecutive CKKS slots. Public LWE masks are read one diagonal at a
time, without a dense diagonal matrix or an encoded plan. `TFHEToCKKSSlotPackedPhaseRotationKeyIndices` and
`TFHEToCKKSSlotPackedPhaseSparseGaloisKeyGen` derive the corresponding sparse
input/output keys for a selected batch size and BSGS step, including a seeded
key overload. Feed the result to
the same bounded-cosine EvalMod stage to obtain canonical messages, or use
`TFHEToCKKSSlotPackedEvalModKeyGen` and
`TFHEToCKKSSlotPackedEvalMod` for the composed dense reverse switch.
`TFHEToCKKSSlotPackedEvalModAffine` additionally performs the final public
`output_scale * phase + output_bias` map (for example, the requested CKKS
plaintext range), consuming a caller-selected plaintext multiplication level.
`TFHEToCKKSSlotPackedEvalModOpenFHE` is the corresponding convenience wrapper
for OpenFHE's `p`, `pmin`, and `pmax` post-processing convention.
For the sparse-output convention used by OpenFHE's `EvalFHEWtoCKKS`,
`TFHEToCKKSSlotPackedEvalModToSlots` repeats a batch across a logical slot
count that is a power-of-two multiple of the batch. Generate its final-level
rotation keys with `CKKSSchemeSwitchReplicationSparseGaloisKeyGen`; full
packing needs none.

`TFHEToCKKSPackedEvalModKeyGen` and `TFHEToCKKSPackedEvalMod` provide that
following stage: they normalize the encrypted partial-decryption phase by a
caller-declared bound and run the existing bounded-cosine EvalMod circuit. The
bound must exceed `|b/Q-a·s/Q|`; `lwe_dimension + 1` is conservative for a
ternary TFHE secret. The resulting block-leading slots contain canonical
approximate CKKS messages. The EvalMod degree, double-angle count, precision,
and ciphertext modulus are template parameters, so production callers must
select a chain with enough remaining levels for their accuracy target.

### Bounded-memory execution and sparse inputs

The old forward key collector built two dense encoded transforms merely to
learn which rotations were needed. With `lvl6param` (`N=32768`, 112-byte
multi-limb coefficients), those encoded plaintexts alone take
`2 * (N/2) * N * sizeof(T) = 112 GiB`, before raw diagonals, keys, or scratch.
Collecting the schedule now takes O(N log N) time and O(log N) binary-key
metadata, without plaintext allocation.

`CKKSLinearTransformBSGSStreamed` retains one encoded giant-step group and
the used baby-step ciphertexts (or their digit FFTs), instead of O(N^2)
plaintext storage. Its provider fills one unrotated diagonal for each requested offset. Offsets may be
negative or unsorted, but normalized duplicates are rejected. The resident
baby ciphertexts cost at most `bsgs_step * sizeof(TRLWE<P>)`. On the fused
multi-limb backend, the baby FFT cache costs up to
`bsgs_step * (P::k+1) * P::l̅ * sizeof(PolynomialInFD<P>)`; the encoded group
has at most `bsgs_step` plaintexts. Rotation keys, EvalMod keys, and arithmetic
scratch are additional. Full logical packing still requires quadratic
diagonal arithmetic; streaming fixes storage, not
that computational cost.

For OpenFHE-style sparse **input** packing, repeat `logical_slots` values
across all N/2 physical CKKS slots, then call
`CKKSSparseSlotsToTFHEViaRingSwitch`. It decodes using only `logical_slots`
diagonals per real-linear branch, switches the ring key once, and extracts
at `i * N/(2*logical_slots)`. Pass the same final `logical_slots` argument to
`CKKSSlotsToCoeffsSparseGaloisKeyGen` (after its noise argument). The lower-level
`CKKSSlotsToCoeffs` and rotation collectors also accept this optional argument;
its default is full physical packing. The existing `...OpenFHEPacked`
wrappers keep their original selected-coefficient semantics and do not assert
that their inputs have this periodic layout.

From the repository root, build and run the focused checks with a Linux
address-space cap (bytes):

```sh
cmake -S . -B build -DENABLE_TEST=ON -DUSE_BLOCK_BINARY=OFF
cmake --build build --target ckks_streamed_transform ckks_to_tfhe_scheme_switch ckks_to_tfhe_full_smoke ckks_tfhe_roundtrip ckks_tfhe_roundtrip_quick -j2
OMP_NUM_THREADS=4 prlimit --as=134217728 --core=0 -- build/test/ckks/ckks_tfhe_roundtrip --schedule-only
OMP_NUM_THREADS=4 prlimit --as=8589934592 --core=0 -- build/test/ckks/ckks_streamed_transform
OMP_NUM_THREADS=4 /usr/bin/time -v prlimit --as=1073741824 --core=0 -- build/test/ckks/ckks_tfhe_roundtrip_quick
OMP_NUM_THREADS=8 /usr/bin/time -v prlimit --as=34359738368 --core=0 -- build/test/ckks/ckks_tfhe_roundtrip --full
```

`ckks_streamed_transform` checks complex matrix-vector products against a
plaintext oracle, non-dividing BSGS steps, normalized-offset validation,
rotation schedules against the materialized reference, and periodic sparse
inputs including their real and imaginary coefficient layout.
`ckks_tfhe_roundtrip --schedule-only` exercises the formerly memory-heavy
full lvl6 schedule under 128 MiB without generating keys. This schedule check
is also the default; the multi-GiB encrypted run requires `--full`.
`ckks_tfhe_roundtrip_quick` checks the entire encrypted path at N=128,
with an experimental eight-dimensional TFHE secret and eight real values.
It takes about one second and 15 MiB here. Its `--full-precision` mode uses
the large test's modulus/scale chain in the small ring (about six seconds).
Both recover values in [-1,1] with roughly 0.0016 maximum error, dominated by
the sine approximation used for modular reduction.

The opt-in full test uses the precision-oriented `FullCKKS` profile in
`test/ckks/scheme_switch_test_params.hpp`: N=32768, a 16-limb torus, reverse
logQ=1000, 58-bit ciphertext/coefficient scales, and a 60-bit linear plaintext
scale. It uses 760-dimensional lvlhalf TFHE and eight repeated real slots.
It checks the forward TLWE phases, encrypted partial decryption, periodic
EvalMod, and an encrypted integer post-scale. Keys are independently
randomized and CKKS uses centered Gaussian noise with a 3.2-unit floor at
every active level. These are functional configurations, not calibrated
128-bit-security claims. The full roundtrip can take tens of minutes;
large multi-limb key switches dominate the linear stage. Set
`TFHEPP_CKKS_TRACE_LINEAR=1` and `TFHEPP_CKKS_TRACE_EVALMOD=1` for stage timings.

On this Ryzen 9 9950X3D host with eight OpenMP threads, the full encrypted
roundtrip passed under the 32 GiB address-space cap:

| Measurement | Result |
|-------------|-------:|
| Peak resident memory | 16.61 GiB (17,415,804 KiB) |
| Total time, including key generation | 46 min 30.94 s |
| Forward maximum error, original message units | 7.90e-7 |
| Reverse normalized phase maximum complex error | 6.44e-9 |
| Final roundtrip maximum complex error | 0.001902 |

The reverse linear transform took 1768.06 seconds and EvalMod took 604.10
seconds. This is one randomized functional run, with no swapping or OOM;
it establishes that this configuration fits this machine, rather than a
memory guarantee for arbitrary parameters or packing sizes.

The noise choice matters: the four-bit signed uniform test distribution
[-8,7] has mean -0.5, whose bias accumulates through long key-switch schedules.
The optional `ckks_seeded_key_precision --rotation-noise` diagnostic compares
zero, bounded-uniform, and Gaussian noise in one full-ring rotation. Its
default mode checks exact zero-noise seeded-row decryption at N=32768.
`ckks_streamed_transform --full-arithmetic` separately checks the fused
16-limb transform against a plaintext oracle, using transparent ciphertexts
to isolate arithmetic from key noise; it is not an encrypted roundtrip.


This is not yet a drop-in equivalent to OpenFHE's tuned production API. The
remaining work is to add:

1. a schedule-bound fused packed transform and an accelerated RNS/hybrid
   specialization beyond the current exact-NTT compact-cross-key limit;
2. application-specific security, noise, and precision calibration for the
   composed linear, compact cross-key, and EvalMod stages. The functional
   roundtrip test does not constitute a security estimate.

The regression `test/ckks/ckks_to_tfhe_scheme_switch.cpp` covers packed
CKKS-slot decoding/extraction into TFHE through both direct and compact ring
keys, coefficient-domain reverse phase construction, and packed TFHE-to-CKKS
partial decryption followed by EvalMod, with independent CKKS and TFHE secret
keys. `test/ckks/ckks_to_tfhe_full_smoke.cpp` separately exercises public
zero-noise CKKS encryption and the exact compact ring switch at TFHEpp's lvl5
(`N=2^14`) CKKS degree and lvl2 TFHE dimension, including 88-bit-scale
positive and negative multi-limb coefficient encoding. It also checks the
reverse encrypted-secret phase construction with a nonzero TFHE mask entry.

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

## Narrow Double Decomposition refresh

The experimental port of `refs/2026-2119.pdf` exposes narrow blind rotation,
S2C and functional refresh through `ckks/narrow_bootstrapping.hpp` and
`ckks/narrow_refresh.hpp`. `ckks/narrow_frequency_basis.hpp` reconstructs
missing positive frequencies with one additional DD multiplication level,
then applies fixed-point public coefficients. Evaluation uses public keys;
secret keys are used for generation and accuracy checks.

The benchmark rejects nonfinite slot-error and amplitude metrics before
updating their maxima. This prevents a NaN from being hidden by a prior
finite maximum. `--validation-finiteness` checks NaN and both infinities in
both operand positions, including the long-double switching diagnostics.
This validation change does not establish the cause of the earlier full-size
sigmoid failure.
A separate O1 mixed-source UBSan diagnostic passes the small N=512 plain37
sigmoid check, including fixed 23-bit and square gates, with no sanitizer
diagnostics. Full N=16384 sigmoid used the same instrumented benchmark/header
templates; existing libraries and FFT assembly are outside this instrumentation scope (`narrow-dd-ubsan-diagnostic.json`). The full-ring
O1 sigmoid run subsequently exited 1 after 2:44:05 with target/evaluation error
1.2028444e20 and no UBSan diagnostic. Fresh input, S2C and first-factor errors
remain small. Selector cleanup completed. The independent Python importer
also rejects the report. This recurrence makes an O3-only explanation
insufficient; it does not identify the root cause. Existing static libraries
and FFT assembly are outside the sanitizer scope, which also does not detect
data races. Failed random keys were not saved. New plain37 diagnostics can
retain the complete keys, secrets, input and selectors before evaluation:

```bash
build/test/ckks/ckks_narrow_benchmark --ll14-plain37-fbs-checkpoint sigmoid /artifact/new-checkpoint report.json
build/test/ckks/ckks_narrow_benchmark --ll14-plain37-fbs-replay sigmoid /artifact/new-checkpoint replay.json
build/test/ckks/ckks_narrow_benchmark --ll14-plain37-fbs-replay-trace sigmoid /artifact/new-checkpoint trace.json
```

Prefix these modes with `--small-` for N=512 checks. The checkpoint directory
contains private keys and is created with owner-only permissions. It must be
new for generation; a completed archive is published before refresh starts.
Replay reuses native coefficient keys, input and selectors, and rebuilds derived
caches. Reports have checkpoint flags and an adjacent `.ciphertext` file for
exact refresh output comparisons; post-refresh product keys are generated
fresh. The small checkpoint and both replay modes pass and produce identical
refresh ciphertext bytes. The retained selector reader checks every group
header and fixed length, retains files on destruction, and preserves normal
owned-cache cleanup. Its standalone UBSan test passes. Evidence is in
`narrow-dd-checkpoint-replay-validation.json`. The full N=16384 checkpoint
run completed successfully, and its traced replay produced the same refresh
ciphertext byte for byte (`narrow-dd-full-checkpoint-completion-validation.json`
and `narrow-dd-full-traced-replay-completion-validation.json`). Two further
retained-key trials passed, including
`narrow-dd-ll14-plain37-checkpoint-reproduction-2.json`. None reproduced a full
failure. A separate original-style full sigmoid trial uses automatic
post-evaluation failure capture and completed successfully without a
pre-evaluation checkpoint. Its maximum refresh error is 4.30e-9 and square
error is 9.41e-9; normal exit, selector cleanup and independent report
validation passed (`narrow-dd-original-style-capture-completion-validation.json`).
The original failure remains unresolved. The failed report is preserved in
`narrow-dd-ll14-plain37-sigmoid-ubsan.json`; the separate Python validation
artifact preserves this failure without overwriting either earlier aggregate.

Packed selectors can also be reconstructed directly into a new selector
cache directory with `CKKSNarrowReconstructSelectorsToFile`. Extraction emits
the same ordered public terms used by the in-memory adapter. A temporary
portable record file accumulates them with one selector in memory; final
conversion holds one selector group. The temporary file is removed before
the completed source is published. Existing directories are rejected, and
an extraction/write failure leaves the previous key source intact. Temporary
records and final groups together require roughly twice the expanded
selector payload on disk during finalization.

`--small-high-packed-bounded-file NEW_DIRECTORY REPORT` exercises packed-key
archive round-trip, two prepared substitution depths, disk reconstruction,
selector decryption checks, and refresh. `--ll14-packed-bounded-file` selects
N=16384 with the same schedule. The N=512 encrypted run passes the fixed
16-bit refresh and square gates: maximum refresh error 1.13e-6, square error
3.99e-7, selector reconstruction error 1.84e-9, and peak RSS 394,048 KiB
(`narrow-dd-small-high-packed-bounded-file.json`); exit and selector cleanup
were verified. Both modes use compact DD/S2C storage and staged
FFT evaluation with chunks of 40 selectors and bounded tail rotations.
This is the standard narrow refresh profile, rather than a full-domain
functional profile. The bounded N=16384 one-shot run completed with exit0:
refresh maximum error3.3470660980609124e-6 and square maximum error
1.093324794545107e-6 pass the fixed2^-16 gates. Strict report and owned-cache
cleanup validation is in `narrow-dd-full-packed-completion-validation.json`.
That run discarded its keys and output, so independent saved-phase checking
requires the new retained checkpoint run. The original launch used a 64GiB
virtual limit and a compiled conservative envelope of47.60GiB; its measured
peak RSS was17,563,328KiB. Exact payloads and original launch conditions are recorded in
`narrow-dd-packed-bounded-storage.json` and
`narrow-dd-packed-bounded-memory-budget.json`. The focused regression command is
`build/test/ckks/ckks_narrow_bootstrapping --packed-selectors`. It passes
exact equality for every head/tail selector with raw and prepared substitution
keys, including final-block padding, source preservation on directory rejection,
and owned-file cleanup (`narrow-dd-packed-file-validation.json`).

Packed extraction reuses one pair of child ciphertexts per DFS depth, a
shared scratch pair, and leaf/product buffers. This reduces payload-buffer
allocations from `16*N-8` to `nbit+3` while retaining the ordered sink terms.
The sink must consume or copy each term during its callback. The regression
options `--packed-workspace-fixture-write NEW_DIRECTORY` and
`--packed-workspace-fixture-check DIRECTORY` save a public packed key and its
ordered terms, then compare a later implementation against those ciphertexts
with zero and two cached substitution depths. A fixture captured before the
reuse and plaintext FFT changes matches all 4,096 terms exactly at N=512;
the independent packed-selector file/in-memory regression also passes.

Plaintext multiplication omits torus rows above the active coefficient modulus
and all-zero unsigned plaintext digit columns. The default
`PolyMulTorusByUnsignedBits` template retains full native-width semantics;
an explicit `ActiveBits` promises equality modulo `2^ActiveBits`. CKKS callers
reduce at their active level before rounded rescaling. The standalone
`ckks_plain_mul_exact` test uses arbitrary-precision direct negacyclic
convolution at N=16384, with positive and negative 64-bit magnitudes, wraparound,
native high-bit data, carry/sign boundaries, and 81-, 413-, and 448-bit moduli.
Every coefficient matches, including the default unsigned path. Replaying a
retained N=512 sigmoid state also passes its accuracy gates and produces the
same refresh ciphertext byte for byte. Evidence is in
`narrow-dd-packed-workspace-validation.json`. These coefficient/workspace
checks are separate from the completed full-size packed report validation and
the ongoing retained full-size checkpoint/independent-phase run.

The benchmark can evaluate sigmoid, GeLU and windowed tanh with one public
setup and one encrypted input. The selector directory must be new; its owner
removes generated files on normal completion. Multi-limb selector groups use
bulk portable word I/O, with the same bytes as the legacy records and
per-word endianness conversion. For a small validation run:

```bash
OMP_NUM_THREADS=4 build/test/ckks/ckks_narrow_benchmark \
  --small-lean-basis-fbs-file all /tmp/narrow-basis-selectors \
  docs/narrow-dd-small-lean-basis-all.json
```

This profile uses h=31, Gaussian sigma=3.2, component domain [-4,4], primary
digits of 18 bits, secondary digits of 8 bits and 576-bit storage. At N=512,
all three functions pass the unchanged 23-bit target and evaluation gates,
retained-square checks and dense scale/modulus normalization. Normalizing a
square does not establish that it lies in the domain for another nonlinear
refresh. The baseline uses `--small-basis-fbs-file all` and 640-bit storage.

`--lean-basis-storage REPORT` reports compiled N=32768 payloads without
allocating the full keys. `--lean-basis-fbs-file all NEW_DIRECTORY REPORT`
runs that full profile. The conservative process envelope is 41.6 GiB,
including a 10 GiB reserve, with a proposed 48 GiB virtual-memory limit;
file selectors need another 167.4 GiB on disk. These are storage estimates,
not measured full-size peaks. Use a sufficiently large artifact volume.
The full N=32768 run passed the fixed 23-bit accuracy gates for sigmoid,
GeLU, and tanh, including retained-square and dense normalization checks
(`narrow-dd-lean-basis-all.json`). It exited successfully and removed its
selector files; measured peak RSS was 29,744,128 KiB. Generic dense LWE
screens do not qualify the structured sparse secret or circular/KDM
assumptions. Current evidence and outstanding work are tracked in
[2026-2119-status.txt](2026-2119-status.txt).

The HT diagnostics `--small-ht14-lean-staged-file` and
`--small-ht15-lean-staged-file` use wider DD digits and check nine retained
products with fresh encrypted unit operands. Both pass at N=512; the full
commands omit `small-`. CMake disables GNU IPA-ICF for this benchmark:
otherwise identical GCC 13 builds reproducibly fail its wider-DD S2C path
with that optimization enabled. The setting applies during compilation and
LTO linking. The failed controls are retained in the status record.
For full HT runs, set `ulimit -s 65536` before launching the benchmark;
large template call frames can exceed the default 8 MiB stack. Fixed-size
benchmark buffers and the factor polynomial now use heap storage, and
both small HT regressions pass. The full HT14 retry uses this 64 MiB stack
allowance within its unchanged 24 GiB process limit. Full N=16384 HT14
now passes the fixed 4-bit refresh, square and nine retained-product gates
(`narrow-dd-ht14-lean.json`), exits successfully and cleans up its selectors.
Peak RSS was 7,084,060 KiB. Full HT15 reached 19.77-bit refresh precision
and completed six retained-product checks, then crashed with signal 11
before its seventh check, without a final report or selector cleanup.
Compiled nested checker frames exceed the 64 MiB stack limit
(`narrow-dd-ht15-retained-stack-frames.json`). The benchmark now performs
each product check in a separate call that returns before the next level;
both small HT profiles pass all nine products after this change. Compiled
full-size product helpers now use only 256--320 bytes for their own frames;
deeper arithmetic frames are released before the next product. The full
N=32768 HT15 retry now passes the fixed 16-bit fresh-input, refresh, square,
and all nine retained-product gates (`narrow-dd-ht15-lean.json`). Its refresh
maximum error is 1.1599258904e-6 (19.72 bits); the worst retained-product
error is 1.1532248680e-6. The chain ends at Q42/Delta39, with exact raw/prepared
DD agreement enforced at every level. The process exits successfully and
cleans up its selectors; peak RSS is 46,809,684 KiB. The original crash
artifacts remain preserved. Structured-secret security qualification and
GPU performance remain incomplete.

`--small-ll14-plain37-fbs-file FUNCTION NEW_DIRECTORY REPORT` tests a further
LL14 candidate after the full-size GeLU baseline narrowly missed the 23-bit
gate. It uses selector385/rootDelta48/plain37, keeps outputQ108, and needs
an optional packed-key modulus of 436 bits, within the 437-bit LL14 budget.
The full command is `--ll14-plain37-fbs-file`. Full N=16384 GeLU passes
the fixed 23-bit target and evaluation gates, with errors of 1.99e-8 and
2.04e-8, and its propagated retained-square gates
(`narrow-dd-ll14-plain37-gelu.json`). The run exited successfully and
removed its selectors; peak RSS was 14,389,284 KiB. The full sigmoid run
failed with error 1.79e20 (`narrow-dd-ll14-plain37-sigmoid.json`), despite
small fresh-input, S2C and first-factor errors. Small controls pass, but the
full failure remains unresolved. A separate full sigmoid reproduction uses
`--ll14-plain37-fbs-trace-file` to decrypt and compare each weighted root;
the traced reproduction passes all 23-bit accuracy and propagated-square
gates, with target error 4.39e-9, and all 24 roots have errors below 1.9e-9
(`narrow-dd-ll14-plain37-sigmoid-root-trace.json`). This does not explain the
earlier failure. A full untraced repetition also passes the fixed 23-bit
gates and retained square, with target error 4.66e-9 and evaluation error
4.00e-9 (`narrow-dd-ll14-plain37-sigmoid-untraced-repeat.json`); it exited
successfully and cleaned its selectors, with peak RSS 14,387,848 KiB. The
original failure remains preserved and unexplained. Trace diagnostics add
work to the measured refresh time. Full tanh now passes the fixed 23-bit
refresh, evaluation, complex uniform-bound and retained-square gates:
maximum error 2.0260509755e-8 (25.56 bits), evaluation error 7.08e-9, and
square error 5.70e-8 (`narrow-dd-ll14-plain37-tanh.json`). It exited
successfully and removed its selector directory; peak RSS was 14,387,584 KiB.
The independent aggregate using the successful sigmoid repeat, GeLU and
tanh passes all three functions. The aggregate retaining the original
sigmoid report fails; that failure remains unexplained. Both are preserved
in Parameter-Selection as `narrow_ckks_ll14_plain37_{repeat,original}_validation.json`.
Accuracy does not qualify structured-secret security or GPU performance.

`--small-ll13-fbs-file FUNCTION NEW_DIRECTORY REPORT` tests a low precision
candidate on the full componentwise domain [-4,4]. Its fixed target precision
is 9.04 bits for sigmoid, 7.17 for GeLU, and 7.89 for tanh, from the paper's
low precision results. It uses selector217/rootDelta27/plain24 and outputQ58.
It normalizes by three bits before a sigmoid/tanh square and six before a
GeLU square, leaving productQ31. This profile uses expanded file-backed
selectors: its optional packed setup would require 254 bits at N=8192,
exceeding the 218-bit evaluation-key budget. All three small encrypted tests
pass, including retained squares. Sigmoid also passes at N=8192
(9.41 bits versus the fixed 9.04-bit gate),
including its retained square. GeLU passes at N=8192 with8.18 bits against
its fixed7.17-bit gate, also with a retained square. Full tanh passes at
8.45 bits against its fixed 7.89-bit gate, including its retained square.
All three full-ring accuracy checks pass; structured-secret security remains
unqualified.
The full command is `--ll13-fbs-file`. Parameter-Selection records
both the used key budget and the unsupported packed budget separately.
The standard `--storage REPORT` command includes this profile. The sibling
candidate script accepts `--storage REPORT --memory-report OUTPUT` to
reproduce its conservative12.9 GiB process envelope with a16 GiB limit.

CUDA is outside the requested scope. The experimental CUDA sources, headers,
tests and build targets have been removed. The validation JSON files and
chronological status entries retain historical results for removed code.
The current CPU scope, source hashes, full-profile evidence and remaining gates
are recorded in `narrow-dd-current-cpu-scope-audit.json`.

The full plain37 sigmoid retained-checkpoint diagnostic completed successfully
at N=16384/h=31: maximum refresh error 3.95e-9 and post-refresh product error
9.78e-9, with exit zero and peak RSS 14387808 KiB. Its new keys/input and public
selectors remain in the private checkpoint directory. An independent parameter
importer accepts this report alongside the existing GeLU/Tanh reports. Original
and O1 sigmoid failures remain unresolved. The completed traced replay is
described below. Evidence for the retained run is in
`narrow-dd-full-checkpoint-completion-validation.json`.

The full N=16384 same-state traced sigmoid replay completed with exit zero,
24 finite root diagnostics and a byte-identical public refresh ciphertext to
the retained untraced PASS run. Refresh maximum error is 3.95e-9; an independent
importer accepts this sigmoid report with the existing GeLU/Tanh reports. Earlier
unsaved failure cases remain unresolved. A new untraced checkpoint run is
seeking a retained failing case after a fresh joint memory/disk audit. Evidence
is in `narrow-dd-full-traced-replay-completion-validation.json` and
`narrow-dd-checkpoint-reproduction-memory-budget.json`.

The CPU regression completed successfully after CUDA removal, including both
DD tree orders, compact/prepared keys, every h=31 packed selector, full complex
refresh and nonlinear/basis evaluation. The second retained full sigmoid trial
also passed: refresh error 3.91e-9, RMS error 1.40e-9 and post-refresh product
error 9.88e-9. Independent prefix inspection reproduces its refresh metrics
exactly. The run completed normally in 1:25:41 with peak RSS 14387444 KiB;
evidence is in `narrow-dd-second-reproduction-completion-validation.json`.
Earlier unsaved failures remain unresolved, and security remains unqualified.

`ckks_narrow_checkpoint_inspect CHECKPOINT.bin REPORT.json.ciphertext` reads
only the saved metadata/secrets/input prefix and checks the refreshed ciphertext
against the saved target slots. It supports the N512/N16384 plain37 profile,
rejects malformed profile/sidecar data and nonfinite metrics, and returns 1
when the fixed 23-bit refresh gate fails. It does not validate the evaluation-key
suffix or retained multiplication. Three saved cases reproduce the benchmark's
fresh, refresh and RMS errors exactly without loading evaluation keys; full-ring
inspection uses about 18 MB of memory. Evidence is in
`narrow-dd-checkpoint-prefix-validation.json`. The second trial budget is in
`narrow-dd-checkpoint-reproduction-2-memory-budget.json`.

`ckks_narrow_checkpoint_inspect --phase-spectrum CHECKPOINT.bin OUTPUT.ciphertext`
checks full N16384 plain37 sigmoid output phases with exact integer NTT secret
multiplication. It reports normalized coefficient moments, deviation from a
uniform coefficient CDF, and Parseval agreement with decoded slot norms for
the saved dense key, the saved sparse key and a dense key with one coefficient
changed. The retained passing trial still meets the 23-bit gate (maximum error
3.91e-9); the incorrect-key controls produce RMS errors 4.24e19 and 4.27e19.
For uniform Q108 phase coefficients at scale48, Parseval predicts slot RMS
near 4.26e19. The two historical failures have RMS ratios 1.00067 and 0.99849
to that prediction. This scale is consistent with lost phase correlation;
it does not prove an incorrect-key cause. The original failing ciphertexts
and keys were not retained. The mode reads only the checkpoint prefix and
output, rejects small-ring profiles, and returns 1 if its correct-key accuracy
or incorrect-key control checks fail. Evidence is in
`narrow-dd-phase-spectrum-validation.json`.

`ckks_narrow_checkpoint_inspect --selector-phase CHECKPOINT.bin GROUP.bin`
checks the full N16384 plain37 sigmoid checkpoint's fresh expanded selector
groups. It reads only metadata and the dense/sparse secrets from the checkpoint,
then streams eight public selectors at a time. Each selector is decrypted with
independent exact integer NTT multiplication and compared with the head/tail
indicator mask encoded at scale 2^48. The fixed 39-unit coefficient-noise gate
is a diagnostic 12-sigma envelope for freshly encrypted Gaussian selectors;
this mode does not validate reconstructed packed selectors. Unsupported
metadata, lengths and sparse geometry return 2; indicator mismatch returns 1.
The evaluator and evaluation-key archive suffix remain outside this diagnostic.
All 36 groups from the second retained full sigmoid trial pass: 19044 selectors
and 312016896 phase coefficients, including the final block's 15 padding
candidates and all unused tail regions. Maximum coefficient noise is 18.
The sequential inspections take 214.59 seconds in total and peak at 37728 KiB
RSS. Malformed metadata is rejected, an all-zero group fails the accuracy gate,
and the existing refresh-prefix mode still reproduces its archived metrics
exactly. Evidence is in `narrow-dd-retained-selector-phase-validation.json`.
The full second-trial refresh passed independently; the earlier unexplained
sigmoid failures remain unresolved.

`build/test/ckks/ckks_narrow_factor_checkpoint CHECKPOINT.bin SELECTOR_DIRECTORY`
checks the full plain37 factor stage using the retained sigmoid trial's final
head and tail groups. It generates only the factor conjugation key under the
saved dense secret, without reading the checkpoint evaluation-key suffix.
Twelve columns use frequencies 1 through 12 with different normalized-mask
tables. Chunk sizes 40 and 17 give identical ciphertexts, and columns 0/5/11
match single-column evaluations exactly. All 196608 decoded slots match the
expected signed table rotations and indicator windows within the fixed 23-bit
gate; maximum error is 1.26e-9. The run takes 150.39 seconds and peaks at
8095236 KiB RSS under a 12 GiB address-space limit, including unused normal key
container storage. It uses production table generation and slot decryption,
and a fresh conjugation key. It does not validate the saved derived cache or
resolve the earlier end-to-end failures. Evidence is in
`narrow-dd-full-ring-factor-checkpoint-validation.json`. No arguments skip
this opt-in diagnostic.

Current LL14 plain37 key-exposure accounting confirms dense selectors at
385 bits, factor automorphism keys at 376 bits, DD tree keys at
348/300/252/204/156 bits, and dense S2C keys at 113/76 bits. Packed provisioning
raises the largest dense exposure to 436 bits. The dense-to-structured switch
uses 76 bits. Native 448-bit storage is separate from these active moduli.
The pinned dense 436-bit estimator screen passes at 128.67 bits. The original
four 76-bit sparse proxy attacks returned nonfinite results at the default
search limit. These screens do not
model public block constraints, random sign counts, ring attacks or correlated
evaluation keys; security remains unqualified. Details and source/report
provenance are in `narrow-dd-current-cpu-key-exposure-audit.json`.

The estimator wrapper now retains raw nonfinite optimizer results. Revalidation
of the LL14 Q76 sparse proxy returns `+Infinity` from all four attacks, rather
than estimator exceptions; the pinned configuration's global maximum block
size is 1754. These outputs remain incomplete screens and do not establish a
security bound. The dense Q436 positive control reproduces all four earlier
finite costs exactly, including its 128.67-bit weakest cost. Evidence is in
`narrow-dd-sparse-estimator-raw-result-validation.json`; public block structure,
ring attacks and correlated evaluation-key assumptions remain unqualified.

The wrapper accepts `--max-beta` for a search-limit sensitivity check and
`--attack` for explicit attack coverage. It updates and verifies the runtime
limits in both primal and dual modules: changing `conf.max_beta` alone would
leave their imported aliases unchanged. At a limit of 4096, the original four
LL14 Q76 sparse-proxy attacks return finite costs (857.23--899.18 bits).
Adding the plain BDD hybrid attack lowers this proxy to 153.31 bits; direct
MITM returns 189.88 bits, and BDD MITM hybrid returns 134.67 bits. The
default dense control retains its four earlier costs exactly. A selected
subset is labeled separately from completion of the original four-attack
screen, and every report retains `security_validated: false`. These are
uniform-sparse, fixed-sign-count proxy estimates; the public block constraint,
random signs, ring structure and correlated-key assumptions are still outside
the model. Evidence is in `narrow-dd-sparse-search-limit-validation.json`.

`Parameter-Selection/python/estimates/narrow_ckks_block_hybrid.py` evaluates one
effect of the public block constraint. It computes exact cumulative guessing
support sizes and hit probabilities from the implementation's independent
uniform positions and signs. Exhaustive toy enumeration checks 32 partitions,
including fully guessed blocks and unequal widths. At the uniform proxy's
fixed lattice parameters, balanced coordinate guessing gives 152.93 bits for
plain BDD hybrid and 133.79 bits for BDD MITM hybrid; prefix guessing is more
expensive at those points. The uniform baselines reproduce the earlier costs
exactly. These are fixed-parameter sensitivity results, not full attack
optimization: GSA uses average unconditional marginal variance, and conditional
and per-block variances, Babai/MITM heuristics, ring structure and correlated
keys remain unqualified. Exact counts alone do not establish security.
Evidence is in `narrow-dd-block-hybrid-probe-validation.json`.

The same script computes exact residual second moments conditional on the
accepted guessed-weight class. Exhaustive toy checks compare the residual
norm and each remaining coefficient's second moment against enumerated secrets.
`--geometry conditional-mean --baseline-report SAVED_PROBE` feeds the resulting
average residual variance into a fixed-parameter GSA sensitivity check. At the
saved LL14 points, balanced guessing gives 163.79 bits for plain BDD hybrid
and 134.48 bits for BDD MITM hybrid. The default unconditional mode reproduces
all earlier fixed-point costs exactly. This scalar variance adjustment does
not model heterogeneous variances, conditional tails or correlations, and
does not optimize lattice parameters. It cannot qualify actual security.
Evidence is in `narrow-dd-block-conditional-geometry-validation.json`.

`Parameter-Selection/python/estimates/narrow_ckks_block_hybrid_grid.py` searches
specified offsets in BKZ block size, guessing dimension and lattice dimension
with exact guessed-weight strata. For guessed weight w, the residual squared
norm is exactly h-w; candidate counts, hit probabilities and per-block residual
second moments pass 32 exhaustive toy partitions. Three LL14 BDD MITM hybrid
grids evaluated 300 points across balanced and prefix layouts. The lowest
searched heuristic cost is 134.29 bits at beta182, zeta14688, d3415 and weight10.
Independent grouped-binomial sums verify each grid's best candidate count and
hit probability, and shared points reproduce the anchor control exactly.
The expanded grid covers beyond the earlier boundary minima and reproduces
12 shared points exactly. These costs do not establish a global optimum or a
lower security bound: the lattice model still uses scalar
variance and omits ring and correlated-key attacks. Evidence is in
`narrow-dd-block-hybrid-grid-validation.json`.

The grid script's `--saturated-blocks` option fully guesses specified numbers
of blocks and balances the remaining coordinates, with prefix and
shortest-block-first orders. At beta182/zeta14688/d3415, 56 layouts include
the original controls and 1..27 saturated blocks. The controls match exactly;
an independent grouped-binomial oracle verifies the best count and probability.
The lowest unrotated heuristic cost is 134.24 bits for five saturated prefix
blocks. This remains a fixed-point scalar-GSA estimate. Evidence is in
`narrow-dd-block-layout-validation.json`.

Ring rotations are a separate missing factor in the earlier LWE probes.
The [2026 ring-hybrid research](https://eprint.iacr.org/2026/279) and its
[author implementation](https://github.com/TabOg/mlwe-hybrids) provide a
rotation-amplification model. The adapter
`Parameter-Selection/python/estimates/narrow_ckks_ring_hybrid_probe.py` pins
author revision a92f335 and its estimator submodule 6019056. At the saved
uniform-support LL14 point, BDGL16 reproduces the old LWE cost 134.67 bits;
the ring estimator-MITM model gives 131.37 bits, and the ring square-root-MITM
model gives 125.88 bits. MATZOV gives 125.87 bits for the latter. Exact uniform
counts pass 27 exhaustive cases, all 192 evaluated costs are finite, and
512/1024-bit probability calculations give identical best costs. These
uniform-proxy results fall below 128, but do not prove an attack on the port:
its public one-per-block support distribution is not rotation invariant.
Joint rotation-event probabilities require a separate public-block model;
the checks below address one guessed-weight event, while compatible guessing
dictionaries remain unconstructed. The earlier unstructured probes
cannot support a 128-bit claim. Evidence is in
`narrow-dd-ring-hybrid-probe-validation.json`.

`Parameter-Selection/python/estimates/narrow_ckks_block_rotation_probability.py`
computes exact per-rotation probabilities and labelled candidate counts for
the actual block distribution. Its independent posterior draws estimate the
probability that at least one rotation has the specified guessed weight,
using `Pr[K>0] = sum_r Pr[event_r] E[1/K | event_r]`. This avoids sampling the
rare unconditional event. Exact marginals, signed projections, posterior path
probabilities and the union identity pass 256 exhaustive toy cases. FFT-based
integer rotation counts are checked against direct counts during sampling.
At the five-saturated-block layout and guessed weight11, two 4000-sample seeds
give union probabilities 1.95e-10 and 1.92e-10; treating rotations as independent
overestimates these by about six times. Their difference is 0.52 combined
standard errors. Approximate Monte Carlo and iid Hoeffding intervals are
retained in `narrow-dd-block-rotation-validation.json`.

`Parameter-Selection/python/estimates/narrow_ckks_block_ring_cost.py` uses those
probabilities with the exact rotation-labelled candidate count and an assumed
full square-root guessing speedup. It reproduces the unrotated control exactly.
Three bounded beta/d grids evaluate 394 finite costs across BDGL16/MATZOV,
including probability endpoints and controls. The lowest searched central cost
is 128.96 bits; the attacker-favorable iid probability endpoint gives 128.78.
This is one zeta/layout/weight, not a global optimum or security qualification.
The probability check measures secret-weight events only: joint nearest-plane,
error and MITM success across simultaneous matching rotations remains
unmodeled. Compatible isometric guessing dictionaries, heterogeneous GSA
geometry and correlated evaluation keys are also unqualified. Evidence is in
`narrow-dd-block-ring-cost-validation.json`.

`Parameter-Selection/python/estimates/narrow_ckks_block_ring_joint_decode.py`
extends that sensitivity model using the posterior histogram of matching
rotations. For decoder marginal `p` and `K` matches, it compares a common
decoder outcome, conditional independence with success `1-(1-p)^K`, and the
optimistic union bound. The independent case uses the exact posterior identity
`sum_r Pr[event_r] E[(1-(1-p)^K)/K | event_r]`; it does not treat the secret
rotation events as independent. Exact Bernoulli enumeration and monotone
summation-by-parts checks pass. A one-sample
[DKW-Massart empirical-CDF bound](https://arxiv.org/abs/1107.5356) covers all
decoder marginals searched on the same fixed-layout histogram simultaneously:
`epsilon = sqrt(log(40)/(2*samples))`, with expectation error at most
`p*epsilon` for the independent decoder function and `epsilon` for `1/K`.
It does not cover additional layouts, weights or separately sampled seeds.
These remain hypotheses about decoding correlations and guessing work, rather
than established attacks or a security confidence interval. Compatible
rotation-labelled MITM dictionaries and correlated evaluation keys remain
unqualified. Two initial grids evaluate 960 finite costs across the two cost
models and three correlation scenarios. The updated fixed-layout grid adds
840 finite costs, including simultaneous 95% posterior uncertainty endpoints.
MATZOV gives minima of 128.96 bits for the common decoder case, 128.56 for
conditional independence, and 128.08 for the optimistic bound; their
attacker-favorable uncertainty endpoints are 128.78, 128.44 and 127.97 bits.
The optimistic endpoint below 128 prevents qualification; it does not establish
an attack. The pinned rerun reproduces every cost from the earlier confidence
grid, and a separate estimator copy gives the same costs. Evidence is in
`narrow-dd-block-joint-decode-validation.json`.

The joint-decoder adapter also accepts the posterior sampler's exact weight
instead of fixing it to the saved weight11 optimum. Adjacent weight9/10/12
searches, a second weight10 seed, a weight9 boundary expansion and a default
weight11 regression produce 4212 finite costs. Weight10 gives conditional
independence estimates of 127.70 and 127.67 bits for the two seeds; common
decoder outcomes give 128.42 and 128.32 bits. The seed estimates differ by
0.68 combined standard errors. Weight9's optimistic central bound reaches
126.88 bits; this remains an unconstructed attack sensitivity. An independent
bit-mask/grouped-binomial oracle verifies exact event and signed labelled
candidate sums at all 16384 rotations for weights9/10/11/12 and passes 576
exhaustive toy cases. Independent Decimal calculations check 2106 joint
probabilities; all 12 default-stratum control costs are unchanged. Each sampled
histogram has its own simultaneous CDF coverage; no combined 95% confidence
claim is made across weights or seeds. Other masks/zeta, actual decoding
correlations and compatible MITM dictionaries remain unqualified. Evidence is
in `narrow-dd-block-weight-joint-decode-validation.json` and
`narrow-dd-block-weight-independent-count-validation.json`.

The optimized CPU executable has now passed a fresh full-ring LL14 plain37
sigmoid trial without a pre-evaluation checkpoint. Maximum error is 4.24e-9
(27.81 bits), RMS error 1.41e-9, and the retained-square maximum error is
9.69e-9. Refresh took 3885.64 seconds; total wall time was 1:24:43 with
14387264 KiB peak RSS. The independent Python importer checks the unchanged
23-bit gates and square bounds, together with archived GeLU/tanh reports.
Source and executable hashes match the launch record; all report numbers are
finite, stderr is empty, and the selector directory was removed. No failure
snapshot was attempted. This fresh pass does not explain the historical
catastrophic failures. Evidence is in
`narrow-dd-current-cpu-capture-completion-validation.json`.

An opt-in LL14 plain37 candidate reduces the auxiliary switching lift from
28 to 16 bits. Its selector, DD tree, input/output levels, Gaussian noise,
domain and fixed 23-bit gates are unchanged. Sparse-key exposure falls from
76 to 64 bits; factor automorphism keys use Q364 and S2C input keys Q101.
The commands are `--small-ll14-plain37-aux16-fbs-file` and
`--ll14-plain37-aux16-fbs-file`, with the same function/cache/report arguments
as the original file mode. The matching `-checkpoint`, `-replay` and
`-replay-trace` modes use a distinct checkpoint magic so that different gadget
row counts cannot be deserialized as the original profile. Prefix inspection
accepts both known formats without reading evaluation-key suffixes. Automatic
failure capture also uses the matching format. The Python candidate importer
requires `--auxiliary-switch-bits 16` to validate this profile; its default
continues to validate the original aux28 profile.
The average-case noise model predicts nearly identical root RMS, but excludes
FFT error and full functional-series behavior. The bounded sparse64 decoder
screen checks 3696 finite costs for weights9/10/11/12; its lowest optimistic
uncertainty endpoint is 130.93 bits. The default sparse76 control costs are
unchanged. These are parameter studies, not actual structured-secret security
qualification; the candidate requires full-ring accuracy validation before
adoption. Evidence is in `narrow-dd-aux16-security-candidate-validation.json`.
`Parameter-Selection/python/estimates/narrow_ckks_auxiliary_lift_scan.py`
reproduces the plain37 factor/tree screen for auxiliary lifts 0 through 32
in steps of four. Auxiliary8 has model RMS 1.0971e-8, about 1.012 times
auxiliary16; nonlinear coefficient weights, maximum-slot tails and FFT error
are excluded. Auxiliary8 is an experimental parameter choice whose full-ring
encrypted accuracy remains unverified.
A separate q56 scalar-GSA study expands the block-size and dimension search
for the existing zeta14560/saturated5/weight12 posterior, then refines its
interior minimum at beta310/d3671. The lowest optimistic estimate is 153.35
bits. This is a single-mask sensitivity, with unproved decoding correlations,
guessing work and correlated-key security. Independent Decimal checks reproduce
2352 saved probabilities and all minima; 52 unavailable decoder points in the
wider grids remain explicitly unvalidated. Evidence is in
`narrow-dd-aux8-joint-grid-validation.json`. Neither study qualifies security
or changes the running auxiliary16 candidate.
The separate `ckks_narrow_auxiliary_probe NEW_DIRECTORY` diagnostic uses N512,
the same plain37 levels, Gaussian sigma3.2, component domain [-4,4], fixed
23-bit gates and one retained DD square. It passes sigmoid/GeLU/tanh with
maximum errors 1.28e-9, 2.85e-9 and 1.55e-8. The three functions share one
fresh key set; this is a small-ring diagnostic, not full-ring qualification.
`--inject-failure` exercises retention under a distinct auxiliary8 probe
archive magic. The state includes secrets, input/output and keys in a private
directory; it is incompatible with CP01/CP16 benchmark checkpoint readers.
An independent Python prefix inspector uses integer negacyclic convolution
and direct canonical polynomial evaluation to reproduce its fresh/refresh
metrics. It checks the ternary and one-per-block secrets, but excludes the
evaluation-key suffix and retained square. Both Python validators reject six
corruption controls. Evidence is in `narrow-dd-aux8-small-validation.json`
and `narrow-dd-aux8-capture-prefix-validation.json`.
The probe now accepts `NEW_DIRECTORY --full` to run all three functions at
N16384 with the same fixed gates and shared fresh key set. Its reports expose
source/selector/S2C/DD/conjugation levels and retain failed states under the
probe format. GNU compilation and linking disable IPA-ICF, matching the
benchmark's guard for parameter types sharing coefficient storage. The extended
executable passes another small all-function run and injected-capture prefix
check (`narrow-dd-aux8-small-v2-validation.json` and
`narrow-dd-aux8-capture-v2-prefix-validation.json`). The full run used four
threads and a 16 GiB address-space limit, under
`narrow-dd-aux8-full-memory-budget.json`, and exited 1 after sigmoid failed
its accuracy gates. It retained the complete diagnostic state and all 36
selector files; GeLU and tanh were not evaluated. This empirical limit is above the
measured auxiliary16 analog VM peak, with fewer gadget rows and worker stacks;
it is not an allocator worst-case proof. The seven generic/uniform-sparse
q56 proxy attacks complete with weakest estimate143.99 bits, recorded in
`narrow-dd-aux8-uniform-security-validation.json`; public block/ring structure
and correlated evaluation keys remain excluded from that proxy claim.
The independent auxiliary8 prefix inspector now accepts both N512 and N16384
probe archives. At full size it uses 15-bit coefficient limbs and a separate
31-bit scratch-prime NTT for exact integer negacyclic decryption, followed by
a Python canonical FFT. The signed-convolution bound fits half the scratch
prime even for a fully nonzero ternary key. Fifteen toy cases match scalar
integer references; full algebraic fixtures check32768 phase coefficients and
16384 slots against closed-form convolution and analytic polynomial evaluations.
A synthetic full prefix also checks the reader. These fixtures do not execute
FBS or inspect evaluation-key suffixes. The real small captured prefix still
matches; nine additional malformed report/state controls are rejected. An
injected large phase error remains a failed23bit case while its metrics compare
using absolute/relative1e-11 tolerance. Original captured secrets/ciphertexts
remain unchanged. Evidence is `narrow-dd-aux8-prefix-full-arithmetic-validation.json`
and `narrow-dd-aux8-prefix-small-regression.json`.
The real full auxiliary8 sigmoid failure is now retained. Independent exact
phase arithmetic confirms maximum error 1.73317e20 and RMS error 4.23590e19;
fresh-input error remains about 6.1e-11. The run completed normally with exit 1
in 1:55:58, with peak RSS 13740352 KiB. The 6449139279-byte state and original
report are hashed, and launch source/executable hashes match. The independent
prefix inspection does not validate the evaluation-key suffix or retained
square, and same-state replay and failure diagnosis remain pending. Evidence
is `narrow-dd-aux8-full-sigmoid-failure-validation.json` and
`narrow-dd-aux8-full-failure-prefix-validation.json`. This failure does not yet
explain the earlier unsaved sigmoid failures.
The exact diagnostic reader accepts this probe archive with
`--aux8-checkpoint`, `--aux8-checkpoint-tree`,
`--aux8-checkpoint-output-conjugation`, and
`--aux8-checkpoint-factor-galois`. It validates the appended retained-square
archive boundary before seeking the refresh-key suffix. The actual failing
state's five DD keys pass 5160960 key-phase coefficients, 163840 exact
arithmetic coefficients and 81920 switched-phase coefficients. All 15 Q356
factor keys pass 21872640 phase coefficients; the Q116 output-conjugation key
passes 475136 phase coefficients and three encrypted checks below 1.4e-10.
Three malformed-profile/suffix controls are rejected, and the older CP01
output-conjugation regression passes. These checks do not establish the failure
cause or validate S2C/encapsulation keys, selector phases, the retained-square
key phase, or same-state functional replay. Evidence is
`narrow-dd-aux8-retained-oracle-validation.json` and
`narrow-dd-aux8-retained-output-conjugation-validation.json`.
`build/test/ckks/ckks_narrow_auxiliary_replay --prepare STATE SELECTORS NEW_DIRECTORY`
loads the probe archive and measures dense S2C, sparse encapsulation, and the
final scale-down separately. `--replay` reevaluates the saved input with the
saved keys and selectors; `--trace` additionally records positive-frequency
root errors. Both replay modes save a private ciphertext prefix for independent
phase inspection and report whether the saved output is reproduced byte for
byte. Exit zero means the diagnostic completed; accuracy is reported separately.
The N512 retained capture passes all three modes. Plain and traced replay match
the original ciphertext exactly, independent phase arithmetic verifies both,
and the original state and all 35 selector files remain unchanged. Three
malformed inputs are rejected before creating output directories. Evidence is
`narrow-dd-aux8-small-replay-validation.json`. Full retained-failure preparation
completed in 12:12.86 with 8427584 KiB peak RSS. Actual S2C, sparse switching,
and scale-down errors are below 2.80e-9. Independent Python integer phase
arithmetic verifies all 16384 prepared coefficients and agrees exactly with
the reported maximum error. Four malformed ciphertext/report/profile controls
are rejected. Evidence is `narrow-dd-aux8-full-prepare-validation.json` and
`narrow-dd-aux8-full-prepared-phase-validation.json`.
The exact selector phase check passes groups 0 through 18, but fails one
ciphertext in group 19. The native `--aux8-selector-phase` reader confirms
10816 erroneous coefficients with errors reaching 152 bits. Separate Python
integer NTT arithmetic finds a high error equal to a negacyclic rotation of
the dense secret times 2^151. Correcting one mask bit in memory leaves a zero
message phase with noise at most 14; original selector files remain unchanged.
This describes the error shape, without establishing where the bit error
originated or its effect on the saved functional output. All 36 original groups
have now been checked: 19044 ciphertexts and 312016896 phase coefficients,
with exactly one failing ciphertext. A private diagnostic selector directory
changes only that one mask bit; all 529 phases in its corrected group pass.
The original group hash is unchanged, and other groups are linked without
modification. A full traced replay using unchanged selectors has completed;
a controlled functional replay with the diagnostic correction now passes against
the reproducible current baseline, as described below.
Evidence is `narrow-dd-aux8-native-selector-failure-validation.json` and
`narrow-dd-aux8-retained-selector-corruption-validation.json`, with complete
group coverage in `narrow-dd-aux8-all-selector-phase-validation.json` and the
private correction check in
`narrow-dd-aux8-selector-correction-preparation-validation.json`. Same-state
archived-output reproduction and the selector error's origin remain unresolved.
The Python diagnostic runner
`Parameter-Selection/python/estimates/narrow_ckks_corrected_selector_replay.py`
queues the controlled one-bit comparison after the original trace. It waits
for the original process identity to terminate, requires a successful report
with byte-identical saved output, and rechecks source/executable, state and
corrected-group hashes before launching under a fresh 16 GiB budget. The saved
keys and input are reused. Its output prefix is independently checked after
completion; exit zero alone does not establish accuracy. Three baseline controls
(failed run, nonmatching ciphertext and wrong profile) reject before loading
keys. The queue rejected the completed trace because its ciphertext differed
from the archived output. Results are in
`narrow-dd-aux8-corrected-replay-validation.json`. The untraced original-selector
replay completed and has the same
5,636,193-byte prefix as the trace, while both differ from the archived output.
An explicit `--matching-baseline-validation` option accepts this reproducible
current pair after checking distinct prefix files, byte identity, profiles and
independent metrics. Its five rejection controls pass. The default gate still
requires archived-output identity. The corrected untraced replay completed against
that current baseline with unchanged saved keys/input; results are in
`narrow-dd-aux8-current-corrected-replay-validation.json`. Current reproducibility
does not explain the archived mismatch or establish the selector error's origin.
The corrected N16384 sigmoid passes the fixed 23-bit gate: independently decoded
maximum refresh error is 4.1792640e-9 and RMS error is 1.4105323e-9, compared with
the original current baseline's maximum error of 1.3212046e20. The single-bit
intervention's observed effect is established for this retained fixture; the
archived output is still not reproduced and the bit error's origin is unknown.
The original retained Q108 multiplication key also produces a passing Q60,
scale-2^48 square: independent maximum error 9.2169691e-9, RMS error 2.4462852e-9,
under the propagated product bound 4.5031960e-7. Evidence is
`narrow-dd-aux8-current-corrected-square-validation.json`. This validates this
diagnostic sigmoid fixture, not the full aux8 profile across all functions or
its structured-key security.

`build/test/ckks/ckks_narrow_auxiliary_functions FUNCTION STATE SELECTORS NEW_DIRECTORY`
applies sigmoid, GeLU or tanh to the original retained aux8 sigmoid input using
the saved refresh and multiplication keys. GeLU uses the existing degree-12
trigonometric approximation and tanh the degree-40 windowed approximation;
tile size stays 12. The saved output prefix carries the new function tag.
`../Parameter-Selection/python/estimates/narrow_ckks_auxiliary_functions_validate.py`
accepts `--state`, `--directory`, and `--output`; it compares the original
secrets, slots and input ciphertext bytes, independently decodes refresh and
square phases, and checks the reported approximation/evaluation gates separately.
The square phase reader also accepts `--function sigmoid|gelu|tanh`, with sigmoid
as its default. All three N512 retained functions pass, and the sigmoid prefix
is byte-identical to the earlier small replay. Thirteen reader/profile controls
and eight full-launcher controls reject invalid inputs. Evidence is
`narrow-dd-aux8-retained-functions-build-small-validation.json` and
`narrow-dd-aux8-retained-functions-controls-validation.json`.

The full retained-function launcher
`../Parameter-Selection/python/estimates/narrow_ckks_auxiliary_functions_run.py`
requires `--repo`, `--state`, `--destination`, `--build-validation`,
`--retained-validation`, `--correction`, `--budget`, and `--output`. It runs GeLU
then tanh by default; `--functions` selects an explicit sequence. Live memory
reservations use `--reserved-job PID:START_TICKS:CAP_GIB`. It checks the passing
small diagnostic and full corrected sigmoid evidence, rechecks state/program
and corrected-group hashes, and obtains a fresh 16 GiB child budget before each
function. Full GeLU/tanh validation is running; results are in
`narrow-dd-aux8-retained-functions-full-validation.json`. This extends function
coverage for the saved fixture without regenerating keys or identifying the
original selector bit error's origin.
Gate checks are in
`narrow-dd-aux8-corrected-replay-gate-validation.json`.

The completed full trace has maximum refresh error 1.32120e20 and RMS error
4.26451e19, independently confirmed by integer prefix-phase arithmetic. All
24 observed weighted roots exceed the 23-bit threshold. These observer metrics
do not independently validate root phases or establish the selector error's
causal role. Evidence is `narrow-dd-aux8-full-trace-prefix-validation.json` and
`narrow-dd-aux8-full-root-trace-validation.json`.

`build/test/ckks/ckks_narrow_auxiliary_square STATE REPLAY_PREFIX NEW_DIRECTORY`
uses the original aux8 sigmoid probe's retained multiplication key to square
the replayed output. It requires identical secrets, target slots, input and
profile in both prefixes, exact replay-prefix termination, and a checked
compact-key suffix. It writes a private `square.bin` at Q60 with scale 2^48
and reports the existing propagated 23-bit product-error gate. Diagnostic
exit zero is separate from accuracy. Run the independent phase checker from
the TFHEpp directory:

```sh
python3 ../Parameter-Selection/python/estimates/narrow_ckks_auxiliary_square_validate.py \
  --state STATE --square NEW_DIRECTORY/square.bin \
  --report NEW_DIRECTORY/report.json --output VALIDATION.json
```

The retained N512 fixture passes with maximum square error 3.16608e-9;
independent integer phases agree. Three C++ archive/input controls and four
Python reader controls reject malformed inputs. Applying the same retained
key to the failed N16384 traced output gives square error 437621.8432 and
fails the product gate, independently confirmed as a negative control. This
does not validate a corrected full output or every multiplication-key row.
Evidence is `narrow-dd-aux8-retained-square-build-small-validation.json` and
`narrow-dd-aux8-full-trace-retained-square-validation.json`.

The corrected-selector square launcher
`../Parameter-Selection/python/estimates/narrow_ckks_corrected_selector_square.py`
accepts `--repo`, `--state`, `--destination`, `--replay-validation`,
`--square-provenance`, `--budget`, and `--output`. Use the completed corrected
replay evidence for `--replay-validation` and the retained-square build evidence
for `--square-provenance`. Each live memory reservation is supplied as
`--reserved-job PID:START_TICKS:CAP_GIB`. The destination must be new. The launcher
can wait for a matching replay launcher with `--wait-for-replay PID:START_TICKS`;
it checks the completed independent validation after that process terminates.
The waiting parent has a 128 MiB address-space limit and starts no square work.
The launcher
checks independent replay metrics and provenance before evaluation, uses a
2 GiB child limit, then independently decodes the resulting square ciphertext.
It records refresh and product accuracy separately from process exit status;
a failed refresh remains usable as a negative control. Its terminal, profile,
and metric rejection controls are recorded in
`narrow-dd-aux8-corrected-square-launch-gate-validation.json`. These controls
do not establish the corrected full-ring square accuracy.

`build/test/ckks/ckks_narrow_dd_exact --aux8-checkpoint-square STATE`
checks the actual retained N16384 Q108 multiplication key. Independent integer
phase arithmetic verifies all 442368 encrypted key-row coefficients against
the gadget-scaled squared secret; maximum observed noise is 15 under the
39-unit diagnostic envelope. The FFT switch matches all 32768 exact NTT
coefficients, and all 16384 switched-phase coefficients satisfy the observed
key-noise convolution bound. A one-bit change to a top secondary body digit
passes archive checks but fails the phase check. The reader skips the middle
refresh-key archive. These checks do not establish the full refresh failure's
cause or qualify security. Evidence is
`narrow-dd-aux8-retained-square-key-oracle-validation.json`.
`Parameter-Selection/python/estimates/narrow_ckks_auxiliary_root_trace.py`
checks the 24 weighted roots emitted by the degree-12 complex sigmoid replay.
It enforces frequency order, finite metrics, identical real/imaginary weights,
and the amplitude bound implied by each reported error. With `--partial`, it
reads only complete lines and records the first root exceeding the fixed
23-bit threshold; an empty or incomplete trace is never reported as a full
pass. The retained N512 trace passes all 24 observations. Six malformed
controls, a truncated live line, and a deliberately inaccurate root exercise
the checker. This validates observer record consistency, not independent root
decryption or final refresh accuracy. Evidence is
`narrow-dd-aux8-small-root-trace-validation.json` and
`narrow-dd-aux8-root-trace-controls-validation.json`.
The full retained-failure baseline has emitted its first completed root:
frequency 1, weight magnitude 0.26998, maximum error 1.26074e20. The trace
checker accepts the record and flags the fixed precision failure. This observes
an error before final output combination, while independent prepared-coefficient
checks pass below 2.80e-9. Root decryption is still the C++ observer path;
independent root phase validation, completed baseline byte equality and the
one-bit corrected functional comparison remain pending. Evidence is
`narrow-dd-aux8-first-root-failure-validation.json`.
The small N512 candidate passes all three fixed accuracy/square gates:
sigmoid 29.63 bits, GeLU 28.40 bits, and tanh 25.95 bits. Injected failure
capture retains an independently inspectable checkpoint, and replay produces
the saved ciphertext byte for byte. Both mixed-profile checkpoint reads are
rejected; an older aux28 checkpoint still replays byte for byte. Evidence is
in `narrow-dd-aux16-small-validation.json`. Sequential full-ring candidate
validation is running under `narrow-dd-aux16-full-memory-budget.json`. Full
GeLU and sigmoid pass their independent fixed 23-bit refresh/square gates;
full tanh remains pending.

The sparse64 study now also covers weights13/14/15 and four neighboring public
guessing masks: zeta14560/14816 with five saturated blocks, and zeta14688 with
four/six saturated blocks. The descriptor generator
`Parameter-Selection/python/estimates/narrow_ckks_block_guessing_geometry.py`
preserves actual secret block widths and supplies only an anchor/mask, with no
optimality claim. Independent grouped-binomial verification checks all rotation
mass and labelled support sums for each mask. Across the original and extended
sparse64 grids, 12924 cost evaluations are finite and 12 are nonfinite. The
nonfinite outputs occur at beta174/d3671 for zeta14560 with a decoder marginal
about 7.57e-324; their raw outputs are retained as missing estimates, never
security bounds. The lowest searched optimistic posterior endpoint is now
130.86 bits. Separately, all seven conventional uniform-sparse proxy attacks
complete with the pinned estimator and block-size limit4096; the weakest is
BDD MITM hybrid at 139.94 bits. These proxy/global optimizer outputs do not
replace the actual-block ring studies or establish security. Per-histogram
coverage is not a joint confidence claim across the searched masks and weights.
Evidence is in `narrow-dd-aux16-extra-weight-validation.json`,
`narrow-dd-aux16-mask-security-validation.json`, and
`narrow-dd-aux16-uniform-security-validation.json`.

The 12 missing sparse64 costs above were reproduced as a Sage coercion issue:
multiplying the positive MPFR joint probability by a Real Double Field decoder
factor rounds it to zero. The adapter preserves the raw backend result and,
only when that multiplication is zero, reconstructs the same pinned BKZ/SVP
formula with a 512-bit joint probability and stable `log1p` amplification.
All 12 cases are now finite (1199–1202 bits), giving 12936 finite original or
repaired evaluations; the lowest searched endpoint stays at 130.86 bits.
Six finite backend controls and independent 2048-bit references pass. The
pinned estimator is unchanged. This numerical repair does not qualify the
structured-secret attack model or security. Evidence is in
`narrow-dd-joint-probability-numerical-validation.json`; reproducible validation
is `Parameter-Selection/python/estimates/narrow_ckks_joint_probability_validate.py`.

The exact retained-key inspector supports the opt-in 16-bit auxiliary profile
with `--aux16-checkpoint`, `--aux16-checkpoint-tree`,
`--aux16-checkpoint-output-conjugation`, and `--aux16-checkpoint-factor-galois`.
These require the distinct CP16 magic and use Q124/Q364 key payload sizes;
existing modes continue to require CP01 and Q136/Q376. Positive reading of a
retained full CP16 checkpoint remains pending. The standalone
`build/test/ckks/ckks_narrow_dd_exact --aux16-output-conjugation` validates all
507904 full-ring Q124 key-phase coefficients with an independent integer NTT
oracle and exercises random complex, unit-circle and sparse-slot ciphertexts
at the fixed 23-bit gate. Maximum error is below 8.1e-11. The saved CP01 Q136
output-key regression passes and the CP16 reader rejects its incompatible
magic. Evidence is `narrow-dd-aux16-output-conjugation-validation.json`.

`build/test/ckks/ckks_narrow_dd_exact --aux16-factor-galois` checks all
15 fresh full-ring Q364 factor automorphism keys with the independent integer
NTT key-phase oracle, then evaluates each key on random complex, unit-circle
and sparse-slot ciphertexts at Q348/scale48 using the 16-bit auxiliary lift.
The expected result comes from direct substitution of encoded plaintext
coefficients. All 22364160 key-phase coefficients pass (maximum noise18
against the diagnostic envelope39), and all 368640 encrypted slots pass the
fixed 23-bit gate (maximum error1.008e-10). The saved CP01 Q376 factor keys
also pass and reproduce their earlier phase diagnostics exactly. These fresh
candidate keys are independent of the still-running full refresh; its accuracy
and live-key diagnosis remain pending. The evidence validator in
`Parameter-Selection/python/estimates/narrow_ckks_factor_oracle_validate.py`
rejects missing/duplicate cases, nonfinite errors and signal termination.
Evidence is `narrow-dd-aux16-factor-oracle-validation.json`.

The projection-dictionary construction in
`Parameter-Selection/python/estimates/narrow_ckks_block_mitm_dictionaries.py`
now provides exact signed partial dictionaries at rotation zero for weights9–15
and the four neighboring masks. It assigns forced blocks to disjoint sides,
splits optional nonzero entries by weight, and filters pairs sharing a block.
Every full projection has exactly binomial(w-s,k) compatible decompositions.
Suffix-DP rank/unrank constructs individual entries without materializing the
full dictionaries. Thirty exhaustive signed toy strata prove coverage and
multiplicity; an independent grouped-binomial oracle matches all 11 actual
full/partial counts, and 44 full-size sampled projections decompose correctly.
For zeta14560/saturated5/weight12, the selected dictionaries enumerate
2^73.59 entries, 4.36 bits above sqrt(single-rotation support). This specific
construction's entry overhead is not a security margin or a lower bound on
other attacks. The existing cost screen is unchanged: lattice-compatible
Babai collision search and isometric reuse across rotations remain unconstructed.
Evidence is `narrow-dd-block-mitm-dictionary-validation.json`.

`Parameter-Selection/python/estimates/narrow_ckks_block_rotated_dictionaries.py`
extends the signed projection dictionaries to every public rotation. Selected
position ordinals map to actual original-block positions, then to X^r*s with
the negacyclic sign at wraparound. All 11 studied weight/mask geometries cover
16384 rotations each. An independent doubled-mask interval oracle and grouped
binomial factors verify all 180224 support counts; direct scalar substitution
checks 88 full-size coordinate samples. Exhaustive signed toy enumeration
covers 170 nonempty strata, including pointwise sign checks so uniform signed
sets cannot hide a sign error. The dictionaries remain implicit rank/unrank
representations. For zeta14560/saturated5/weight12, separately enumerating these
dictionaries requires 2^87.75 entries, while sqrt(labelled support) is2^76.38
and the sum of per-rotation square roots is2^83.35. These counts describe this
construction and do not bound shared algorithms or confer a security margin.
Lattice-compatible Babai collision search and isometric reuse remain unverified;
the existing security screen is unchanged. Evidence is
`narrow-dd-rotated-dictionary-validation.json`.

A two-thread N512 packed-file comparison uses private executable snapshots
of the original live packed job and the current build, with the same fixed
identity profile and accuracy/square gates. Both complete normally and clean
up owned selectors. Reconstruction takes160.22 seconds for the old build and
132.99 seconds for the current build (observed1.20x ratio). One sample per
build with independent keys and concurrent full CPU load does not establish
an isolated optimization effect or full-ring speedup. The original full job
continues with its original mapping; positive CPU/I/O deltas confirm activity.
For N16384 and two cached extraction depths, source-derived tree counts are
49152 cached and16380 uncached automorphism calls, with6002573312 bytes of
derived substitution FFT payloads. Additional depths trade payload memory for
fewer raw-key preparation calls; this comparison makes no cache-policy change
or runtime prediction. The independent validator is
`Parameter-Selection/python/estimates/narrow_ckks_packed_cpu_compare.py`.
Evidence is `narrow-dd-packed-cpu-comparison-validation.json`.

The full N16384 plain37 auxiliary-16 GeLU run completed normally in1:22:14
with maximum refresh error2.1881e-8 (25.45 bits), RMS6.6239e-9 and retained
square error1.6173e-7, passing the independent fixed 23-bit gates. Peak RSS is
14134536 KiB. Owned selectors were cleaned up, no failure snapshot was needed,
and the launch source/executable hashes still match. The parameter importer
accepts this function; the later sigmoid and tanh completions also pass.
The sequential runner has completed all three functions; the candidate is not yet
adopted or security-qualified. Evidence is
`narrow-dd-aux16-full-gelu-completion-validation.json`.

The full auxiliary16 sigmoid run completed normally in1:29:11 with
14135112 KiB peak RSS. Maximum refresh error is3.6567e-9 (28.03 bits), RMS
error1.4046e-9, and retained-square error9.4852e-9. The independent importer
checks GeLU, sigmoid and the later tanh report against the fixed23bit and square gates.
Source/executable hashes match launch;
stderr is empty, owned selectors were cleaned, and no failure snapshot was
needed. The sequential runner has completed full tanh. Historical unsaved
sigmoid failures and actual structured-key security remain unresolved. Evidence
is `narrow-dd-aux16-full-sigmoid-completion-validation.json`.

The full auxiliary16 tanh run completed normally in 4:12:26 with 14135244 KiB
peak RSS. Maximum refresh error is 1.98736e-8, RMS error is 1.09920e-8, and
post-refresh square error is 5.42810e-8. The independent report importer checks
the fixed 23-bit refresh threshold and propagated product thresholds for all
three functions. This checks reported metrics, without independently decrypting
the cleaned-up evaluation keys. The completion evidence is
`narrow-dd-aux16-full-tanh-completion-validation.json`; the combined parameter
report is `Parameter-Selection/python/estimates/narrow_ckks_ll14_plain37_aux16_all_validated.json`.
The profile remains opt-in and its structured-key security is unqualified.

`build/test/ckks/ckks_narrow_packed_cache_profile --create DIRECTORY` saves a
public N512 Q422 packed key, and `--run DIRECTORY DEPTH` compares depths2/3/4
against the ordered depth-two ciphertext terms. Using the same key, extraction
takes111.53/78.63/67.92 seconds respectively; depths3/4 match all2097152
coefficients exactly. Cache payloads are184090624/276135936/368181248 bytes.
An isolated changed term is rejected with exit2. These are single small-ring
trials under concurrent load with different reference write/read I/O; they do
not prove a full-ring speedup. The existing conservative full profile payload
sum grows from47.60 GiB at depth2 to53.19 GiB at depth4. No production cache
policy changed and the original full packed run continues at depth2.
Evidence is `narrow-dd-packed-cache-profile-validation.json`.

`build/test/ckks/ckks_narrow_tensor_exact --full` checks the N16384 tensor
arithmetic at every LL14 DD tree level (348/300/252/204/156), the retained-square
level 108, and unequal 348/300 inputs. Sixteen zero/48-bit rescale cases compare
786432 coefficients (5505024 native words) exactly against an independent
single-prime integer NTT oracle using 16-bit chunks and checked signed capacity.
The oracle uses its own wide integer assembly/rounding and preserves the
balanced 4-bit wrap threshold at ties. Inputs include unused high native bits,
negacyclic wraparound, positive/negative half-way values and values on both
sides of the balanced wrap threshold. The full-ring check passes with no stderr.
This verifies tensor arithmetic; DD key switching is checked separately below,
and the unexplained end-to-end sigmoid failures remain unresolved. Evidence is in
`narrow-dd-full-ring-tensor-exact-validation.json`. No arguments skip this
opt-in workload.

The `--packed` mode of both exact diagnostics covers the full packed profile's
49-bit scale drops and DD relinearization after rescale. Tensor cases use
345/296/247/198/149-bit inputs, including unequal 345/296 inputs and explicit
partial-digit sign/rescale boundaries. All 20 cases match 983040 coefficients
exactly. The independent oracle uses centered active-Q integers when the top
4-bit digit is partial, and the asymmetric balanced wrap threshold when Q is
a multiple of four. The DD check covers generated keys at
296/247/198/149/100 bits: 4079616 key-row phase coefficients, 163840 exact
FFT/NTT output coefficients and 81920 switched-phase coefficients pass.
These checks validate full-ring arithmetic for the packed modulus chain;
they do not validate the live packed selector reconstruction or its complete
refresh. Evidence is in `narrow-dd-packed-profile-exact-validation.json`.

DD FFT preparation and switching enforce the accumulated convolution bound
`Bg + Bbar + nbit + bit_width(primary_rows) + 3 < 53`. The compact-key factory
and the switch accepting a directly supplied FFT key now use the same bound
as the raw-key and chain factories. Checking only an individual row product
could admit a parameter set whose accumulated FFT coefficients exceed the
precision budget. Compile checks accept 4/4 and 18/8 digit profiles and reject
the N16384/Q432 18/16 profile: its individual-product bound is 51 bits, while
its accumulation bound is 56 bits. This compile-time guard does not establish
the cause of the earlier 4/4 sigmoid failures, which had ample headroom.
Evidence is in `narrow-dd-fft-accumulation-headroom-validation.json`.

`build/test/ckks/ckks_narrow_dd_exact --full` compares the production FFT DD
switch against an independent exact integer NTT oracle at N16384 and active
moduli 348/300/252/204/156/108. Each case uses the same generated encrypted
compact key and polynomial in both paths, with 4-bit primary/secondary digits,
dense ternary secrets and Gaussian noise standard deviation 3.2. The oracle
extracts digits and assembles the result independently, streaming secondary key
rows under a checked single-prime convolution bound. Random native-width inputs
include high inactive bits and explicit carry/wrap boundaries. All 196608
coefficients (1376256 native words) match exactly. The initial arithmetic-only
run is recorded in `narrow-dd-full-ring-dd-exact-validation.json`.

The same command now independently verifies the encrypted messages in every
DD key row: exact NTT decryption of 5603328 coefficients recovers the expected
gadget multiple of the squared secret, with maximum noise 14 to 18 across the
six keys. It then checks 98304 switched phase coefficients against
`input*(s*s)` modulo the active modulus, using the observed key-noise convolution
bound. Secret products, phase subtraction and reduction use independent integer
arithmetic. All phase checks and FFT/NTT equality checks pass in 90.87 seconds,
with peak RSS 2462284 KiB and no stderr. The 39-unit key-noise gate is a
diagnostic 12-sigma envelope. These generated-key cases do not resolve the
unsaved failing sigmoid cases or qualify security. Current source hashes and
results are in `narrow-dd-full-ring-dd-phase-validation.json`.
No arguments skip this opt-in workload.

`build/test/ckks/ckks_narrow_dd_exact --checkpoint CHECKPOINT.bin` checks the
actual stored DD keys in a full plain37 checkpoint. Its pre-rescale tree keys
use moduli 348/300/252/204/156. The reader retains the portable archive's
endianness, reads the dense-secret prefix, then seeks to the fixed-size compact
DD chain and final cache-mode field. It validates every compact-key metadata
record and the suffix termination; it does not validate the skipped middle
archive or load the linear/automorphism keys. All five keys from the second
retained sigmoid trial pass: 5160960 key-message coefficients, 81920 switched
phase coefficients and 163840 FFT/NTT coefficients. The run takes 45.25 seconds
with peak RSS 2896636 KiB. Truncated/wrong-magic prefixes are rejected. The
generated-key `--full` mode also passes after the refactor. Evidence is in
`narrow-dd-retained-dd-exact-validation.json` and
`narrow-dd-generated-dd-revalidation.json`. These checks use fresh FFT
preparation; the live derived cache and full product tree remain separate
checks, and the original/O1 sigmoid failures remain unresolved.

`build/test/ckks/ckks_narrow_dd_exact --checkpoint-tree CHECKPOINT.bin`
checks the complete full-ring tree with that stored DD chain. Each case has
31 freshly encrypted unit-root leaves and one transparent complex beta leaf,
followed by five DD multiplication levels with relinearization before rescale.
Head, tail and full-slot cases include beta weights `0.125-0.25i`, `1e-6+2e-6i`
and `1+0.5i`. All 24576 output slots pass the fixed 23-bit gate, with maximum
error 5.45e-10. The 93 multiplication nodes produce Q108/scale48 outputs.
The run takes 116.17 seconds and peaks at 6687348 KiB RSS under an 8 GiB
address-space limit. The DD suffix digest matches the independently checked
stored chain. Leaves and FFT caches are generated for this diagnostic;
it does not capture the live evaluator buffers or establish full FBS stability.
Evidence is in `narrow-dd-retained-dd-tree-validation.json`.

`build/test/ckks/ckks_narrow_dd_exact --checkpoint-output-conjugation CHECKPOINT.bin`
checks the stored Q136 output-conjugation key from a full plain37 checkpoint.
The reader checks the adjacent DD header and final cache mode, then reads the
fixed-size automorphism payload without loading the intervening evaluation
keys. Independent integer NTT phase decryption of all 557056 key coefficients
recovers the gadget multiple of the conjugated secret; maximum noise is 15.
The explicit secret transformation uses `s[0]` and `-s[N-n]`, independently of
production automorphism code. Three fresh complex ciphertext patterns then
use the saved key's FFT cache and production 28-bit auxiliary-scale evaluation
at Q108/scale48. All 24576 slots pass the fixed 23-bit gate; maximum error is
8.44e-11. The run takes 2.17 seconds and peaks at 686016 KiB RSS, including the
unused compact-chain container used to calculate the suffix length.
Malformed prefixes, adjoining DD headers and cache modes are rejected, as is
a synthetic key coefficient exceeding the diagnostic 39-unit noise envelope.
The skipped archive and live evaluator buffers are outside this check;
earlier end-to-end failures remain unresolved. Evidence is in
`narrow-dd-retained-output-conjugation-validation.json`.

`build/test/ckks/ckks_narrow_dd_exact --checkpoint-factor-galois CHECKPOINT.bin`
streams all 15 stored Q376 factor automorphism keys, including the tail-rotation,
identity and factor-conjugation keys. It independently derives `X -> X^d` by
permuting secret coefficients and applying negacyclic signs, then verifies
every encrypted gadget row with exact integer NTT phase subtraction. All
23101440 coefficients from the second retained sigmoid checkpoint pass the
39-unit diagnostic noise envelope. This validates the saved key messages;
FFT cache preparation, rotation composition and live factor buffers are
outside this check. The shared output-conjugation mode still passes after the
reader extension. Evidence is in
`narrow-dd-retained-factor-galois-validation.json`.

## Practical Notes

- `build/test/ckks/ckks_narrow_memory_replay CHECKPOINT_DIRECTORY` is an opt-in
  check for a retained small-ring plain37 sigmoid checkpoint. It reuses the
  saved input and keys with a deterministic degree-12 periodic series whose
  positive-frequency weights are complex, evaluates with frequency tiles of
  12 and 1, and requires identical output ciphertexts, finite errors below
  the fixed 23-bit gate, and no FFT caches retained in the public key. It reads
  the checkpoint and selectors without changing them. This curve differs from
  the saved sigmoid approximation, so the check does not compare the saved
  sigmoid output or establish full-ring correctness. The focused source can
  be built with ASan/UBSan; prebuilt libraries and FFT assembly need separate
  instrumentation, and these sanitizers do not check data races. The normal
  CMake build passed on the retained small fixture: all 256 complex slots had
  maximum error 2.84e-11, tile outputs matched exactly, and all fixture file
  hashes were preserved. The focused ASan/UBSan build also passed with no
  diagnostics, in 17:15 with 536020 KiB peak RSS. Evidence and
  instrumentation limits are in `narrow-dd-memory-replay-validation.json`.
- Plain37 file-selector benchmarks automatically retain a failed accuracy run
  or an exception during functional evaluation/validation, unless an explicit
  checkpoint was already requested. A new private directory named
  `SELECTOR_DIRECTORY.failed-state-N` holds replay keys, secrets, input/target
  slots, `output.ciphertext`, and a failure note. Known public selector files
  are retained with hard links (or copied across filesystems); the caller's
  directory and earlier snapshots are not moved or replaced. `checkpoint.bin`
  is published after all replay artifacts have been written. An exception may
  leave a partial output; `failure.txt` records whether refresh completed.
  Snapshot errors are reported without turning the failed run into a pass.
  This does not capture failures before the functional evaluation phase,
  process signals, evaluator scratch buffers or the fresh post-refresh square
  key. The small-ring `--small-ll14-plain37-fbs-failure-capture-test` command
  explicitly injects an accuracy failure. Its saved ciphertext replays
  byte-for-byte after normal selector-cache cleanup; an earlier snapshot stays
  untouched. A controlled reporting exception also publishes an independently
  inspectable snapshot. These are capture tests, not reproductions of the
  historical sigmoid failure. Evidence is in
  `narrow-dd-failure-capture-validation.json`.
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

Historical aux8 observation, before the original packed run completed: full retained aux8 GeLU completed runner11898 subtask exit0: independent refresh max2.0248871947261327e-8 RMS6.603142031048314e-9 and retained Q60/delta48 square max1.5583841990389361e-7 RMS3.0487194317651916e-8 PASS bound1.4678657411332307e-6. Original retained secrets/slots/input byte identity verified; one corrected selector bit, no fresh keys. Wall1:08:37 peak13708844KiB. Tanh sequentially launched after fresh provenance and memory checks; packed reconstruction remains live. Archived mismatch/bit-error origin/security remain unresolved. Evidence narrow-dd-aux8-full-retained-gelu-completion-validation.json.

At 2026-10-09 20:05 UTC, retained tanh runner11898 and its parent/time/evaluator processes were absent. No tanh report or saved exit status exists; the interruption cause is unknown. Completed GeLU evidence remains valid. A tanh-only retry, runner55609 (evaluator1865842/start222290519), launched after unchanged source/executable/state/corrected-selector hashes and a fresh memory budget passed. It reuses the same original keys and input; no fresh keys are generated. Independent refresh and square validation remain pending. Packed reconstruction1756575/start211884551 remains live. Evidence: narrow-dd-aux8-retained-tanh-interruption-validation.json, narrow-dd-aux8-retained-tanh-retry-validation.json and narrow-dd-aux8-retained-tanh-retry-live-provenance-validation.json.

Full retained aux8 tanh retry runner55609 completed exit0. Independent N16384 refresh max2.0237842266785237e-8 RMS1.1008497132442232e-8 and Q60/delta48 square max5.39755485425234e-8 RMS2.691988812787784e-8 pass fixed23bit and square bound4.5615794863043815e-7. Saved output hashes and original secret/slot/input identity were checked. Together with completed sigmoid and GeLU, all three corrected retained functions independently pass refresh and square. This is the original-key/input fixture with one corrected selector bit; clean fresh production qualification, archived mismatch origin, bit-error origin and structured-key security remain unresolved. At the time of the tanh completion, original packed reconstruction was still live; it later completed exit0. Evidence: narrow-dd-aux8-full-retained-tanh-completion-validation.json and narrow-dd-aux8-corrected-retained-all-functions-validation.json.

Selector integrity snapshots are available through `Parameter-Selection/python/estimates/narrow_ckks_selector_integrity.py`. Use `snapshot --directory SELECTORS --manifest NEW_MANIFEST --n 512` (or `16384`) immediately after generating a completed aux8 selector cache, then use `verify --directory SELECTORS --manifest MANIFEST` before and after evaluation. Store the manifest outside the selector directory. The tool checks all expected group headers, exact sizes and SHA-256 payload hashes with bounded memory, rejects symlink files and refuses to replace a baseline manifest. Fourteen controls pass, including a same-length single-bit payload error, and snapshot/verification of all 35 groups in the retained encrypted N512 cache passes. A manifest proves identity from its observation time only; it does not validate encryption phases, authenticate the manifest, or explain historical corruption. Evidence: `narrow-dd-selector-integrity-controls-validation.json` and `narrow-dd-selector-integrity-retained-small-validation.json`.

Full corrected retained N16384 aux8 selector cache: all36 group headers, sizes and SHA256 hashes recorded and independently re-read successfully in capped128MiB runner96884, exit0. The corrected group matches the previously phase-checked payload hash. The retained-function launcher now records a baseline by default (or accepts --selector-manifest) and verifies every group before and after each function, with a fresh memory budget after pre-launch scanning. All13 manifest/provenance controls pass. These are observation-time integrity checks; historical bit-error origin and encryption phase correctness are separate. Evidence: narrow-dd-selector-integrity-retained-full-validation.json and narrow-dd-selector-integrity-launch-controls-validation.json.

Tightened the packed identity report validator: both square-error gates now use the exact fixed2^-16 C++ identity threshold; evaluation error is explicitly checked, numeric substitutes for boolean profile flags are rejected, and degree/retained-level path consistency is checked. Both actual completed N512 old/current runs still pass strict profile, reported accuracy, normal exit and cleanup checks. All20 controls PASS. Added Parameter-Selection/python/estimates/narrow_ckks_packed_report_validate.py for N16384 completed report validation, including phase timings, storage geometry, /usr/bin/time terminal exit0 and owned-cache cleanup. Synthetic full-profile controls only exercise the reader. The actual full packed one-shot report subsequently passed strict validation with terminal exit0; see narrow-dd-full-packed-completion-validation.json. Evidence: narrow-dd-packed-report-strict-controls-validation.json. Report checking is distinct from independent ciphertext phase validation and security qualification.

### Reusing the packed benchmark setup

Use `--ll14-packed-bounded-checkpoint new-checkpoint-directory [report.json]`
for the first packed identity test and
`--ll14-packed-bounded-replay checkpoint-directory [report.json]` for subsequent
compatible evaluations. The N512 controls use
`--small-high-packed-bounded-checkpoint` and
`--small-high-packed-bounded-replay` with the same arguments.

The checkpoint retains the dense and structured secrets, input slots and
ciphertext, refresh setup keys, and reconstructed public selector files. Its
directory is private; keep it outside the repository. A completed
`checkpoint.bin` is published before refresh evaluation, so the retained setup
can also be replayed after an evaluation failure. The ring size, ciphertext
moduli, scales, storage width, archive version, and function are checked before
loading the key payload. Plain37 and packed checkpoints have different archive
magic values.

Replay skips dense/structured secret generation, refresh setup key generation,
packed key generation, and selector reconstruction. It reruns all selector and
refresh/square accuracy checks. Reports mark `checkpoint_replayed` and
`refresh_keys_reused`; `packed_setup_seconds` and
`selector_reconstruction_seconds` are zero on replay because those phases were
skipped. `key_setup_seconds` includes checkpoint loading and selector validation.
The post-refresh square diagnostic saves its multiplication key as `square-key.bin`
in the private checkpoint directory and loads it on subsequent runs. Existing
checkpoints acquire this sidecar on their first replay. The loader checks the
profile and dense secret; malformed or incompatible saved keys fail the run.
Reports mark `post_refresh_square_key_reused`. Checkpoint evaluations also save
an adjacent `.square.ciphertext` for comparing square outputs across replays.
Use the original creation report for packed setup/storage measurements; a replay
report does not qualify as a fresh setup performance trial. Checkpoint modes
retain selector files and are therefore separate from the one-shot report
validator's cache-cleanup contract.

The completed original N16384 one-shot test has no reusable checkpoint. These
modes support future runs; they cannot recover keys already discarded by that
process.

Validate a retained packed identity checkpoint independently with:

```sh
python3 ../Parameter-Selection/python/estimates/narrow_ckks_packed_checkpoint_validate.py \
  --checkpoint /private/checkpoint/checkpoint.bin \
  --refresh report.json.ciphertext \
  --square report.json.square.ciphertext \
  --report report.json --output phase-validation.json
```

Run this from the TFHEpp directory. The validator reads only the checkpoint's
metadata, secrets, slots and input prefix. It computes integer negacyclic secret
products and decodes slots independently of TFHEpp for the input at Q85/scale41,
refresh at Q100/scale49, and square at Q51/scale49. All three must pass the fixed
2^-16 identity threshold; the independently calculated metrics must agree with
the C++ report within 1e-11. The prefix hash identifies the inspected bytes;
evaluation-key payloads and selector phases are outside this check.

Both retained N512 creation/replay outputs pass, with independent refresh error
1.1839592284083153e-6 and square error 4.185404876517042e-7. Evidence is in
`narrow-dd-packed-checkpoint-independent-phase-validation.json`. The reader
controls in `narrow_ckks_packed_checkpoint_controls.py` cover little/big-endian
N512 and little-endian N16384 synthetic archives and 16 rejected malformed,
profile, ciphertext and report cases. Synthetic controls establish reader
behavior; they do not establish encrypted full-ring packed refresh accuracy.
The completed original full-ring one-shot run has no saved checkpoint/output
pair for this independent check, and security remains unqualified.

`narrow_ckks_packed_checkpoint_run.py` automates retained setup and replay:

```sh
python3 ../Parameter-Selection/python/estimates/narrow_ckks_packed_checkpoint_run.py \
  --exe build/test/ckks/ckks_narrow_benchmark \
  --n 16384 --root /private/new-packed-run --output run-validation.json
```

The output root must not exist. To reuse a compatible completed checkpoint,
add `--checkpoint /private/existing-checkpoint`; then both evaluations use
replay mode. The runner pins a copy of the benchmark, independently checks
saved input/refresh/square phases after each successful evaluation, compares
refresh and square ciphertext bytes, and hashes retained files before and after
replay. It uses four OpenMP threads, a 64GiB address-space limit at N16384
(4GiB at N512), and retains private artifacts. Provision sufficient RAM and
scratch disk: full packed reconstruction writes about112GiB cumulatively
by updating a fixed accumulator of about32.55GiB. Accumulator and finalized
selector files together occupy at most65.09GiB during finalization; checkpoint,
pinned executable and report files require additional space. File hashes
establish identity from the first completed evaluation; they do not establish
selector generation-time integrity or key security.

The N512 runner control reused the existing checkpoint and passed without
fresh setup. A retained N16384 setup/replay is running; its terminal results
will be recorded in `narrow-dd-full-packed-checkpoint-run-validation.json`.
This remains separate from the completed original one-shot report.

The packed identity profile's security bookkeeping is separate from the
functional auxiliary-switch profiles: its structured switch exposes Q48,
not Q56/Q64/Q76. The archived generic dense screens at Q379 and Q427 report
weakest costs150.01 and131.60bits under the pinned BDGL16 model. The Q48
uniform-sparse proxy has only nonfinite attack-cost errors and remains
incomplete. Raw results are now retained in
`narrow-dd-packed-q48-raw-optimizer-validation.json`: all four attacks return
positive infinity at both search limits1754 and4096. These are incomplete
optimizer outcomes, not finite security bounds. That proxy also fixes
16positive/15negative signs instead of the actual random signs with one secret
per public block. Generic screens exclude ring-specific and correlated-key
attacks. No128-bit structured-key claim follows from these numbers, and the
higher-modulus sensitivity results do not qualify this packed profile.
Evidence: `narrow-dd-packed-actual-key-exposure-validation.json`.

For aux8 prefix diagnosis, `canonical_component_decode` first applies the
real/imaginary polynomial projections with integer arithmetic, then decodes
each real component separately. This preserves a small component when the
other is so large that ordinary complex FFT decoding causes cancellation.
Ten analytical constant/cosine controls at N16/N512/N16384 cover both dynamic
range directions and direct/FFT decoding; retained encrypted N512 component
checks also pass. Aggregate metrics used by existing report validation retain
the ordinary decoder and its original gate semantics.

The archived failed aux8 output has isolated real maximum error
4.165432110347922e-9 (fixed23-bit PASS) and imaginary error
1.7331694525133953e20 (FAIL). The current always-bad-selector replay fails
both isolated components. Ordinary complex decoding had reported a spurious
millions-sized real error from cancellation against the huge imaginary signal.
The archived complex refresh remains failed. Evidence:
`narrow-dd-aux8-archived-component-error-validation.json` and
`narrow-dd-component-phase-decode-controls-validation.json`.

The diagnostic command
`ckks_narrow_auxiliary_replay --selector-halves STATE REAL_SELECTORS IMAGINARY_SELECTORS NEW_DIRECTORY`
loads original retained keys/input and uses the first source for the real
12-root tile and the second for the imaginary12-root tile. It verifies the
complete ordered selector passes and rejects an unexpected load sequence.
This intentionally models different selector states between component passes;
it does not mutate caches. N512 same-source control gives the original
ciphertext bytes with70group loads and unchanged retained files. The full
N16384 corrected-real/original-imaginary replay completed with exit0 and72
ordered group loads. It reproduces the archived output ciphertext exactly.
Independent integer-phase validation also confirms byte-identical saved
prefixes: isolated real error4.165432110347922e-9 passes the fixed23-bit
bound, while imaginary error1.7331694525133953e20 fails. The retained state
and both selector sources are unchanged. No new keys were generated.

Different selector states between the real and imaginary passes are sufficient
to reproduce the archived failure. This does not establish when or how the
single selector bit changed. The diagnostic is successful; encrypted accuracy
is failed. Evidence: `narrow-dd-aux8-selector-halves-control-validation.json`,
`narrow-dd-aux8-selector-halves-full-validation.json`,
`narrow-dd-aux8-selector-halves-independent-phase-validation.json` and
`narrow-dd-aux8-selector-halves-full-completion-time.log`.

CPU performance remains below the paper's implementation. Section 8 uses
SEAL with AVX2/OpenMP; Table 5 reports LL14 real bootstrapping at3.18s
on4cores and410ms on48cores. Section 8.2 excludes one-time packed-key
reconstruction and does not give its duration. The completed TFHEpp packed
N16384 run measured137255.75s of reconstruction and628.36s of refresh
with4threads. Different parameters, arithmetic and hardware prevent a
controlled comparison; these measurements do not reproduce paper latency.
Key reuse removes repeated setup, but does not reduce that measured refresh
latency. The CPU implementation remains experimental.

Uncached packed extraction now uses CKKSEvalAutoStreamedFFT: transform input
digits once, transform key digits on demand, sum primary-row products before
one inverse FFT per secondary row/component, then recompose the torus.
The full quadratic key spectrum cache is not retained. Existing prepared
depths and extraction order are unchanged. The pinned live full checkpoint
job uses the earlier executable and is unaffected by this source change.

Same-key N512 Q422 validation compares all2048 ordered terms (2097152
coefficients) against the pre-change reference and passes with exit0.
Both fixture files retain their original SHA256 hashes. A current raw-path
baseline also passes: extraction70.410s versus57.784s streamed, peak RSS
249080KiB versus249096KiB. These are single sequential observations under
concurrent full-ring load, not a full-ring speedup claim. Evidence:
`narrow-dd-packed-streamed-fft-validation.json`.

Factor evaluation now reconstructs CRT digits directly into accumulators
owned by table/component/column. The chunked LL14 12-column path removes
1056964608bytes of torus digit scratch, replacing it with1048576bytes of
worker inverse spectra at4threads. The prepared path additionally retains
301989888bytes of modular sums for those dimensions; its pointwise batch
sharing remains unchanged. This is payload accounting, not measured full-ring
RSS. Retained N512 replay uses existing refresh/square keys, preserves all37
files, and gives byte-identical refresh/square ciphertexts; independent
integer-phase accuracy passes. Focused N512 one/two-prime tests match the
coefficient reference for1/3columns and all five selector chunk policies.
Replay refresh took3.333s; no controlled timing comparison or full-ring speedup
is established. Hardware sampling with perf was unavailable under the host
setting perf_event_paranoid=4. Evidence: `narrow-dd-direct-recomposition-validation.json`.

`ckks_narrow_factor_checkpoint --retained-aux8 STATE SELECTOR_DIRECTORY`
loads compatible saved auxiliary8 keys instead of generating a conjugation
key. It rejects foreign metadata and reports optional caller-owned factor
stage counters. Full N16384 final head30/tail4 checks pass: chunk40 versus17
and single columns0/5/11 are coefficient-identical,196608 decoded slots
pass2^-23 with maximum1.354e-9, and state/both loaded group hashes are
unchanged. Four-thread first-call per-group timings are about11.7s matrix
multiplication,7.8s selector NTT preparation,2.3s conjugation,0.4s inverse
NTT/CRT/recomposition,0.14s table preparation, and0.34-0.41s selector load.
These timings exclude other checks and key setup, use the retained aux8
profile, and do not reproduce the whole-refresh paper benchmark. Evidence:
`narrow-dd-retained-factor-profile-validation.json`.

A table-spectrum Shoup quotient experiment passed coefficient-reference
controls and both retained full-ring factor groups, but its observed matrix
stages were slower (head13.936s versus11.721s; tail12.594s versus11.749s).
The change was rejected, the exact validated header restored, and the factor
diagnostic rebuilt. These are single trials under concurrent load, not a
statistical throughput comparison. Evidence:
`narrow-dd-table-multiplier-experiment-validation.json`. Integer arithmetic
bounds permit proposed48-bit auxiliary NTT selector limbs for both actual
N16384 Q385/plain37 andQ379/plain34 profiles with the existing two-prime
CRT. That proposal is not yet implemented or runtime-validated; evidence:
`narrow-dd-wide-ntt-limb-feasibility.json`.

Auxiliary selector NTT limbs are now adaptive:16bits belowN4096, otherwise
48bits when the complete worst-case convolution fits the existing two-prime
CRT range, with32/16bit fallbacks. Compile-time range screens use division
to avoid overflowing while rejecting a candidate, and runtime checks bound
the actual encoded tables. Ciphertext moduli/keys are unchanged; NTT caches
are derived and rebuilt after loading. Explicit16-bit baseline evaluation
is available through the optional digit-width template argument.

Retained full N16384 Q385/plain37 aux8 head30/tail4 comparisons pass: all
786432 factor coefficients exactly equal16-bit evaluation; chunk40/17 and
single0/5/11 comparisons also pass, and196608 slots meet2^-23. No full keys
were generated, and state/loaded group hashes are unchanged. Paired matrix
timings improve12.122to4.157s and11.745to4.093s; selector preparation improves
7.870to2.826s and7.724to2.821s. These are sampled factor-stage observations,
not whole-refresh/paper performance claims. N512 forced48-bit and default
16-bit coefficient-reference controls pass for both prime counts and five
chunk sizes. A valid N65536/plain50/h4 configuration statically exercises
the32-bit fallback. Evidence: `narrow-dd-wide-ntt-validation.json`.

The optional factor observer reports completed public head/tail groups and
accumulated factor-stage seconds. It is excluded from key archives and
cleared on load. N512 retained replay with observation preserves the saved
ciphertext and complete unobserved prefix, reports all70groups, and passes
accuracy without newkeys. A vector-to-span template conversion issue found
by the full replay build was corrected by making the prepared-cache span
a non-deduced parameter; explicit digit-width overrides remain available.
Evidence: `narrow-dd-factor-observer-control-validation.json`.

One full N16384 optimized retained sigmoid replay is running with the original
keys/input and the previously verified corrected selector fixture. Its runner
requires equality to the prior passing prefix, unchanged state/all36selector
groups/baseline prefix, and independent integer-phase validation. It skips
key generation and packed reconstruction. Whole-refresh performance remains
unmeasured until actual completion. Evidence:
`narrow-dd-wide-ntt-full-refresh-validation.json`.

`makeCKKSCompactDDRelinKeyWithSecondaryBase<NewBbar>` changes only the
public compact DD secondary representation. It reconstructs each encrypted
primary row moduloQ, then decomposes it at the new base. It does not need
a secret or encrypt new rows, and rejects out-of-range source digits. The
original key stays unchanged. N512 controls atQ17/85/379/448 from4-bit
secondary digits, plusQ97 from12/24-bit digits, all pass16-bit conversion
round trips and coefficient-identical scalar/FFT switches across3random
inputs per case. The production compact DD refresh now prepares its derived
FFT chain using16-bit secondary digits when the existing mantissa guard and
storage width permit; otherwise it retains the original width. Serialized
keys and primary encrypted rows are unchanged. Conversion prepares one key
level at a time. Full-ring performance remains unverified. Evidence:
`narrow-dd-secondary-representation-validation.json`.

The retained N512 full sigmoid replay with the converted chain preserves
the prior ciphertext and complete output prefix exactly. All36input files
are unchanged, no keys are generated, and independent integer-phase decoding
passes the fixed23-bit gate. Its33.99s refresh versus the earlier33.87s
control does not establish a whole-refresh speedup. The full-ring replay
already running is pinned to the earlier executable and cannot validate this
cache change. Evidence: `narrow-dd-secondary-cache-replay-validation.json`
and `narrow-dd-secondary-cache-independent-phase-validation.json`.

Bounded tail alignment now gathers all tail leaves in a frequency tile before
calling the binary rotation batch once. Each ciphertext still receives its
same ascending-bit chain and auxiliary-scale rounding. The staged factor
conjugation cache is released before alignment; only one derived rotation
key is resident. Original tail vectors already retained all leaves, so the
additional allocation is only the public rotation-step vector. Structural
key-preparation counts fall from12to3 per N512 tile and28to6 per N16384
tile. These counts do not imply a measured latency ratio.

The retained N512 full sigmoid replay matches the previous ciphertext and
complete prefix exactly, leaves all36input files unchanged, generates no
keys, and passes independent23-bit phase decoding. Its42.17s timing overlaps
a build and two running full benchmarks; speedup and full-ring correctness
of this change remain unverified. Both replay and broader bootstrapping
targets build successfully. Evidence: `narrow-dd-tail-shared-replay-validation.json`
and `narrow-dd-tail-shared-independent-phase-validation.json`.

`ckks_narrow_auxiliary_replay --dd-cache-control STATE SELECTORS NEW_DIRECTORY`
compares4/16-bit derived DD FFT caches using the same retained encrypted keys.
At N16384, bothQ348 andQ156 pass exact secondary-conversion round trips and
three switch probes (uniform polynomial, centered boundaries, negative impulse).
All196608switch coefficients match; the original6.45GBstate hash is unchanged,
and no encryption or new keys are used. The control exits0 under a16GiBcap.

AtQ348, the FFT payload falls from1,984,167,936to501,743,616bytes; three
switches take1.247s versus0.434s. AtQ156, payload falls from398,721,024to
102,236,160bytes; switches take0.245s versus0.086s. Conversion itself costs
2.891s/0.573s, and new FFT preparation adds0.407s/0.082s, so one-time cache
preparation is slower than simply preparing the old key. These single
component observations under concurrent CPU load do not establish complete
refresh latency. Intermediate levels, complete converted-chain refresh, and
shared full-ring tail alignment remain pending. Evidence:
`narrow-dd-secondary-full-ring-validation.json`.

`--tail-cache-control` computes a frequency1 factor from retained tail0
selectors, then compares per-pack and shared alignment for31leaves in each
of three distinct public column multiples. At N16384 all93ciphertexts
(3,047,424coefficients) match exactly, including the partial final pack.
The retained state and consumed selector group hashes are unchanged; no
keys or encryption are generated. Both methods retain a single2,076,442,624
byte rotation-key payload. Alignment takes129.45s versus94.57s in one paired
observation under concurrent CPU load. This compares the alignment stage,
not the complete refresh; factor generation and archive loading are excluded.
N512 controls also pass. Evidence: `narrow-dd-tail-full-ring-validation.json`.

The combined current-source N16384 retained sigmoid refresh uses adaptive
48-bit selector NTT limbs, the16-bit derived DD chain, and shared tail
alignment. It generates no keys or input and consumes the previously
verified corrected selector fixture. The runner pins its executable, checks
future memory growth against existing job caps plus a10GiBhost reserve,
and requires equality to the previous passing prefix, unchanged hashes of
all38retained files, and independent23-bit phase validation before success.
The earlier wide-NTT full replay is a separate pinned baseline. Combined
whole-refresh accuracy and latency remain pending. Evidence:
`narrow-dd-combined-cpu-full-refresh-validation.json`.

Prepared automorphism FFT caches accept an explicit secondary digit width:
`makeCKKSPreparedAutoKeyFFT<P,LogQ,16>(raw)`. The primary raw gadget rows
and offsets stay unchanged. Secondary digits use their own row count and
physical shifts, including unalignedQ, with an FFT accumulation guard and
storage-divisibility check. The default now chooses16bits when the storage
width and FFT mantissa guard permit, otherwise the originalP::Bbar width.
SpecifyP::Bbar explicitly to prepare the original baseline representation.

At N512 Q422, three retained packed substitution keys (depths0/3/8) and
four encrypted indicator inputs each produce coefficient-identical outputs
with4/16-bit caches: all12ciphertexts/12,288coefficients pass. The original
public-key fixture hash is unchanged; no keys or encryption are generated.
Cache payload falls from92,045,312to23,445,504bytes per key. These controls
do not establish full-ring accuracy, broader parameter compatibility, or
performance. Evidence: `narrow-dd-auto-secondary-representation-validation.json`.

Additional N512 auto-cache controls publicly reduce the retained Q422
primary rows toQ17/85/356/416, preserving their physical gadget shifts.
At each reduced modulus, four encrypted indicator inputs compare both
4/16-bit FFT formats to raw switching exactly. Combined with the original
three Q422 depths, all28ciphertexts/28,672coefficients pass; the saved public
key hash stays unchanged. No keys or encryption are generated. Low-modulus
projections test arithmetic equality, not cryptographic accuracy. The
repository's existing prepared-auto test also now checks explicit16-bit
caches against its raw reference; broader target execution remains pending.
Evidence: `narrow-dd-auto-secondary-modulus-validation.json`.

The earlier full N16384 adaptive48-bit NTT sigmoid replay has completed and
passed all gates: exact equality to the prior passing corrected-selector
prefix, unchanged hashes of all38retained files, and independent23-bit
phase decoding. It generates no keys/input and performs no packed-key
reconstruction. Measured whole functional refresh latency is3218.265s
(53.64minutes), with72factor groups. This is the12-frequency complex sigmoid
workload, not the single identity refresh or the paper's real-input benchmark.
It uses the original4-bit DD FFT chain and per-pack tail alignment; the
combined newer-cache replay is separate and remains running. Source and
executable provenance are pinned in
`narrow-dd-wide-ntt-full-refresh-validation.json`.

At N16384, the retained factor/output conjugation keys atQ356/116 pass
secondary4/16-bit comparisons on a lifted saved input, uniform arithmetic
coefficients, and signed boundaries: all196,608coefficients match. No keys
or encryption are generated, and the retained state hash is unchanged.
AtQ356, cache payload falls from2,076,442,624to536,608,768bytes, and three
switches take0.937s versus0.385s in a component observation. Evidence:
`narrow-dd-auto-full-ring-validation.json`.

Adaptive secondary16 auto FFT caches are now enabled in production where
the guard permits, with unchanged primary keys and archives. Both the full
replay and broader bootstrapping targets build. A complete retained N512
sigmoid replay preserves the previous ciphertext and prefix exactly, leaves
all36input files unchanged, and passes independent23-bit phase decoding.
Its42.88s timing overlaps a build and other benchmarks, so it does not
establish whole-refresh speedup. The running combined N16384 replay remains
pinned to its earlier auto4 executable; full current-source auto16 refresh
validation is pending. Evidence: `narrow-dd-auto-integrated-replay-validation.json`.

The current-source full N16384 sigmoid replay with all three CPU cache
changes is now running: adaptive48-bit selector NTT limbs, secondary16 DD
and auto FFT caches, and shared tail alignment. It retains the original
keys/input and verified corrected selectors, under a20GiBaddress-space cap.
Its runner pins source/executable provenance and requires exact equality to
the known passing prefix, unchanged hashes for all38retained files, and
independent23-bit phase decoding before success. No whole-refresh gain is
claimed until completion. The earlier auto4 combined run remains separate.
Evidence: `narrow-dd-adaptive-auto-full-refresh-validation.json`.

Unsigned plaintext polynomial multiplication now selects16-bit computational
digits when the digit FFT backend, storage divisibility and mantissa guard
permit, with the original gadget width as fallback. Explicit fourth template
arguments to `PolyMulTorusByUnsignedBits` select the baseline or wider radix.
Encrypted gadget rows, secrets and archives are unchanged. Plaintext input
coverage retains the original gadget-rounded bit boundary, including unused
upper bits, so widening cannot include coefficients previously ignored.

At N16384, independent arbitrary-precision convolution checks signed sparse
products atQ81/93/413/448, zero, full-width unsigned and maximum-amplitude
inputs, plus legacy plaintext coverage17/20/37bits. All163,840oracle
coefficients pass; another49,152dense coefficients match the4-bit baseline
atQ93/413/448. Sample component timings improve0.116to0.012s,0.516to0.052s,
and0.820to0.093s under concurrent CPU load. These are component observations,
not whole-refresh or paper-performance claims. No keys or encryption are
generated. Evidence: `narrow-dd-plaintext-radix-controls-validation.json`.

The default full-ring plaintext oracle and all three affected build targets
pass. A complete retained N512 sigmoid replay preserves the saved ciphertext
and entire prefix exactly, leaves36input files unchanged, and passes
independent23-bit phase decoding. It takes34.98s under concurrent load;
whole-refresh speedup is unestablished. Both running functional N16384
replays are pinned to executables preceding this plaintext-radix change, so
full current-source refresh validation remains pending. Evidence:
`narrow-dd-plaintext-radix-integrated-validation.json`.

Monitoring scripts update shared status files and must be called sequentially.
A concurrent-write trailing fragment in the monitor JSON was repaired by
retaining its complete object and revalidating both live process identities;
computation and cryptographic artifacts were unaffected.

The zero-rescale ciphertext tensor product now uses balanced 16-bit arithmetic
digits when the multilimb FFT precision bound permits. This radix is independent
of the encrypted DD gadget and changes neither serialized keys nor key generation.
Changing the integer representatives only adds multiples of the input moduli, so
these differences vanish modulo the smaller input modulus. Nonzero rescaling
retains the original gadget radix because its balanced-digit wrap threshold
controls signed rounding ties.

The full-ring control compared the original 4-bit and explicit 16-bit paths with
an independent exact integer NTT oracle across 18 random and boundary cases
(884,736 coefficients). The production oracle then checked both zero and nonzero
rescaling across full and packed modulus schedules, including dense maximum
16-bit digits at Q448 (2,359,296 coefficients). Both controls passed; see
`narrow-dd-tensor-radix-controls-validation.json` and
`narrow-dd-tensor-radix-production-validation.json`. Component timings were taken
under concurrent CPU jobs and do not establish whole-refresh speedup. Running
full-ring refresh binaries predate this change; full current-source validation
remains pending.

The complete retained N512 sigmoid refresh also passed with identical saved
ciphertext prefix, all 36 original input files unchanged, and an independent
fixed 23-bit phase gate. No keys were generated. This run included both the
plaintext radix16 optimization and the zero-rescale tensor radix16 change; see
`narrow-dd-tensor-radix-integrated-validation.json`. Its 40.24-second timing under
concurrent builds and full jobs does not establish whole-refresh speedup.

The full N16384 retained refresh combining secondary16 DD preparation and shared
tail rotations completed with actual exit0, identical passing ciphertext prefix,
38 unchanged retained files, and independent 23-bit phase checks. It took
3221.11 seconds (53.7 minutes), versus 3218.27 seconds for the earlier wide-NTT
run. These concurrent runs do not establish a whole-refresh speedup. Its pinned
binary predates adaptive16 auto preparation and both plaintext/tensor radix
changes; see `narrow-dd-combined-cpu-full-refresh-validation.json`.

Multilimb tensor rescaling now extracts limbs directly for division by one and
uses signed magnitude shifts for division by a power of two. This replaces the
bitwise general integer divider and preserves ties away from zero. The full-ring
independent NTT checks passed for 2,359,296 coefficients, and 2,992 independent
cpp_int accumulator cases cover drops 0/1/48/49/63/64/65/127/128/256/447. The latter
are reproducible with `ckks_narrow_tensor_exact --rescale`. See
`narrow-dd-tensor-shift-controls-validation.json`.

The complete saved-key N512 sigmoid replay passed with an identical output prefix,
36 unchanged original inputs, and an independent fixed 23-bit phase check; see
`narrow-dd-tensor-shift-integrated-validation.json`. It took 20.41 seconds under
concurrent CPU jobs, which does not establish a controlled whole-refresh speedup.
The ongoing full-ring current-radix binary predates the power-of-two shift change;
full current-source refresh validation is still pending. No new keys were generated
for these tests.

The zero-rescale multilimb tensor path now accumulates modulo the torus storage
width and skips digit products whose shifts are at least the active output Q.
These products are exact multiples of the output modulus. At Q348 this retains
253 of 484 radix16 digit pairs; the accumulator uses 7 rather than 16 limbs.
Nonzero rescaling retains the full signed product and wide accumulator.

Independent full-ring NTT checks passed for 2,359,296 tensor coefficients; the
N16 scalar fallback independently passed 2,304 coefficients. A complete retained
N512 sigmoid refresh preserved the exact saved prefix, all 36 original input
files, and an independent 23-bit phase gate without generating keys. Evidence:
`narrow-dd-modular-tensor-controls-validation.json`,
`narrow-dd-modular-tensor-fallback-validation.json`, and
`narrow-dd-modular-tensor-integrated-validation.json`. Component timing probes are
concurrent observations, not a controlled whole-refresh speedup measurement.
The running full current-radix executable predates both the power-of-two shift
and modular accumulation changes; full latest-source refresh remains pending.

Uncached streamed packed substitutions now use the same adaptive secondary radix
as the prepared auto cache (16 bits when the FFT precision bound allows), while
retaining the original primary gadget. Raw key coefficients are decomposed on
demand; serialization and encrypted keys are unchanged. Explicit original-radix
streaming remains available as a comparison template argument.

N16384 synthetic raw-coefficient controls at Q422/Q356/Q17 compared prepared4,
streamed4, and streamed16 for random, sign/wrap, and dense maximum-digit inputs:
589,824 coefficients matched. These are arithmetic controls, not encrypted packed
accuracy evidence. The retained public N512 packed-key fixture then reconstructed
all 2,048 ordered terms identically (2,097,152 coefficients), with both fixture
files unchanged and no new encrypted keys. That concurrent extraction took
15.75 seconds; individual large-Q streamed calls were approximately three times
faster, without a full-ring reconstruction speedup claim. See
`narrow-dd-streamed-secondary16-validation.json`. The long retained packed run
uses its original pinned executable and does not validate this newer path.

Source inspection also found the fused BSGS plaintext kernel still uses four-bit
arithmetic digits independently of the optimized standalone plaintext multiplier.
The proposed independent arithmetic radix and necessary exact group controls are
recorded in `narrow-dd-fused-bsgs-radix-feasibility.json`; no BSGS change or gain
is claimed yet.

The earlier full N16384 adaptive auto16 refresh completed with exit0, exact
passing prefix, 38 unchanged retained files, and an independent 23-bit phase
check. Its runtime was 2926.96 seconds (48.8 minutes). This pinned binary predates
the plaintext, tensor, and streamed-secondary changes; see
`narrow-dd-adaptive-auto-full-refresh-validation.json`. Concurrent timings do not
constitute a controlled matched speedup or paper-performance reproduction.

The fused multilimb BSGS kernel now uses independent arithmetic digits (adaptive
16-bit by default), retains only the rows covering the active Q, and omits
products divisible by 2^Q. Its existing active-Q reduction precedes signed
rounding, so representative changes and omitted high terms cannot alter ties.
The full signed-magnitude uint64 plaintext range and bounded Fourier batching
are preserved. Rotation keys, plans, serialized gadgets, and group order are
unchanged. An explicit final template argument can select the comparison radix.

The frozen original4/full-width kernel and revised kernel both matched an
independent cpp_int negacyclic convolution plus signed-rounding oracle at Q81,
Q93, Q413, and Q448. Controls included 33-term groups spanning Fourier batches,
signed sparse full-uint64 magnitudes, dense maximum magnitudes, and sign/wrap/tie
boundaries: 786,432 coefficients across both paths passed. At Q93 the revised
kernel uses six torus rows and four plaintext chunks versus the original 112
and 16. The 33-term sparse Q93 sample took 0.067 seconds versus 1.037 seconds;
these concurrent component observations do not establish whole-refresh gain.
See `narrow-dd-fused-bsgs-controls-validation.json`.

A complete retained N512 sigmoid refresh also preserved its exact passing prefix,
36 unchanged original input files, and independent 23-bit phase accuracy, without
new keys; see `narrow-dd-fused-bsgs-integrated-validation.json`. The latest-source
full-ring runner reuses the retained corrected-selector fixture and enables the
existing public linear-stage timing logs. Full latest-source refresh validation
remains pending in `narrow-dd-fused-bsgs-full-refresh-validation.json`.

The fused BSGS independent oracle is now a repository test:
`cmake --build TFHEpp/build --target ckks_fused_bsgs_exact`, followed by
`TFHEpp/build/test/ckks/ckks_fused_bsgs_exact --all-rings`. It checks native FFT
rings N8192, N16384, and N32768 with 448-bit coefficient storage, using direct
cpp_int convolution and signed rounding. All 36 groups and 1,376,256 coefficients
passed. The default command checks N16384 only. No keys or ciphertext encryption
are needed; see `narrow-dd-fused-bsgs-all-rings-oracle-validation.json`. These are
arithmetic backend checks, not complete encrypted parameter-profile validation.

The live latest-source full N16384 refresh completed both S2C transforms, logging
51.321 and 45.965 seconds (97.286 seconds together), with 128 finalized giant
groups per transform. This is stage timing evidence, not terminal refresh
accuracy or whole-refresh latency; see
`narrow-dd-fused-bsgs-live-s2c-timing-validation.json`. The full refresh continues
with retained keys and independently validates its output only after completion.

The current implementation also squared the known passing full N16384 sigmoid
ciphertext using its original retained compact DD square key. The Q60/delta48
ciphertext matched the prior passing square byte-for-byte, and an independent
integer phase decoder confirmed max error 9.21697e-9 within the propagated
23-bit product bound. State, prefix, and prior square files remained unchanged;
no keys were generated. See
`narrow-dd-current-retained-square-control-validation.json`. This validates
current retained multiplication on a previously passing output; it does not
claim the running latest refresh has completed.

A separate latest-output square runner is prepared at
`/home/work-bleaker/kmatsuoka/ckks-diagnostics/fused-bsgs-post-refresh-square-20261010/run.py`.
It requires terminal successful refresh validation, identical passing prefix,
unchanged retained inputs, and an independent 23-bit refresh gate before it can
square the latest output. Its validation file is
`narrow-dd-fused-bsgs-post-refresh-square-validation.json`; the runner has not
started while the full refresh is still live.

The earlier full N16384 current-radix refresh completed with exit0, an identical
passing output prefix, 38 unchanged retained files, and independent 23-bit
accuracy. It took 2381.946 seconds (39.7 minutes). Its pinned executable includes
adaptive plaintext16 and zero-rescale tensor16, but predates power-of-two tensor
rescaling, modular tensor pruning, and fused active-Q BSGS arithmetic. Thus it
qualifies that earlier profile, not the whole latest source. See
`narrow-dd-current-radix-full-refresh-validation.json`. Concurrent timings are
not a controlled matched speedup or paper-performance reproduction.

The latest full N16384 sigmoid refresh completed with exit0 and all gates passing:
identical known passing prefix, 38 unchanged retained files, and independent
23-bit accuracy (max4.17920e-9). The 1274.871-second runtime is 21.25 minutes.
The actual newly completed output was then squared with its original retained
key; independent phase accuracy and exact prior square ciphertext both passed.
No keys were generated. Evidence:
`narrow-dd-fused-bsgs-full-refresh-validation.json` and
`narrow-dd-fused-bsgs-post-refresh-square-validation.json`. These are the retained
corrected-selector fixture; original selector-bit error cause and structured-key
security remain unqualified. Concurrent run timings do not reproduce the paper's
matched CPU latency.

Current-source N512 sigmoid, GeLU, and tanh each passed independent refresh and
retained-square checks using unchanged original keys/input/selectors. Evidence:
`narrow-dd-current-functions-small-validation.json`. The current full N16384
GeLU/tanh runner evaluates those functions sequentially with a 16GiB child cap,
original keys and selectors, independent output/square phases, original-input
identity, and selector integrity checks. Its first evaluator is GeLU; full
function completion remains pending in
`narrow-dd-current-functions-full-validation.json`.

The paired signed-table matrix experiment was not adopted. All public modular
outputs matched, but three interleaved trials showed mixed timings (including a
19.5% slower median for the current 12-column shape). These concurrent timings
are insufficient to justify changing the production loop. Evidence:
`narrow-dd-paired-table-matrix-controls-validation.json`.

Selector preparation now retains immutable public NTT plans per schedule and
digit width. Previously each selector/prime constructed identical twiddle and
twist tables; the shared plans preserve all transforms and selector formats.
No secret material or ciphertexts are cached by this change. Native OpenMP
controls compared every old/new spectrum and inverted every row against direct
active-Q limb extraction: N512 with one/two primes and 16/48-bit digits, and
N16384 with two primes, 48-bit digits, and 40 selectors. All twelve trials
passed; the full-ring component medians were 0.2371s before and 0.1741s after
under concurrent load. This is not an end-to-end or paper speedup claim.
Evidence: `narrow-dd-shared-selector-plans-validation.json`. Existing full
GeLU/tanh and packed-checkpoint jobs use their earlier executables and therefore
do not validate this new plan-sharing change.

The shared-plan source subsequently passed an actual N512 retained encrypted
sigmoid replay: exact passing prefix SHA256, independent 23-bit phase gate,
36 unchanged retained files, and no fresh keys. The target build and both
replay/independent decoder exited0. See
`narrow-dd-shared-selector-plans-small-validation.json`. Full-ring encrypted
validation of this change remains pending; the earlier full refresh evidence
applies to its pinned executable.

A pinned current-source full N16384 shared-selector-plan replay is now running
in `/home/work-bleaker/kmatsuoka/ckks-diagnostics/shared-selector-plans-full-refresh-20261010`.
It reuses the same corrected-selector fixture and original retained keys/input,
compares the full output prefix with the known passing hash, checks all retained
files, and independently decodes the output after exit0. Its 20GiB cap passed a
memory check reserving the packed job's remaining address-space headroom, the
full 16GiB sequential-functions cap, and 10GiB host reserve. Evidence:
`narrow-dd-shared-selector-plans-full-refresh-validation.json`.
Until all gates complete, this is a pending validation; concurrent runtime
cannot establish a controlled speedup or a matched paper benchmark.

The pinned full N16384 GeLU function completed with exit0: refresh1069.223s
(17.8min), independent max error2.02489e-8 within the fixed 23-bit gate, and
independent retained-square max error1.55838e-7 within the propagated product
bound. The original secrets, slots and input matched the retained fixture;
selector payloads matched their pre-run integrity snapshot. No keys were
generated. Tanh has started next in the same sequential controller. Evidence:
`narrow-dd-current-functions-full-validation.json`. This is the retained
corrected-selector fixture and the pinned executable before shared NTT plans;
the selector-bit error cause and actual structured-key security remain
unqualified, and concurrent timings do not reproduce paper performance.

The frequency-basis benchmark now has checkpoint/replay modes for both its
standard and lean profiles at N512 and N32768. A completed checkpoint retains
original setup keys, basis-reconstruction keys, one shared square key, secrets
for diagnostic decryption, and the encrypted input. Its directory is private;
portable metadata rejects incompatible profiles before loading key payloads.
File selectors remain attached through a validated existing cache on replay.
The square validation key is generated once for all requested functions.
Per-function refresh, square and normalized ciphertexts are saved for exact
replay comparison and independent phase checks. Build/replay validation of
these new modes remains pending until recorded in
`narrow-dd-basis-checkpoint-controls-validation.json`.

The independent basis checkpoint decoder uses fourteen-bit limbs for N32768
with a checked signed-convolution bound under its 31-bit scratch prime. Two
public sparse-ternary controls (random residues and maximal residues) matched
direct Python integer convolution at every coefficient, 65536 checks total.
Evidence: `narrow-dd-basis-wide-phase-decoder-controls-validation.json`.
This validates the arithmetic decoder, not a full encrypted basis refresh.

For example, create one reusable small lean-basis fixture and then replay it:

```sh
TFHEpp/build/test/ckks/ckks_narrow_benchmark --small-lean-basis-fbs-checkpoint all /absolute/new-basis-checkpoint /absolute/initial-report.json
TFHEpp/build/test/ckks/ckks_narrow_benchmark --small-lean-basis-fbs-replay all /absolute/new-basis-checkpoint /absolute/replay-report.json
```

The report path is required for these modes so the public output ciphertexts
can be retained beside it. The checkpoint is a private diagnostic artifact
containing validation secrets. The analogous full-ring mode drops `small-`;
full-ring validation of the new reuse path remains pending. The initial small
fixture and its first replay passed all three C++ accuracy gates and independent
integer refresh/square/normalization gates with identical output bytes and
unchanged checkpoint files. The negative-profile control initially found an
uncaught exception; the benchmark now catches basis errors and returns2.
Final corrected-build rejection/replay validation is pending.

Final frequency-basis checkpoint controls passed with the corrected executable:
initial fixture generation once, two encrypted all-function replays without
fresh keys, exact serialized output/square/normalized ciphertexts for sigmoid,
GeLU and tanh, independent integer phase accuracy, and unchanged checkpoint
files. Incompatible profiles and an existing creation directory both return2
with a clear error and preserve all checkpoint files. Both builds and the
corrected controller completed with exit0. Evidence:
`narrow-dd-basis-checkpoint-controls-validation.json`. This validates current
N512 checkpoint reuse; latest encrypted N32768 basis validation remains pending.

`ckks_fused_bsgs_exact --basis-rings` now exercises N32768 native 576- and
640-bit storage. Its eighteen public arithmetic cases include the basis S2C
Q129/drop37 level, high levels518/536 and full storage576/640. Sparse signed
full-uint64 plaintexts, dense maximal magnitudes and signed halfway ties all
matched independent cpp_int convolution/rounding:1179648 coefficients. No
keys were generated. Evidence:
`narrow-dd-basis-storage-fused-bsgs-oracle-validation.json`. This does not replace
full encrypted frequency-basis validation.

The compiled current lean-basis storage reduces the old conservative envelope
to39206846704bytes; a40GiB process cap retains margin. The launch check currently
fails because existing packed, tanh and shared-plan jobs reserve more remaining
headroom than available memory. No full basis job was launched; retained runs
were not restarted. Evidence:
`narrow-dd-current-basis-checkpoint-memory-budget.json`. Disk has room for the
179721732096-byte selector payload, but memory remains the launch constraint.

The retained replay helper adds `--replay-prepared` to retain derived factor,
output-conjugation and DD-tree FFT caches across frequency tiles. Raw keys,
selector chunks, bounded tail rotations, arithmetic and output format remain
compatible. It reports cold cache preparation separately. Its N512 retained
sigmoid replay matched the passing prefix exactly and passed independent
23-bit accuracy with36 unchanged files and no fresh keys. Evidence:
`narrow-dd-persistent-fft-small-refresh-validation.json`. Concurrent timing did
not establish a speedup; it is an optional diagnostic mode. A full retained
runner is prepared at
`/home/work-bleaker/kmatsuoka/ckks-diagnostics/persistent-fft-full-refresh-20261010/run.py`;
it requires the shared-plan full replay to finish every gate before launch and
rechecks its memory budget. Full persistent-cache accuracy remains pending.

The full N16384 shared-selector-plan sigmoid replay completed with exit0:
1540.531seconds (25.7min), identical known passing prefix,38 unchanged retained
files, and independent23-bit accuracy (max4.17926e-9). Evidence:
`narrow-dd-shared-selector-plans-full-refresh-validation.json`. It uses the
retained corrected-selector fixture; the original error cause and structured
key security remain unqualified. Concurrent timing is not a controlled speedup
or matched paper benchmark. All72 factor groups completed before the final
tree/output work, so that counter alone must not be treated as job completion.

`Parameter-Selection/python/estimates/narrow_ckks_basis_checkpoint_run.py` now
controls a full lean-basis N32768 checkpoint and all-function replay. It requires
validated small reuse and wide storage arithmetic, freezes the validated
executable, checks a40GiB child/512MiB parent envelope and210GiB disk reservation,
and preserves reserved caps until each earlier controller finishes. Initial
setup is generated once; replay must use the same checkpoint and produce exact
refresh/square/normalized ciphertext bytes. Each phase is independently decoded
against its function target and product bound; checkpoint/selector hashes must
remain unchanged. A live controller at
`/home/work-bleaker/kmatsuoka/ckks-diagnostics/basis-full-current-checkpoint-20261010`
is waiting for the tanh controller and memory eligibility. No new keys have
been generated yet. Evidence:
`narrow-dd-basis-full-current-checkpoint-validation.json`. Latest full encrypted
basis validation remains pending. Avoid launching other large jobs outside its
reservation list while this automatic controller waits.

Retained basis setup is published as `checkpoint.bin` before any function
evaluation. An evaluation failure therefore need not require key regeneration:
keep that directory and use the corresponding `--[small-][lean-]basis-fbs-replay`
mode with a new report location after diagnosing the failure. Replay still needs
its normal accuracy gates; a retained file is not itself proof of correctness.

`Parameter-Selection/python/estimates/narrow_ckks_basis_checkpoint_inspect.py`
checks the checkpoint's portable public profile and every expected selector
header/fixed length. `--hash-payloads` records all checkpoint and selector file
hashes; `--reference` checks a prior manifest without modifying inputs. The
inspector does not decode or print secrets and does not validate evaluation-key
suffixes or ciphertext phases. For example:

```sh
python3 Parameter-Selection/python/estimates/narrow_ckks_basis_checkpoint_inspect.py --directory /absolute/checkpoint --hash-payloads --output /absolute/new-integrity.json
python3 Parameter-Selection/python/estimates/narrow_ckks_basis_checkpoint_inspect.py --directory /absolute/checkpoint --reference /absolute/new-integrity.json --output /absolute/new-recheck.json
```

Actual retained N512 controls verified all35 selector groups and the complete
byte manifest. A private copied selector with a payload-bit flip was rejected
against the manifest; corrupting its header was rejected without hashing.
Original files were independently rechecked unchanged. No keys were generated.
Evidence: `narrow-dd-basis-checkpoint-inspector-controls-validation.json`.
Hashes describe the observation time and do not establish generation-time
correctness or the cause of a pre-existing error.


Packed primary-row projection diagnostic (2026-10-10):
`ckks_narrow_packed_primary_projection FIXTURE` reads the retained N512
`packed.bin` / `terms.bin` fixture and selects 53 existing automorphism rows
from the original 106 rows, corresponding to primary decomposition bits 8
instead of 4. Indicators are unchanged. This diagnostic does not generate
or serialize encryption keys; it reconstructs the fixture's deterministic
diagnostic dense secret only for decryption comparisons.

All 2048 ordered selectors were checked against the reference by decoding
ciphertext differences. Maximum slot difference was 2.030767446803211e-10,
below 2^-23. Sixteen evenly spaced terms also passed independent direct
`cpp_int` negacyclic phase calculations with zero coefficient discrepancy.
Both public fixture SHA-256 hashes remained unchanged. Preparation took
0.0181 s; extraction plus validation took 12.1577 s under concurrent load.
The combined time is not an extraction benchmark or evidence of a speedup.
Evidence: `narrow-dd-packed-primary8-projection-validation.json`.

This changes primary decomposition and its noise, unlike adaptive secondary
FFT decomposition. Production parameters are unchanged: adoption requires
full-ring encrypted accuracy and controlled timing evidence. No paper latency
or security claim follows from this small-ring diagnostic.

The projection diagnostic also accepts `--baseline` to evaluate the unchanged
primary4 retained key through the same validation callbacks. Each callback's
elapsed time (including reference reads, decryption, and the sampled integer
oracle) is accumulated separately; extraction time excludes that total.

Four sequential trials in order primary4, primary8, primary8, primary4 passed
all accuracy and phase checks. Baseline ciphertexts matched the archived
reference exactly. Projected ciphertexts differed, with the same maximum slot
difference 2.030767446803211e-10 in both trials. Median extraction times excluding
callbacks were 13.4347 s for primary4 and 10.6951 s for primary8 (observed ratio
1.256). No keys were generated; both fixture hashes remained unchanged.
Evidence: `narrow-dd-packed-primary8-timing-validation.json`.

These are exploratory N512 measurements during concurrent full-ring work,
with only two samples per variant. They support a full-ring experiment when
its retained packed key becomes available, but do not establish an isolated
speedup, predict full setup time, or qualify a production parameter change.


Reusable packed-selector mask FFTs (2026-10-10):
`include/ckks/prepared_plain.hpp` prepares signed magnitude digits of a public
small plaintext once. Each negative coefficient contributes negative unsigned
magnitude chunks. Multiplication therefore uses one torus decomposition and
one convolution per nonzero plaintext chunk, rather than separate positive
and negative products. The same active-Q modular arithmetic and signed rescale
are preserved. `CKKSNarrowExtractSelectorTerms` lazily caches each distinct
head/tail mask; non-FFT parameter backends retain the original path. Encryption
key parameters, primary gadget rows, and serialized key formats are unchanged.

`ckks_prepared_plain_exact` checked seven N16384 cases against independent
`cpp_int` negacyclic products: sparse mixed signed 64-bit magnitudes at Q81,
Q413 and Q448; a zero plaintext; dense maximum unsigned magnitude at Q413;
and alternating maximum signed magnitudes at Q413 and Q448. All 114688
coefficients passed. Retained N512 packed extraction matched all 2048 ordered
terms / 2097152 ciphertext coefficients exactly. The primary8 projection
experiment also still passed its decrypted differences and sixteen direct
integer phase checks. No fresh encryption keys were generated, and public
fixture hashes stayed unchanged.

Four old/new timing runs with primary4 in order old,new,new,old gave median
extraction times excluding validation callbacks of 12.8963 s and 10.6796 s
(observed ratio 1.208). These measurements are exploratory under concurrent
load, not an isolated or full-ring speedup. Evidence:
`narrow-dd-prepared-packed-masks-validation.json`. Full encrypted packed
reconstruction with this updated header remains unvalidated. Existing large
jobs execute pinned binaries that predate this mask cache.

Inspection also confirmed that the packed benchmark currently discards its
raw packed key after reconstruction; its checkpoint retains expanded selectors
and refresh inputs/keys. It can replay refresh without key regeneration, but
cannot rerun reconstruction from that checkpoint. Future reconstruction
experiments need a retained raw packed archive. The running old full packed
binary has not been changed or restarted.


Raw packed-key retention (2026-10-10):
`include/ckks/narrow_packed_file.hpp` provides
`CKKSNarrowSavePackedKey<Schedule>(packed, path)` and
`CKKSNarrowLoadPackedKey<Schedule>(path)`. The archive contains the four
public encrypted indicators and raw substitution keys, with nine profile
words; derived FFT caches and secrets are omitted. Loading checks the exact
fixed length and public profile before allocating the packed key. Saving
stages a mode0600 file in a private directory, then publishes it with a hard
link that refuses to replace an existing destination. Failed publication
cleans only its own staging directory. This is closed-file publication,
not a power-loss durability guarantee or a payload authenticity check.

New packed checkpoint benchmark builds save `packed.bin` before preparing
substitution FFTs and reconstructing selectors. This preserves the raw key
for later reconstruction experiments even if reconstruction fails. It does
not save the associated secret keys before reconstruction; an interrupted
setup is not yet a complete refresh replay checkpoint. Existing completed
checkpoints without `packed.bin` cannot recover a discarded raw key.
Already-running pinned binaries retain their original behavior.

To reconstruct again using compatible retained data:

```cpp
auto packed = TFHEpp::CKKSNarrowLoadPackedKey<Narrow>(checkpoint / "packed.bin");
TFHEpp::CKKSNarrowPreparePackedSubstitutions<Narrow>(*packed, 2);
TFHEpp::CKKSNarrowReconstructSelectorsToFile<Narrow>(key.narrow, *packed, new_selector_directory);
```

Keep the checkpoint's original refresh key and input associated with this
archive; no new encryption keys are needed for reconstruction. Independently
validate the reconstructed selectors before treating the new cache as usable.

Retained N512 controls loaded and saved the existing public fixture with
byte-identical SHA-256, confirmed empty derived caches and mode0600 publication,
and rejected overwrite, an existing staging directory, wrong profile,
truncation and trailing data. Original fixture files remained unchanged;
no encryption keys were generated. Evidence:
`narrow-dd-packed-archive-controls-validation.json`. These controls validate
the archive helper, not a fresh/full packed benchmark setup.

The retained full tanh controller completed with exit0: refresh4019.056 s,
independent maximum refresh error2.0237842266785237e-8 (fixed23bit pass),
square error5.39755485425234e-8 below4.5615794863043815e-7. Secrets, slots,
input and corrected selectors were reused unchanged. This remains the
corrected-selector diagnostic fixture; its original error and structured-key
security are unqualified. The full N32768 basis controller subsequently passed
its prerequisite/headroom gates and began its initial checkpoint generation
with its previously pinned executable. Its full encrypted accuracy is pending.

The updated benchmark compiled and linked successfully (exit0), and its
`--storage` diagnostic produced a valid JSON report without encryption-key
generation. Its fresh packed-checkpoint retention branch has not yet been
executed end to end; existing live runs continue using older pinned binaries.


Packed reconstruction replay command (2026-10-10):
`ckks_narrow_packed_archive` supports the N512 and N16384 primary4/448-bit
packed profiles (selectorQ379, scale49, mask scale34, h31). In addition to
`--validate` and `--copy`, it accepts:

```sh
OMP_NUM_THREADS=4 TFHEpp/build/test/ckks/ckks_narrow_packed_archive --reconstruct /absolute/packed.bin /absolute/new-selector-directory 2
```

The optional last argument selects 0 through4 cached substitution depths;
the default is2. The destination must be new. Reconstruction consumes only
the public packed archive, generates no encryption keys, preserves completed
selector files, and does not require loading a full refresh key. The new
`CKKSNarrowReconstructSelectorFile<Schedule>(packed, directory)` factory
returns the selector-file source directly; the existing attachment function
still installs it into a compatible refresh key only after success.

Two actual N512 replays from the original public fixture at cached depths2
and3 completed with exit0. An independent checker summed all2048 original
reference terms with arbitrary-precision integers and one final reduction
mod2^379. All595 reconstructed selectors /609280 coefficients passed per
replay, including unused zero positions. All35 completed selector files had
identical hashes across cache depths. Existing destination and invalid cache
depth were rejected; retained fixture and completed selector files stayed
unchanged. Evidence: `narrow-dd-packed-archive-reconstruction-validation.json`.

This verifies the complete small-ring archive-to-selector-file replay path.
The N16384 CLI path is compiled but has not been executed with a retained full
packed archive. The original running full packed binary predates retention;
its eventual refresh checkpoint alone does not retain that raw archive.

Both reconstruction replay targets were registered by CMake and built
successfully. A further replay using the CMake-built archive command passed
the CMake-built independent checker and produced identical selector-file
hashes to the diagnostic binaries. A copied selector with one ciphertext
payload bit flipped was rejected by that checker; all originals stayed
unchanged. The recorded timing from this additional run is under concurrent
load and should not be compared as an isolated performance result.


Full-ring primary8 noise control (2026-10-10):
`ckks_narrow_auto_primary_projection STATE N PUBLIC_OUTPUT` reads only the
needed retained archive fields and original Q356 conjugation key. Compact
S2C payloads and irrelevant fixed-size key arrays are skipped after bounded
profile checks, avoiding a full refresh-key allocation. It selects45 existing
rows from the original89 to represent primary8 instead of primary4, with
secondary16 FFTs. Deterministic diagnostic masks are given exact zero phase
under the original retained dense secret. No secret/evaluation keys are
created or serialized. The output contains three input/old/projected public
ciphertext triples; reported phases use scale49. This tests key-switch noise,
not fresh ciphertext noise, auxiliary-scale lifting, or a full packed tree.

N512 passed: maximum projected slot error3.551562349152503e-9. N16384 failed
the fixed23bit bound2^-23=1.1920928955078125e-7: projected maxima across three
samples were1.250680448908549e-7,1.3304590910988487e-7 and1.2863148670242208e-7.
Original primary4 maximum4.35990950728044e-8 passed. Independent bounded15bit
integer-NTT phase decoding confirmed all input phases exactly zero and matched
all reported output metrics. Both retained state hashes stayed unchanged.
Evidence: `narrow-dd-auto-primary8-projection-validation.json`.

Consequently primary8 is not qualified for general use at Q356/scale49, and
production decomposition settings remain unchanged. Prior N512 packed
projection results do not justify adoption. This component result does not
prove the outcome of the untested Q427 packed extraction tree or an actual
auxiliary-lifted refresh schedule. Those scopes need their own accuracy gates.
The independent verifier is
`Parameter-Selection/python/estimates/narrow_ckks_auto_primary_projection_validate.py`;
exit2 indicates that the independently observed error exceeds the fixed bound.


Primary8 with auxiliary scaling (2026-10-10):
The direct Q356/scale49 failure above does not describe the existing factor
conjugation path, which lifts Q348/scale48 by8 bits, switches at Q356, then
rounds back to Q348/scale48. New controls use that exact auxiliary operation
on three nonzero complex encrypted inputs with minimum Gaussian sigma3.2
per ring. They reuse the original retained dense and conjugation keys.
Both N512 and N16384 passed independent integer phase decoding against the
known conjugated slot targets. The N16384 explicit-row projection maximum
error was9.88109032130982e-10, below2^-23. No encryption keys were generated,
and the complete retained state hashes remained unchanged. Evidence:
`narrow-dd-auto-primary8-auxiliary-validation.json`.

`CKKSPreparedAutoKeyFFT` and `makeCKKSPreparedAutoKeyFFT` now accept an
optional fourth template argument `PrimaryBits`, defaulting to the parameter's
original Bbar width. An explicit multiple of the raw width selects existing
matching gadget rows directly into the FFT cache, without allocating a second
raw key or changing ciphertext parameter types. The FFT mantissa and gadget
coverage guards still apply. For example:

```cpp
auto projected = TFHEpp::makeCKKSPreparedAutoKeyFFT<P, 356, 16, 8>(raw);
TFHEpp::CKKSEvalAutoWithAuxiliaryScale<P, 348, 48, 8>(output, input, 2 * P::n - 1, *projected);
```

Default prepared cache behavior remains the original primary decomposition;
streamed raw switching is unchanged. Choosing wider primary digits changes
switching noise, unlike repacking secondary arithmetic. The caller must check
its actual scale/lift schedule before adopting a projected cache.

Native primary8 cache controls matched explicit retained-row projection
ciphertexts exactly at N512 and N16384, all three samples per ring. Independent
full-ring maximum error1.0454354049423862e-9 passed the fixed23bit bound.
Evidence: `narrow-dd-auto-primary8-native-cache-validation.json`.
A default-primary4 packed extraction regression also matched all2048 archived
terms /2097152 ciphertext coefficients exactly, with retained files unchanged.
Evidence: `narrow-dd-primary-cache-default-packed-regression-validation.json`.

These qualify an optional component operation, not a complete refresh, the
Q427 packed extraction tree, security, or paper latency. The new option has
not been enabled in production refresh setup. Large in-progress evaluators
continue to execute their original pinned binaries.


Optional projected-conjugation refresh replay (2026-10-10):
`CKKSNarrowPrepareProjectedAutomorphisms<Schedule>(key)` prepares primary8 /
secondary16 caches for factor and output conjugation from existing raw rows.
It requires at least8 auxiliary bits and a raw primary width dividing8.
These derived caches are not serialized and are cleared when loading,
regenerating keys, or changing the cache policy. Existing explicit override
caches retain precedence. The factor-batch evaluator skips redundant staged
primary4 conjugation preparation when the projected cache is present. DD
relinearization and bounded tail rotations retain their existing algorithms.
Defaults are unchanged.

`ckks_narrow_auxiliary_replay --replay-primary8 STATE SELECTORS NEW_DIRECTORY`
enables only those projected conjugation caches; DD caches remain staged.
Complete retained N512 sigmoid refresh passed with independent maximum error
1.356399076362949e-9 below2^-23. Squaring the refreshed ciphertext with the
original retained square key also passed: independent maximum error
3.1450901746204564e-9 below4.50319595316329e-7. All36 original state/selector
files stayed unchanged, and no encryption keys were generated. Ciphertext
bytes differ from primary4 because primary decomposition changes switching
noise. Preparation0.00909 s and refresh8.82784 s were observed under concurrent
load; these are not an isolated speedup or a full-ring latency claim.
Evidence: `narrow-dd-auto-primary8-small-refresh-validation.json`.

This qualifies the complete small-ring sigmoid refresh and post-refresh
multiplication in the optional mode. Full-ring refresh with this mode remains
pending, as do the broader paper performance and structured-key security gates.
The current full basis and original packed setup evaluators remain pinned to
older binaries. Do not launch a competing large replay outside the full basis
controller's reservations while it is active; small bounded controls fit the
host reserve.

The optional primary8 conjugation cache now has a queued retained N16384 refresh and post-refresh square validation (`Parameter-Selection/python/estimates/narrow_ckks_primary8_full_run.py`). The controller freezes the executable qualified by the N512 test and its independent phase validators, waits for successful completion of the full frequency-basis validation, and checks a 20 GiB child memory budget with the existing packed job reservation. It reuses the original auxiliary8 state and diagnostically corrected selectors, compares retained secrets/slots/input, checks fixed 23-bit refresh and propagated product accuracy, and hashes retained files before and after evaluation. Evidence: `narrow-dd-auto-primary8-full-refresh-validation.json`; waiting is not full-ring accuracy qualification. The independent square CLI and bounded parent/larger child limits passed on the retained N512 result (`narrow-dd-primary8-full-runner-preflight-validation.json`). No security or paper timing qualification is implied.

Factor matrix accumulation now shares the selector spectrum entry and Shoup quotient between the two signed plaintext tables (`include/ckks/narrow_matrix.hpp`), preserving separate canonical modular sums and selector order. The two-prime integer reference checked 59,047,936 accumulated coefficients; four N16384 public matrix trials each matched 6,291,456 original coefficients. The native retained N512 default refresh passed independent 23-bit phase checks and retained the exact previous prefix SHA256, with all 36 original files unchanged and no key generation. Evidence: `narrow-dd-paired-selector-quotient-validation.json`. Public matrix timings are exploratory under concurrent load; full encrypted N16384 replay of this source and paper performance remain unverified. Running/queued large tests retain their earlier frozen executables.

Compact Double Decomposition keys now support optional wider-primary derivation through `makeCKKSCompactDDRelinKeyWithPrimaryBase<NewBg>` and the corresponding chain factory. NewBg must be a multiple of the retained primary width; each derived row copies the retained encrypted row with the same gadget exponent. This requires no secrets or fresh key generation and changes switching noise, so default schedules remain unchanged. The retained N16384 auxiliary8 chain projected from primary4 to primary8 passed independent exact key phases (2,621,440 coefficients), switch outputs (163,840 coefficients), and switched-phase convolution bounds (81,920 coefficients) at Q348/300/252/204/156. All 353,894,400 selected compact digits matched the retained rows; the source state hash is unchanged. The chain uses 160 primary rows instead of 315. Evidence: `narrow-dd-primary8-retained-full-dd-validation.json`; peak test RSS 1,562,300 KiB under a 2 GiB cap. Full refresh accuracy and performance with projected DD remain unverified. Run the retained component check with `ckks_narrow_dd_exact --aux8-checkpoint-primary8 STATE`.

The optional retained primary8 DD chain is integrated through `CKKSNarrowPrepareProjectedRelinearization<Schedule>(key)` and replay mode `--replay-dd-primary8`. It derives a primary8/secondary16 FFT cache from the serialized primary4/secondary4 chain, keeps the other stages on their existing cache policies, and skips redundant original DD preparation while the projected cache is present. Loading or key generation clears the derived cache. The retained N512 complete refresh and original-key square passed independent phase checks (maximum errors 1.357797968385218e-9 and 3.1445111756106072e-9), all 36 retained files were unchanged, and default replay retained the exact previous prefix. Evidence: `narrow-dd-primary8-dd-small-refresh-validation.json`. The frozen full N16384 run is queued behind successful completion of the existing full primary8 conjugation test (`narrow-dd-primary8-dd-full-refresh-validation.json`); it has not launched an evaluation child. The earlier full runner report guard was corrected to match the actual replay report fields and resumed with its frozen artifacts intact. A duplicate live-controller resume rejects with exit 2. Full DD refresh accuracy and speed remain unqualified.

Full-ring primary8/secondary16 DD product-tree qualification now passes for retained N16384 keys. The diagnostic encrypts new known test plaintexts with the original key and evaluates 31 encrypted leaves plus transparent beta across five levels for head-only, tail-only and all-slot complex targets. Public outputs contain only known targets and ciphertexts; `Parameter-Selection/python/estimates/narrow_ckks_dd_tree_validate.py` reads the original private secret prefix and independently decodes exact integer phases and canonical slots. All 24,576 slots passed the fixed 23-bit bound, with maximum error 5.728720693816637e-10. Copied wrong-level metadata and a flipped ciphertext bit48 reject. Evidence: `narrow-dd-primary8-retained-full-tree-validation.json` and `narrow-dd-primary8-full-tree-verifier-controls-validation.json`. The initial 2GiB attempt failed during worker-thread creation; releasing a temporary converted coefficient key before evaluation allowed the same four-thread cap to pass at 1,780,816KiB RSS. Failed evidence is retained separately. The source key file hash is unchanged. Run `ckks_narrow_dd_exact --aux8-checkpoint-primary8-tree STATE NEW_PUBLIC_OUTPUT`; this validates the tree, while whole refresh accuracy, matched performance, and security remain pending.

Section7 frequency-basis reconstruction is exposed as `CKKSNarrowEvalSeriesFromBasisRoots<Schedule>` for unweighted roots in the plan frequency order, encrypted under the same key as the basis public key. The all-function benchmark prepares the coefficient halves once and caches roots by frequency for the same encrypted input. Sigmoid evaluates eight roots, GeLU reuses all eight, and tanh evaluates fourteen while reusing four; total root evaluations decrease from 34 to 22. Caches are local to the input, use heap storage, and do not change checkpoint archives. Current native retained N512 lean replay passes all independent refresh/square/normalization checks and exactly matches all three archived ciphertext result files; all 36 checkpoint files remain unchanged and no keys are generated. Evidence: `narrow-dd-basis-root-reuse-validation.json`. Shared coefficient preparation is timed separately; cached GeLU reconstruction time is not standalone bootstrap latency. Full N32768 optimization validation and paper performance remain pending. Separately, the live full basis checkpoint passed read-only metadata and 36 selector-group length checks (`narrow-dd-current-full-basis-checkpoint-geometry-validation.json`); this does not validate encrypted phases or evaluation-key contents.

Full N32768 frequency-root reuse validation is queued through `Parameter-Selection/python/estimates/narrow_ckks_basis_root_reuse_run.py`. The controller freezes the N512-qualified current benchmark and independent basis decoder, waits for successful original full basis validation and termination of both optional primary8 refresh tests, and reuses the original full checkpoint without key generation. It compares all three refresh/square/normalization archives with the original validated replay, independently decodes phases, checks the 22-new/12-reused root counts, and verifies retained files against the original baseline before and after evaluation. The 40GiB child reservation includes a 10GiB host reserve and the original packed job remainder; its full cap remains reserved when the evaluator disappears before packed follow-up validation completes. Evidence: `narrow-dd-basis-root-reuse-full-validation.json`. The waiting controller is live, with no evaluation child launched; full encrypted root reuse and matched paper performance remain unverified. A superseded queue preflight is preserved because both sibling test terminations must be checked explicitly even if one fails before launching.

The retained full-ring DD switch now has a paired diagnostic command: `ckks_narrow_dd_exact --aux8-checkpoint-dd-timing STATE`. At N16384/Q348, both the current primary4/secondary16 and projected primary8/secondary16 paths match independent integer-NTT coefficient references on a uniform input including native-width garbage, signed/carry boundaries and negacyclic wrap. After warming both paths, eight alternated trial groups each run eight switches and verify the resulting ciphertext. The observed median ratio is 1.5354403541167592 under concurrent packed/basis load; it is not an isolated or complete refresh speedup. Cache sizes decrease from 501,743,616 to 253,755,392 bytes. The original retained state hash is unchanged and no keys are generated. Evidence: `narrow-dd-primary8-retained-full-dd-timing-validation.json`. This supports the optional cache optimization; default adoption still awaits full refresh validation.

Explicit selector NTT digit widths now accept up to60 bits when the worst-case two-prime CRT bound permits them. The default remains48/32/16 pending full encrypted qualification. N16384/Q385/plain37 and Q379/plain34 support60; N32768/Q518/plain39 supports56 and rejects60. `ckks_narrow_selector_wide_exact` checks262,144 full-ring coefficients against independent multiprecision negacyclic convolution, including dense signed values near the CRT bound, sparse wrap, native-width garbage and partial final digits. Retained N512 head/tail factors using60 match the baseline16, alternative chunks and single-column ciphertexts exactly and pass23-bit decryption gates; all36 original files are unchanged and no keys are generated. Run `ckks_narrow_factor_checkpoint --retained-small-aux8-wide60 STATE SELECTORS`; `--retained-aux8-wide60` provides the corresponding full-ring control, which has not yet run. Retained N16384 group preparation trials with529 selectors observed median48-bit3.36658s versus60-bit2.42485s (1.38836629x) under concurrent load, with the original group hash unchanged; this measures forward NTT preparation only. Run `ckks_narrow_selector_wide_exact --retained-aux8-timing SELECTORS`. Evidence: `narrow-dd-selector-wide60-validation.json`. Complete refresh performance and full encrypted wider-width factors remain unverified.

Standalone prepared conjugation keys are now accepted by the internal factor/reference kernels, avoiding allocation of unrelated raw rotation and product-tree keys. The existing full-key dispatch uses the same extracted conjugation helper. `ckks_narrow_factor_checkpoint --retained-aux8-wide60-bounded STATE SELECTORS` reads only the original dense/sparse prefix and final factor conjugation key, using checked retained square/DD/cache-mode suffix anchors to locate that payload. It skips unused payloads; this is a known-fixture loader, not validation of the complete key archive. Selector chunks8/5 and three columns keep the comparison under2GiB. Retained N16384 head30 and tail4 factors using60 matched baseline16, default48 alternative chunks and single columns exactly (196,608 unique ciphertext coefficients); all49,152 decoded slots passed23bits with maximum error1.2522476022539664e-9. Peak RSS 1818840KiB; the original private state and both selected groups retain their prior hashes and no keys were generated. Small-ring bounded and legacy-container checks passed; wrong square metadata and truncated suffixes reject. A rebuilt current-source default N512 refresh independently passes23bits, preserves the earlier fb28516d prefix hash exactly, and leaves all36 retained files unchanged. Evidence: `narrow-dd-selector-wide60-full-factor-validation.json`. Wider-selector complete refresh and N32768 encrypted56-bit factors remain unqualified; defaults are unchanged.

Frequency-basis batching now tiles across input boundaries, sharing selector setup between real/imaginary halves and independent ciphertexts. It retains only the current root tile plus one unfinished basis, checks input/frequency count overflow, and preserves finite coefficient validation for empty batches. The benchmark requests up to twice the per-half basis size (maximum18 columns), and its all-function root cache collects missing roots from both halves in one batch. Per-function reports include factor-group counts and the actual tile width. Retained N512 cross-input controls with three distinct coefficient ciphertexts and tile sizes5/7/12 exactly match per-input evaluation, including partial bases split across tiles; selector-group counts are105/70/35. Current all-function replay preserves all three refresh/square/normalization archives exactly and independently passes all phase gates; selector groups decrease140 to70 while roots remain22 new/12 reused. Generic single-function complex refresh also preserves its archive and passes independent phases. All36 checkpoint files retain their previous hashes and no keys are generated. Evidence: `narrow-dd-basis-joint-input-tiles-validation.json`. Compiled N32768 public payloads plus conservative scratch/root/stack allowances total34.28766GiB against the existing40GiB evaluation cap; this is a planning bound, not measured full-ring RSS. The previous waiting root-reuse controller was superseded before launching a child. Its replacement freezes this qualified binary, uses the same original checkpoint and dependency/memory gates, and additionally requires full SIG/GeLU/Tanh factor-group counts36/0/36. Evidence: `narrow-dd-basis-joint-input-full-validation.json`. Full encrypted joint tiling and paper latency remain unqualified; the active original large runs retain their frozen binaries.

Exact active-level CKKS torus multiplication now transforms the small signed polynomial once and chooses up to56-bit unsigned limbs using its actual L1 norm. The strict centered-prime bound is checked before transforming; zero polynomials return exact zero and oversized signed coefficients/bounds reject. No cryptographic modulus or noise distribution changes. Independent multiprecision checks cover245,952 coefficients at N64/N16384/N32768, Q128/385/448/518, dense near-bound signed digits, sparse wrap, native-width garbage, partial final limbs and zero inputs; all match the frozen16-bit reference. Sixteen simultaneous retained corrected-mask products also match262,144 coefficients. Paired serial calls observed median old16 0.041038687500000004s versus adaptive 0.01134028125s (3.6188421252779777x) under concurrent load. This measures the product kernel only. The shared full-ring encryption/scheme-switching smoke target passes. Evidence: `narrow-dd-active-product-adaptive-validation.json`. The retained narrow encryption keys are reused; the independent existing smoke target creates its own small switch-key setup. Full refreshed-input performance and selector error origin remain unqualified.

A later read-only selector identity audit does not reproduce the historical original head19 checksum. Both the original-labelled file in `run/selectors` and the diagnostic corrected file currently hash to69abb7c8..., matching the recorded corrected baseline. Buffered and direct reads agree on the corrected bit151 value; reversing that bit in memory recreates the historical original whole-group SHA256 a60720b0... and record SHA25616250cc8.... These controls do not modify either file. This establishes a current identity discrepancy confined to the known correction bit, without attributing prior writes, cache behavior or hardware faults. Earlier claims that the historical original remains intact cannot qualify its current contents. Evidence: `narrow-dd-selector-original-identity-discrepancy.json`. Corrected-fixture tests remain explicitly diagnostic; full basis/packed work uses separate checkpoints.

The initial fresh N32768 basis sigmoid completed in4717.9384s, excluding2937.0785s setup. Independent integer phase decoding passes refresh, square and normalization gates, with respective maximum errors5.2630583e-9,9.4724187e-9 and2.6925382e-9. The result archive hash is unchanged; this validation reuses the original checkpoint and uses no diagnostic selector correction. Evidence: `narrow-dd-basis-full-initial-sigmoid-independent-validation.json`. GeLU/Tanh, complete replay and the latest joint tiling remain pending. Accuracy passes do not establish practical paper performance: Table7 reports high-precision sigmoid2.79s on48 CPU cores with different parameters/hardware, versus this experimental four-thread run taking78.6min.

`ckks_narrow_basis_batch_checkpoint --phase-sample CHECKPOINT_DIRECTORY` reads only the existing dense/sparse prefix and samples every active indicator record plus one zero record per group. N512 and N32768 controls cover all35/36 groups, use a fixed16-bit integer-NTT phase reference and require the adaptive production product to match exactly. Every sampled coefficient stays within the39-unit noise bound; a body-bit100 change in a separate private copy rejects without changing the original group hash. No keys are generated. Evidence: `narrow-dd-basis-selector-phase-samples-validation.json`. These are sampled selector checks, not validation of every ciphertext or the full evaluation-key suffix, and do not qualify security.

The paired selector-matrix loop now uses an exact AVX2 Shoup kernel for four or more columns when the auxiliary prime is below2^62. Unsigned32-bit partial products recover the full64-bit high product; both lazy products and canonical sums stay below2^63, making signed lane comparisons exact. Other widths/moduli and builds without AVX2 use the scalar path. Independent128-bit references pass196,083,712 modular accumulations in each native/fallback build, including near-bound moduli, boundary selectors, partial vectors and scalar tails. Eight full-degree paired matrix trials retain identical6,291,456-coefficient sums after eight repetitions each; the observed median component ratio is1.539598 under concurrent load, not a complete bootstrap speedup. Retained corrected N16384 head30/tail4 controls now use five columns to exercise vector lanes and scalar tails:327,680 ciphertext coefficients match baseline16, alternative chunks and selected single-column outputs, and all81,920 slots pass23bits. N512 joint-input controls with tiles5/7/12 also pass exact per-input comparisons. All retained fixture hashes are unchanged and no keys are generated. Evidence: `narrow-dd-selector-matrix-avx2-validation.json`. Existing large jobs and the waiting full joint-input executable remain frozen without this kernel; whole full-ring refresh performance remains unqualified.

`ckks_narrow_basis_replay --[small-]lean-basis-fbs-replay all CHECKPOINT_DIRECTORY REPORT.json` provides a dedicated CPU retained-checkpoint runner. It shares the basis benchmark implementation in `test/ckks/ckks_narrow_basis_benchmark.hpp` and instantiates only the N512/N32768 lean profiles. The native optimized build completes in27.21s with roughly1.5GiB peak compiler RSS, versus the much larger all-profile benchmark build. It checks the public checkpoint profile before allocating evaluation keys, rejects existing report/result files, and offers no key-generation command. Function results now report six factor-stage wall-time totals, root-generation and basis-reconstruction times when sharing roots; completed groups flush progress lines. This profiling changes no ciphertext arithmetic.

Current N512 all-function replay with AVX2 and profiling matches all three earlier refresh/square/normalization archives byte for byte, independently passes all phase gates, leaves all36 checkpoint files unchanged and emits70 completed-group progress lines. Wrong-profile, invalid-function and existing-output controls reject under256MiB. Evidence: `narrow-dd-basis-profiled-replay-validation.json`. These are small-ring controls; cached GeLU reconstruction is not standalone FBS latency, and the measurements do not qualify full N32768 performance. The dedicated runner is compiled and tested. Separate legacy all-profile syntax checks exhausted6GiB and9GiB address-space caps; the monolithic benchmark rebuild remains pending. Existing frozen large executables and the queued full joint-input test are unchanged.

The full replay controller accepts `--exe` for the dedicated qualified runner, verifies its qualification hash, and freezes the executable, qualification JSON, controller and independent decoder before waiting. Profiled runs additionally require finite nonnegative factor-stage totals contained by root generation, and root generation/reconstruction contained by refresh time. Preflight controls reject an unrelated executable and a failed original full validation without launching an evaluator. The preceding waiting joint-input controller was verified live with no child, then superseded before evaluation. Its replacement uses the same retained full basis checkpoint, successful original-validation prerequisite, optional-test termination gates and40GiB child reservation. It will compare all three result archives and independently decode refresh/square/normalization phases, with full factor counts36/0/36. Evidence: `narrow-dd-basis-profiled-avx2-full-validation.json`; preflight controls are retained in the private `basis-profiled-controller-controls-20261010` directory. The controller is waiting with no evaluation child, so latest full AVX2/profiled performance remains unqualified. The original packed and basis jobs retain their frozen binaries and keys.

The compiled current N32768 storage profile gives167.37890625GiB raw selectors and409.1484375GiB for their complete default48-bit NTT cache, versus125.58266GiB physical host memory. Paper Section8.2 reports a dual-socket EPYC7K62 server with1TB DDR4-3200; Section5 represents sparse-operator values by stored selector NTT coefficients. The torus port instead recovers exact convolutions from small limbs. Its current file-backed path loads raw groups and prepares NTT chunks on each factor pass, so retaining compatible keys removes key generation without removing repeated transform work. This is a concrete remaining storage/evaluation optimization target. These payload counts do not constitute a measured timing breakdown or attribution of the whole performance gap. Evidence: `narrow-dd-basis-selector-memory-observation.json`.

The pre-cache-extension monolithic native refactor build completed successfully under the reserved14GiB address-space cap. Its finiteness diagnostic passes, and retained N512 all-function replay matches all three independently validated dedicated-runner result archives exactly. Evidence: `narrow-dd-basis-legacy-native-refactor-build-validation.json`. This closes the earlier6/9GiB build-cap limitation for that source snapshot; subsequent optional prepared-cache additions are separately compiled in the dedicated runner.

`CKKSNarrowNTTSelectorFile<Schedule>` adds optional reusable disk storage for exact prepared selectors when BLAKE3 is enabled. `prepare(NEW_DIRECTORY, raw_source, chunk_size)` transforms existing raw ciphertexts in bounded chunks without generating keys, publishes fixed-size portable group files through pending files, and preserves the completed cache. Reopening requires all group profiles/lengths to match; each loaded chunk checks canonical residues, its unkeyed BLAKE3 payload digest and a fingerprint of the corresponding current raw selector. These checks detect accidental corruption and stale-source mismatch, not adversarial authenticity. The original raw selectors remain available for reference mode. The factor loop streams prepared chunks directly and assigns cache I/O/verification to selector-load timing rather than selector-NTT time. Derived caches remain outside key serialization.

Create the tested lean-profile cache once with `ckks_narrow_selector_ntt_file CHECKPOINT_DIRECTORY NEW_CACHE_DIRECTORY`. Reuse it with `ckks_narrow_basis_replay --small-lean-basis-fbs-replay all CHECKPOINT_DIRECTORY REPORT.json PREPARED_SELECTOR_DIRECTORY` (omit `small-` for full degree). Retained N512 creation/reopen with chunks5/7 checks40,212,480 exact spectrum words. Raw-path, cached-path and previous independently decoded all-function refresh/square/normalization archives are byte-identical; all36 original checkpoint files retain their hashes and cached evaluation reports zero selector NTT work. Copied payload damage, incompatible digit metadata, truncation, changed raw selectors and existing output-cache directories reject. Group files are private0600. Evidence: `narrow-dd-prepared-selector-ntt-cache-validation.json`. Small timing observations do not establish a speedup: added cache I/O and integrity checks can offset saved transforms. Full N32768 cached evaluation and current-source monolithic rebuild remain unqualified. Existing large jobs and the queued AVX2/profiled runner retain their frozen binaries and do not use this optional cache.

`CKKSNarrowNTTSelectorFile::prepareChunk` and `loadChunkArchive` support bounded partial prepared-selector archives with a distinct magic and explicit expected range. They share the complete cache's record codecs, raw-source fingerprint and payload checks. A partial archive cannot be attached as a complete evaluation source; invalid reads preserve the caller's output. Raw fingerprints now hash complete contiguous coefficient arrays on little-endian hosts, retaining the prior digest format while avoiding one hasher call per coefficient. Run `ckks_narrow_selector_ntt_file --sample-chunks CHECKPOINT_DIRECTORY NEW_SAMPLE_DIRECTORY` for the retained range control. Current N32768 beginning/middle/end ranges in head0/head30/tail0/tail4 cover22 raw ciphertext records and31,719,424 exact spectrum words. Thirteen wrong-range or partial-attachment controls reject. The public checkpoint profile and selected raw headers/ranges retain their hashes, and no keys are generated or full evaluation-key suffixes loaded. Current N512 replay reopens the earlier complete cache and preserves all three previously independently decoded result archives exactly with zero selector NTT work. Evidence: `narrow-dd-prepared-selector-ntt-full-chunks-validation.json`. This qualifies full-degree chunk transport only, not every selector, complete N32768 cached refresh, paper latency or security.

The factor kernel releases permutation and base-table spectrum buffers immediately after the selector matrix stage, before inverse CRT and conjugation allocate workspaces. At N32768 the permutation table alone occupies 264.5MiB. Current retained N512 all-function replay under a 2GiB address-space cap produces the same three ciphertext archives as the previously validated cached runner, without generating keys. Evidence: `narrow-dd-factor-buffer-lifetime-validation.json`. The full-degree paired-factor retry reuses the head cache left by the earlier allocation failure; complete qualification remains pending.

The current matrix path computes Lemma 5.1 permutation indices from one odd exponent per selector instead of retaining a degree-sized index array. Natural-order indices use the low bits of the exponent product; at N32768 this reduces index storage from277,348,352 bytes to8,464 bytes. `ckks_narrow_permutation_exact` compares against independent coefficient-domain substitution followed by fresh NTTs across N512/8192/16384/32768 and23 primes (15,968,256 points). Current retained N512 all-function ciphertext archives remain byte-identical. Evidence: `narrow-dd-implicit-ntt-permutation-validation.json`. Latest-source full refresh and performance remain unqualified. The earlier frozen full-factor retry has completed one fresh and one cached head trial with exact ciphertext identity, but cache loading/verification on shared storage took202.88s and made that cached trial slower (229.41s vs43.93s fresh). Disk caching remains optional; it has no demonstrated full-refresh speedup.

The matrix kernel constructs a `SelectorMultiplierFactory` once per prime and uses its exact invariant-divisor reciprocal to produce Shoup constants without division at each selector point. The reciprocal estimate underestimates the quotient by at most one; an exact remainder comparison corrects it. `ckks_narrow_selector_quotient_exact` checks3,000,180 quotients against independent128-bit division and234,014,040 modular products across all23 transform primes and small, power-of-two and near-boundary moduli. Six paired N16384 matrix trials preserve2,097,152 output sums each, and retained N512 complete all-function replay archives remain identical. Evidence: `narrow-dd-selector-reciprocal-validation.json`. Component timings are not isolated; the compile-time-divisor control was slower with the factory and runtime-divisor observations are noisy. No complete-factor or refresh speedup is established. Current-source full-degree encrypted refresh remains pending.

The full prepared/fresh factor retry now completes successfully under2GiB using its frozen early-buffer-release binary: eight paired trials preserve1,048,576 unique ciphertext coefficients,262,144 head/tail slots pass23bit accuracy, and selected complete raw groups plus the checkpoint secret/public prefix retain baseline hashes. The initial allocation failure remains in the evidence. Cold cache trials take189-229s on shared storage, so no complete-refresh speedup is established. Evidence: `narrow-dd-prepared-selector-ntt-full-factors-validation.json`. This snapshot predates implicit indices, reciprocal generation and point blocking.

Current accumulation traverses256-point destination blocks while retaining selector order at every point; it does not reconstruct or round at block boundaries. Twenty-four paired N512/N32768 controls compare75,003,904 unique modular sums, with122,880 independent small partial-block oracle checks. Raw and cached retained N512 all-function refresh archives remain identical. Four-selector/eight-column and forty-selector/eighteen-column component observations favor blocking, but concurrent jobs make the timings non-isolated and no whole-refresh speedup is qualified. Evidence: `narrow-dd-matrix-blocks-validation.json`.

The initial full GeLU result now passes independent refresh/square/normalization phase checks, with maximum errors1.55412e-8,1.26216e-7 and2.15174e-9. Its result hash remains unchanged. The older frozen evaluator took5945.51s for GeLU and has proceeded to Tanh; this remains far from the paper CPU latency2.71s. Evidence: `narrow-dd-basis-full-initial-gelu-independent-validation.json`. All-function replay, latest-source full refresh, paper performance and security remain incomplete.

`ckks_narrow_basis_factor_cache_checkpoint --full-fresh CHECKPOINT_DIRECTORY NEW_DIRECTORY` validates current full-degree factors with retained keys while avoiding prepared-cache creation and shared-storage cache reads. Current N32768 head30/tail4 with eight columns and chunk4 completes exit0 under2GiB in2:36.47. Two repeated trials per group compare1,048,576 ciphertext coefficients;262,144 decoded slots pass23bit accuracy, with the same maximum errors as the older full factor control. Complete selected raw groups and the checkpoint secret/public prefix retain baseline hashes. No keys or prepared-cache payloads are generated. Peak RSS1788912KiB includes key preparation. Evidence: `narrow-dd-current-full-factor-blocks-validation.json`. The prior control did not retain ciphertexts for cross-version coefficient comparison, so only repeated current-trial equality and cross-version decoded maximum-error equality are claimed. Complete current-source refresh, matched latency and security remain unqualified.

`CKKSNarrowBlindRotateBasisPlansBatch` evaluates the distinct frequency union of several series for each coefficient input, returning separate root maps per input. Duplicate plans do not regenerate roots; empty input/plan collections perform no selector evaluation and a zero tile rejects. The caller sets the column tile. The all-function runner uses22 columns for the eleven frequencies across real/imaginary halves, reducing selector groups from70 to35 at N512 (and planning72 to36 at N32768). Native and CMake raw/cached N512 all-function archives remain exact, and three distinct-input controls compare33,792 root coefficients across22/33-column tiles. Evidence: `narrow-dd-basis-frequency-union-validation.json`. The 22-column N32768 conservative payload estimate is35.9GiB under the40GiB evaluator cap; this is not a measured full-run RSS. Full encrypted union evaluation and matched latency remain pending. Inspect the profile without generating keys using `ckks_narrow_basis_replay --lean-basis-storage NEW_REPORT.json`.

The22-column full-factor control completed its first head factor but exceeded the2GiB allocation limit on the repeated trial while retaining the first ciphertext vector. Its revised control stores private0600 portable reference ciphertexts, reads one reference at a time for coefficient comparison/decryption, and retains the same mathematical work and memory cap. The N51222-column reference control passes. The streamed-reference full retry also exhausted the2GiB address-space cap on its second head trial, so a new retry uses the same22-column calculation with a3GiB diagnostic cap. Both failed attempts are retained; no production memory qualification is inferred from these component runs. The initial failure remains preserved in `narrow-dd-basis-frequency-union-full-factors-validation.json`. These reference files contain test factor ciphertexts, not prepared selector caches or regenerated keys.

The22-column N32768 factor control now completes exit0 under a3GiB address-space cap, with peak RSS1,943,208KiB. Two head30 and two tail4 trials compare2,883,584 ciphertext coefficients exactly;720,896 decoded slots pass23bit accuracy (head max8.26952e-10, tail max7.82567e-10). Selected complete raw groups and the checkpoint prefix retain baseline hashes. Both failed2GiB attempts remain preserved in `narrow-dd-basis-frequency-union-full-factors-validation.json`; this does not qualify a2GiB cap. The current CMake N512 archives also pass independent phases for all three functions.

The idle older profiled replay queue was superseded after verifying its process identity, waiting state and absence of children. A new frozen combined-frequency-union full replay controller is live and waiting, with six artifact hashes verified. Evidence: `narrow-dd-basis-frequency-union-full-validation.json`. It requires completed original validation, preceding tests to terminate, the40GiB child allowance, packed-process reservation and10GiB host reserve. Its ciphertext identity, independent phase, source hash, retained-file integrity and stage timing gates remain enforced, with22/0/0 new roots and36/0/0 group counts. All-function timing includes joint root generation in the first function and reconstruction alone thereafter; it is not a standalone paper latency measurement. No full union refresh result is available yet.

For scratch primes below2^61, the selector matrix retains sums below2p and combines the Shoup lazy product with one reduction by2p, then canonicalizes before returning each block. Products plus sums stay below4p<2^63; larger moduli retain canonical accumulation. Sixty-four native and64 scalar controls cover N512/N32768, nonzero initial sums, maximal canonical residues, SIMD tails and the modulus boundary/fallback, comparing161,542,144 unique modular sums per binary. Current CMake raw/cached N512 all-function ciphertext archives remain byte-identical, independent phases pass and36 retained files retain hashes. Evidence: `narrow-dd-lazy-matrix-validation.json`. Same-blocking component timing controls favor lazy accumulation for active primes, but concurrent jobs prevent an isolated speedup claim.

Current lazy-matrix full22 factor control is live under3GiB with retained keys, frozen executable and repeated trials. Its first head reference archive is byte-identical to the preceding canonical build across1,441,792 ciphertext coefficients; both head trials pass and360,448 head slots pass23bit accuracy. Tail qualification is pending in `narrow-dd-lazy-matrix-full-factors-validation.json`. The waiting full union replay keeps its previously qualified frozen executable; this subsequent kernel optimization is not included in that snapshot. Complete current-source refresh and paper latency remain unqualified.

The current lazy-matrix full22 factor control has completed exit0. Both head30 and tail4 portable archives match the preceding canonical build byte-for-byte across2,883,584 ciphertext coefficients. Repeated current trials also match, and all720,896 decoded slots pass23bit accuracy; selected complete raw groups and checkpoint prefix retain hashes. Peak RSS1943076KiB under3GiB. Evidence: `narrow-dd-lazy-matrix-full-factors-validation.json`. This qualifies encrypted full-degree factors, not complete current-source refresh or matched paper performance. The queued full replay still uses the earlier frozen snapshot.

The dedicated runner now exposes both create-once and replay modes without rebuilding the large legacy benchmark. Current CMake N512 checkpoint creation/replay passes all three independent refresh/square/normalization phases, exact ciphertext identity, one creation and zero replay key generations, and36 unchanged retained files. Seven N32768-mode preflight controls under256MiB reject path/function conflicts before full key allocation. The fresh fixture is N512 only. Evidence: `narrow-dd-basis-runner-lifecycle-validation.json`. Original full evaluators and the queued frozen full union replay keep their existing keys and executables; current-source full refresh and matched latency remain unqualified.

Raw selector files now expose bounded reads through the public selector-source interface. Staged factor evaluation reads and transforms one raw chunk instead of loading a complete group. For N32768 and40-record chunks, raw vector capacity falls from4,992,270,336 to188,743,680 bytes (4,803,526,656 bytes saved); this is a static buffer bound, not a measured complete-refresh RSS reduction. Sources without raw-chunk capability keep the original whole-group path, and prepared-cache sources retain their existing path. Read and selector-NTT times are reported separately.

Current CMake N512 streamed/legacy input-owner controls pass, including1,260 bounded reads with full-group loading forbidden and malformed short-chunk rejection. Raw and cached all-function ciphertexts remain byte-identical and independently pass, with36 unchanged files. Current N3276822-column head30/tail4 control completes exit0 under3GiB:720,896 slots pass23bit accuracy and both portable factor archives match the preceding build across2,883,584 coefficients. Complete selected raw groups and checkpoint prefix retain hashes. Evidence: `narrow-dd-raw-chunk-streaming-validation.json` and `narrow-dd-raw-chunk-streaming-full-factors-validation.json`.

The former idle full union replay queue was superseded by the current qualified raw-streaming build; six frozen artifact hashes were verified and no active evaluator was interrupted. Evidence: `narrow-dd-raw-chunk-streaming-full-validation.json`. Its original-validation, preceding-test termination,40GiB child and10GiB reserve gates remain unchanged. Complete current-source refresh, peak RSS, paper latency and security remain unqualified.

The original frozen N32768 initial all-function evaluator completed exit0 with no selector correction. Independent integer-phase decoding passes the saved input, refresh, retained square and normalized square for Sigmoid, GeLU and Tanh. Refresh maximum errors are5.26306e-9,1.55412e-8 and2.07608e-8; all satisfy the fixed23bit gate. Tanh square max5.73711e-8 and normalized-square max2.37425e-9 also pass. Tanh refresh took6774.21s (112.9min), excluding setup; this older snapshot does not reproduce the paper CPU4.51s. Evidence: `narrow-dd-basis-full-initial-all-independent-validation.json`. The parent controller is hashing retained files before replaying with the same keys; retained replay and latest-source full refresh remain incomplete.

The original full controller has saved a37-file retained baseline and started all-function replay with the same checkpoint keys and encrypted input. Replay setup takes72.24s versus2937.08s for initial creation; it does not generate new keys. Replay accuracy/identity and latest-build end-to-end results remain pending.

An external30-bit AVX2 matrix prototype passes80 independent modular-oracle controls (48,664,576 unique sums),400,016 quotient comparisons and modulus/boundary checks. Its four-prime22-column matrix is slower than the current two-prime64-bit kernel. Padding to24 columns with two zero columns yields only an observed6.5percent matrix gain on the busy host, before additional NTT/CRT costs. It is not integrated and has no encrypted-factor or complete-refresh qualification. Evidence: `narrow-dd-rns32-matrix-prototype-validation.json`. Production sources and the queued frozen build remain unchanged.

For long validation runs, `../Parameter-Selection/python/estimates/narrow_ckks_live_progress.py` records read-only process CPU/RSS and selected file metadata. Supply `--job NAME:PID:START_TICKS` from controller evidence, a new `--output` JSONL history and new `--summary` JSON path; optional `--artifact PATH` tracks file updates. The observer guards PID identity, preserves existing outputs, writes a private atomic summary, and defaults to30-second samples for30minutes. It reads no key payloads and does not control jobs or establish accuracy/completion. Live/stale/absent-process, overwrite-preservation and actual sampling controls pass in `narrow-dd-live-progress-observer-validation.json`. The current observer watches packed reconstruction, original same-key replay and the idle latest-build queue while production remains frozen.

Deferred factor rescaling (2026-10-11): the N16384 retained functional
Sigmoid candidate with per-selector rescaling failed the unchanged 2^-23
accuracy target (independent error 6.607062634e-7, evaluation 418.633 s).
Its retained files and saved keys were unchanged. The sampled head30/tail4
factor control independently confirmed that rounding selectors before the
weighted sum amplifies error; rescaling completed factors instead gives
1.85e-8 to 2.54e-8 sampled error. `CKKSNarrowDeferredRescaleFactorSource`
now computes at retained379/49 selector and345/49 factor levels, then publicly
rescales completed factors to337/41. The target DD tree is unchanged.
The optional runtime factor source is excluded from key archives and reset
on loading or generating keys. All196608 ciphertext coefficients in six
sampled integrated factors match the independently checked control exactly
(`narrow-dd-ll14-deferred-factor-integration-validation.json`). This is a
component result, not full-refresh, security or paper-performance qualification.
The corrected sequential full replay uses the same saved keys and public input
(`narrow-dd-ll14-deferred-fbs-validation.json`); it stops on an accuracy failure.

The corrected full Sigmoid replay is terminal and still fails the fixed2^-23
accuracy gate: independent error1.877998030e-7, remaining-level square
error3.656644741e-7 versus3.533396362e-7 tolerance. Evaluation273.370 s.
All37 original files, both saved keys, sources and reused input are unchanged.
No new keys were generated. GeLU/Tanh were not launched after this failure.
The default retained N512 regression passes all3 exact archives and independent
phase/product checks (`narrow-dd-deferred-factor-default-small-validation.json`).
A frequency1 diagnostic now separates decoded-factor product error from
rotation/tree error (`narrow-dd-ll14-root-noise-isolation-validation.json`).

The corrected independent frequency1 root diagnostic finds3.230619952e-7
error in the product of decoded factors before rotation/tree evaluation,
versus9.284311166e-8 additional evaluated-root/factor-product difference.
Its initial decoder used an incorrect tail width; that receipt is archived,
and only independent decoding was repeated after correction.
A separate sampled plain37-versus34 table control uses the same retained
379/49 selectors and345 conjugation key, with public auxiliary3 at factor342/49.
After factor rescaling to334/41, sampled errors are1.19e-8 to1.56e-8,
versus1.87e-8 to2.54e-8 for the original plain34 factors. Independent phases
and integrity checks pass (`narrow-dd-ll14-factor-plain37-isolation-validation.json`).
No keys were generated. This higher-precision profile is not yet integrated
or qualified by a full functional refresh.

Public reconstruction-mask normalization (experimental,2026-10-11):
Increasing table precision to37 and root scale to42 alone leaves the complete
frequency1 error at3.296103738e-7. The original plain34 reconstruction masks
have a public gain which repeats over31 factors; its product deviates from1
by3.226885374e-7. Independent decoding of saved outputs predicts7.658078511e-8
root error after dividing by this public gain. The new opt-in
`CKKSNarrowPublicMaskNormalizedFactorSource<Schedule,SelectorSchedule>` wraps
an existing factor source and corrects the existing transparent beta leaf.
It uses only public masks recreated with the original selector reconstruction
parameters, incurs no extra multiplication level, and generates no keys.
It applies only to selectors produced by `CKKSNarrowExtractSelectorTerms`
with those masks, not directly encrypted selectors. The default path remains
unmodified unless this wrapper is explicitly attached.
The deferred factor adapter now permits a separate selector archive schedule
so original selector bytes can be reused with different table precision.
The encrypted normalized frequency1 root and default regression are pending;
plaintext prediction is not full encrypted/functional qualification.

The integrated encrypted frequency1 normalized root passes2^-23:
independent error7.670712721e-8 (`narrow-dd-ll14-normalized-root-validation.json`).
All37 consumed original files, saved DD chain and switch key are unchanged;
no keys generated. The default retained small all-function regression remains
byte-identical and passes independent phase/product/integrity checks
(`narrow-dd-public-mask-default-small-validation.json`). Full normalized
Sigmoid/GeLU/Tanh replay is now running at original fixed gates
(`narrow-dd-ll14-plain37-scale42-fbs-validation.json`), reusing saved keys
and the exact prior input. The frequency1 result does not prove full FBS,
security or paper-performance qualification.

The normalized plain37/root42 full replay is terminal: Sigmoid passes
(independent6.131495673e-8,23.959 bits; remaining-level square1.127581070e-7,
273.897 s), but GeLU still fails the fixed2^-23 accuracy target
(independent approximately1.84554e-7,277.180 s). Tanh was not launched.
All37 original files, supplied saved keys, pinned sources and reused input
are unchanged. No keys generated. The repository experimental runner uses
the normalized profile; it remains unqualified for all three functions.
The next component checks plain40/root43 with56-bit matrix scratch,
publicly projected DD keys333/290/247/204/161 and the same345 conjugation key
(`narrow-dd-ll14-plain40-scale43-root-validation.json`). It generates no keys
and reuses the same structured input.

The plain40/root43 normalized root does not materially improve the retained
packed-selector residual (independent8.096710579e-8). A direct expanded
selector control using the SAME retained dense/structured private keys and
conjugation345 key gives sampled factors approximately4e-10 before rescale
and2.1e-9 to2.6e-9 at root43 (`narrow-dd-ll14-direct-selector-pilot-validation.json`).
Two groups took18.0 and17.7 s to generate and are retained for reuse.
`CKKSNarrowEncryptSelectorGroup` now exposes bounded expanded-group setup
without regenerating automorphism/DD/private keys. The remaining34 groups
are being generated once; subsequent tests load them. Direct selectors do
not use public packed-mask normalization. Full GeLU/Tanh/Sigmoid evaluation
is queued after setup qualification, in that order to fail early on GeLU.
This expanded Section5 control does not qualify high-precision coefficient
packing. Full-function accuracy, security and paper-speed gates remain open.

2026-10-11: Direct expanded selector setup completed and saved: 34 new public groups plus 2 reused pilot groups; no new private/automorphism/DD keys. Full retained-key replay independently passes GeLU (6.43148e-8,23.890 bits,271.456 s) and Tanh (6.97954e-8,23.772 bits,710.183 s), including the remaining-level square. Sigmoid is running in the same original controller. This is accuracy progress, not paper performance or coefficient-packed high-precision qualification. Paper Table7 CPU latencies are GeLU2.71 s/Tanh4.51 s; hardware differs. Evidence: narrow-dd-ll14-direct-fbs-validation.json.

2026-10-11 update: Sigmoid also independently passes (3.09960e-8,24.943 bits,272.317 s) and the remaining-level square. All3 direct expanded-selector functions pass accuracy; final retained-file integrity check is still running in the original controller. Repository ckks_narrow_ll14_retained_fbs.cpp now uses this tested direct profile, with exact body identity after equivalent include normalization (narrow-dd-ll14-direct-repo-runner-validation.json). CLI additionally requires the saved direct-selector-directory. No setup/key generation on replay. Accuracy is comparable to paper Table7, while CPU runtime is approximately100x GeLU/98x Sigmoid/157x Tanh the reported times, on different hardware. Compact coefficient-packed high-precision path and paper speed/security qualification remain open.

2026-10-11: Direct full GeLU/Tanh/Sigmoid retained replay is terminal/all_checks_passed, including all consumed original files, all36 saved direct selector groups, supplied saved keys, sources and reused input integrity. Bounded coefficient extraction added with default full-range behavior; N512 reused fixture passes2120 exact term comparisons and2 invalid-range controls. Evidence: narrow-dd-ll14-direct-fbs-validation.json and narrow-dd-packed-range-extraction-validation.json. Next component is packed435/selector387,delta57 with primary8 and all14 depths cached; bounded samples precede any full reconstruction. No claim of compact high-precision or security/performance completion.

2026-10-11: High-precision coefficient-packed pilot qualified: public packed435 key (primary8), selector387/57, masks34, all14 substitution depths cached. Same retained dense/structured private keys;14 required new public substitution keys and4 indicators generated once and saved (1,420,296,265 bytes). Setup20.7795 s, cache preparation6.0592 s,16 sampled extraction terms15.6124 s. Independent full-ring phase comparison against encoded public masks max5.31991e-12 across262144 coefficients. This is bounded component evidence, not full group/factor/FBS qualification.

Bounded extraction API defaults to all coefficients and matches2120 saved N512 terms exactly. Public archive tool accepts the435/57/primary8 profile and explicit cache depths through nbit, reports spectrum bytes and extraction progress. N512 all9-depth default reconstruction preserves all35 group files exactly. Optional new scratch directory places the random-access accumulator on local disk; positive controls preserve all35 files and clean scratch, existing-scratch control rejects without changing unrelated contents. Evidence: narrow-dd-packed-range-extraction-validation.json, narrow-dd-packed-archive-high-precision-tool-validation.json, narrow-dd-packed-local-accumulator-validation.json.

Full reconstruction network-accumulator attempt measured512/65536 terms in81.3712 s (approximately2.9 h projection), exceeding the2 h cap; owned evaluator was intentionally stopped and its temporary accumulator removed. Saved packed archive was retained. Full reconstruction now uses local scratch /tmp/ckks-packed-ll14-local-20261011 and network-disk final selectors with all14 depths cached, session72994,2 h wall cap, no key generation. Evidence: narrow-dd-ll14-packed-precision57-reconstruction-validation.json. No full reconstruction or high-precision packed functional claim until actual completion and independent validation. Paper-speed/security gates remain open.

2026-10-11: Local accumulator alone did not remove the dominant reconstruction cost:768/65536 in116.277 s still projected about2.76 h. Owned second attempt intentionally stopped; owned scratch removed and saved public packed archive retained. Kernel profile then isolated Auto435 primary8 at57.54 ms/call and serial public mask multiplication at82.88 ms/call. Parallelizing independent torus digit rows reduces mask kernel to32.10 ms (2.582x) with all12 complete mask ciphertexts plus automorphism control byte-identical; original serial body is retained for n<8192. Full N16384 saved-key16-term extraction replay is byte-identical to prior serial pilot, independently max5.31991e-12, no new public or private keys. Default N512 all35 selector group files are byte-identical and scratch cleanup passes.

Evidence: narrow-dd-packed-mask-parallel-profile-validation.json, narrow-dd-ll14-packed-precision57-parallel-pilot-validation.json, narrow-dd-packed-mask-parallel-default-small-validation.json. Full reconstruction now uses the separately compiled qualified parallel archive tool, all14 cached depths, local accumulator /tmp/ckks-packed-ll14-parallel-20261011, network retained final selectors, session88587 and2 h cap. Evidence: narrow-dd-ll14-packed-precision57-reconstruction-validation.json. Do not alter pinned reconstruction sources while it runs. No paper online-speed or full high-precision packed FBS qualification claimed.

2026-10-11 measured parallel full-reconstruction progress: first256 terms23.612 s (approximately100.74 min linear extraction estimate), compared with38.9765 s using serial masks/local accumulator. This is a progress estimate, not a completed reconstruction benchmark. Subsequent root44/table42 component sources are prepared in ll14-packed-precision57-root44-20261011, using sourcefactor345/57, public post-factor rescale13 to332/44, public DD projections332/288/244/200/156, output112/44 and the original345 conjugation key. No new private, Galois, DD or input generation required for that component. Do not launch it before full selector reconstruction qualification.

2026-10-11: Opt-in integer NTT Double Decomposition is implemented in prepared_dd_ntt.hpp. The N16384 retained-key component uses primary64/secondary32 and a100-bit auxiliary switch at data333/43, projecting/rescaling the nonlinear term through key433. Generic CKKSMultRelinThenRescale independently passes all output slots (max2.29084e-9); all7 raw DD encryption rows and signed integer convolution/rounding bounds are checked. One new public DD433 was generated once under the retained dense key and saved; repository replay generates zero private/evaluation keys. Switch timing0.01436 s versus0.08064 s FFT (5.61x) is preliminary under concurrent reconstruction, and does not imply a full FBS speedup. Evidence: narrow-dd-wide-primary-ntt-component-validation.json. The original default backend remains in use pending full-tree/integration qualification.

The actual integer noise sampler currently truncates a continuous normal sigma3.2 toward zero. Independent raw DD phases measured integer RMS2.82691 over114688 coefficients, consistent with analytic truncation RMS2.82830. Security accounting must use this actual distribution; configured continuous sigma is not the integer sigma. No sampler or saved-key changes have been made. Evidence: narrow-dd-integer-noise-model-observation.json.

Full high-precision packed FBS runner ckks_narrow_ll14_packed_fbs.cpp compiled with actual exit0. It requires retained input, saved DD chain and switch key, and reconstructed selectors; it never encrypts input or generates keys. Profile: selector387/57, table42/scratch54, public mask-normalized beta, factor332/44, DD332/288/244/200/156, output101/44, remaining-level square57/44. Independent full-phase verifier checks all slots at fixed2^-23 accuracy gate. Session87335 waits for the existing root44 controller2058209 (FBS controller2060333); execution is conditional on original reconstruction and root qualification, stopping on failed gates. Authoritative receipt: narrow-dd-ll14-packed-precision57-fbs-validation.json. Compilation/queue is not full functional validation.

2026-10-11: prepared_dd_ntt_chain.hpp prepares all product-tree levels by public projection from one saved DD433/64/16 master, with secondary32 NTT caches and auxiliary100. Repository ckks_narrow_dd_ntt_tree.cpp compiled and completed existing generic five-level productTree integration at data333/290/247/204/161, output118/43. Cache367001600 bytes, preparation2.78857 s, NTT tree2.93351 s versus FFT4.56757 s (1.557x), measured under concurrent reconstruction. Independent all32-factor/product integer-phase verification is running in original session70820; no full-tree accuracy claim until terminal receipt. All private/evaluation keys reused, no generation. Evidence: narrow-dd-wide-primary-ntt-tree-validation.json. This is not full FBS integration; existing full-function runners still use the qualified FFT DD backend.

2026-10-11: Five-level NTT DD tree validation is terminal/all_checks_passed, including557056 independently computed phase coefficients (all32 inputs plus both full outputs), retained source/key/ciphertext hashes unchanged. NTT maximum error1.50069e-8; baseline FFT1.29374e-8, both below fixed2^-23. Timings remain preliminary,1.557x tree speedup, not full FBS improvement. Evidence: narrow-dd-wide-primary-ntt-tree-validation.json.

Security tool now accepts --noise-stddev and rejects nonfinite/nonpositive inputs. Historical default3.2 is preserved; each report explicitly records the discrete-Gaussian RMS proxy and omits exact sampler-distribution qualification. Four completed dense screens atN16384/q435 with actual integer RMS2.828301687834186 give weakest128.9583 classical generic bits. This closes incorrect variance accounting for that proxy only; actual structured/ring/correlated/circular keys and exact noise distribution remain unqualified. No overall128-bit security claim. Evidence: narrow-dd-actual-noise-dense-security-validation.json and Parameter-Selection/python/estimates/narrow_ckks_ll14_actual_noise_dense_security.json.

2026-10-11: Componentwise complex packed FBS runner ckks_narrow_ll14_packed_complex_fbs.cpp compiled, current source hash recorded with actual exit0. It checks f(real)+i*f(imaginary), complex amplitude-scaled square bounds and independent full phases. Both halves share selector transforms with tile2*degree and28GiB evaluator cap. Saved input, private/evaluation keys and reconstructed selectors reused; no input encryption or new keys. Original session79163 waits for original real packed FBS controller2060333's terminal qualification; failed gates stop all dependent work. Receipt: narrow-dd-ll14-packed-precision57-complex-fbs-validation.json. No full complex accuracy claim from compilation/queue.

Opt-in NTT DD integration is compiled and running in an isolated1.1MiB copy of include headers at ckks-diagnostics/narrow-ntt-integration-20261011. Production narrow_bootstrapping.hpp is unchanged, preserving live packed jobs' pinned sources. The proposed optional product-tree evaluator is derived nonserialized state, cleared on load, bypasses FFT DD cache setup when active, and receives raw keys without capturing the key object's address. Helper publicly projects one saved DD433/64/16 master into shared auxiliary100/secondary32 NTT caches for all five333/43 tree levels. Root integration compares unchanged input+36 factors byte-for-byte to already independently qualified FFT reference and independently decrypts only changed root plus original sparse input. No regenerated keys or repeated decryption of unchanged factors. Receipt: narrow-dd-wide-primary-ntt-root-integration-validation.json, original session36119. Not yet production integration/full FBS qualification.

2026-10-11: Saved packed435 worker profile completed4/8/12/4 controls; all13 full ciphertext archives are byte-identical to qualified reference, zero key generation and unchanged archive. Mask times37.5948/30.6774/31.3170/36.6460 ms, automorphism70.0537/63.1721/62.0512/64.5652 ms. Eight-worker mask speedup1.210x against the two4-worker controls' mean, preliminary under concurrent reconstruction/integrity reads. Twelve offers little further benefit. No live reconstruction settings or production defaults changed. Evidence: narrow-dd-packed-workers-profile-validation.json.

Isolated NTT root integration's first attempt returned1 at the diagnostic CLI before any evaluation because the old argc6 guard had not been updated. Preserved failure: narrow-dd-wide-primary-ntt-root-integration-cli-failure-validation.json. Corrected argc8 wrapper compiled; session48196 is evaluating. The retained prehash from the failed CLI-only attempt is reused, with mandatory all-consumed posthash; zero regenerated keys. Production core headers remain unchanged.

Isolated full NTT packed real FBS compiled successfully and is queued in session48159. It uses saved DD433, public auxiliary100 projections at data332/288/244/200/156 and post101, primary64/secondary32 NTT, primary8 conjugation caches and8 workers. It waits for original complex packed validation termination, original real packed qualification and corrected NTT root qualification. Independent all-slot refresh/square gates and all consumed artifact integrity remain mandatory. It does not regenerate input/private/evaluation keys. Receipt: narrow-dd-ll14-packed-precision57-ntt-fbs-validation.json. Compilation/queue are not full NTT FBS qualification or paper speed/security reproduction.

2026-10-11: Corrected isolated NTT root integration evaluator and independent validator both returned0. All37 prepared-input/factor ciphertexts are byte-identical to qualified FFT reference. Independent sparse input and dense final phases check all8192 output slots, max7.99902e-8 below2^-23. Root-stage timings factor114.847 s, alignment13.2294 s, DD-cache4.55e-7 s (prepared once outside root), trees5.81534 s. Final retained integrity is still running in original session48196; this does not yet establish terminal/all_checks_passed. No new private/evaluation keys or production core edits. Evidence: narrow-dd-wide-primary-ntt-root-integration-validation.json.

2026-10-11: Isolated NTT narrow root integration is terminal/all_checks_passed, including independent root7.99902e-8 and all saved key/input/selector/source integrity. No production core changes. Evidence: narrow-dd-wide-primary-ntt-root-integration-validation.json.

Backend move/load lifecycle and73,917,741-byte portable archive compatibility pass with synthetic public-only N512 data. Audit additionally found that explicit FFT preparation and evaluation-key regeneration must detach the optional callback. Second isolated revision fixes staged/default/projected FFT selection and generation reset, constructs the replacement default FFT cache before committing the backend switch, and passes all three selection/released-ownership controls. Regeneration reset placement is source-checked, not runtime key-generation tested. Arithmetic is byte-for-byte unchanged relative to the qualified isolated root after removing lifecycle-only edits. Original production and queued NTT headers are unchanged. Reviewable promotion patch: docs/narrow-dd-ntt-integration.patch, not applied while original tests pin production headers. Evidence: narrow-dd-ntt-backend-lifecycle-validation.json and narrow-dd-ntt-backend-lifecycle-v2-validation.json. No key generation in these checks.

2026-10-11: Full high-precision packed reconstruction is terminal/all_checks_passed with actual child exit0: all65536 terms, all36 retained selector group headers/lengths/hashes, source integrity and unchanged packed435 archive. New key generation0. Extraction5830.4 s (97.173 min), reconstruction including final selector writing6129.8906 s (102.165 min); cached14-depth preparation6.01072 s. Retained selector bytes34,945,894,980. Canonical root44 controller started automatically, reusing the reconstructed groups and existing keys/input. Root evaluation stages completed factor56.4929 s, alignment7.19428 s, DD-cache6.59481 s and trees5.44675 s; independent full factor/root phases are still running. This is reconstruction completion, not full FBS accuracy or paper runtime/security qualification. Evidence: narrow-dd-ll14-packed-precision57-reconstruction-validation.json and narrow-dd-ll14-packed-precision57-root44-validation.json. The lifecycle-v2 production promotion patch passes git apply --check and remains unapplied while original tests pin headers.


### CPU frequency-basis replay with retained selector spectra

The experimental `ckks_narrow_ll14_basis_ntt_fbs` runner evaluates real-input
GeLU, Tanh, and Sigmoid with the Section 7 frequency basis and exact NTT
products. One retained-input run per function measured **40.93, 77.55, and
42.58 seconds**, respectively, on an i9-12900K with `OMP_NUM_THREADS=8`.
S2C and all input-dependent root/series work are timed. Public key loading,
automorphism preparation, and selector-spectrum preparation are setup costs.
These measurements do not establish a controlled speedup against older O2 or
other-machine results, and the remaining 10x performance goal is not achieved.

`narrow-dd-basis-cached-cpu-validation.json` records the runs and independent
integer-phase checks over every coefficient and slot. Each output passes the
fixed 23-bit error gate and a subsequent-square gate. The final ciphertext is
93/44; its square is 49/44. Private keys and the input were reused.

This profile trades memory and new public setup for online time. The selector
spectrum payload is 89,860,534,272 bytes (about 83.7 GiB), plus vector overhead
and temporary workspace; observed process RSS was about 95.7 GiB. The cache
constructor rejects an insufficient explicit budget before loading a group.
It retains only key-dependent data. Online selector reads and transforms are
zero; evaluated factors and roots are recomputed for each function.

The direct selector bank contains 19,044 ciphertexts at 440/57 and uses
continuous Gaussian sigma 102.4 (integer RMS 102.00190833992261), rather than
extending the old low-noise bank to a larger modulus. Generation and writing
took 277.37 seconds once. Its DD433 primary20/secondary20 key uses the retained
dense private key. The four generic dense-LWE screens have minimum estimated
cost 128.968 bits; this is not an overall security qualification. Structured,
ring, correlated-key, and circular assumptions remain unqualified. The saved
packed435 public key is also required for public automorphism derivation.

Build the two opt-in targets, then supply the retained artifacts:

```sh
cmake --build build --target ckks_narrow_ll14_basis_setup ckks_narrow_ll14_basis_ntt_fbs -j2
# Once only; NEW_SELECTORS must not exist. Reuses CHECKPOINT/checkpoint.bin keys.
OMP_NUM_THREADS=8 build/test/ckks/ckks_narrow_ll14_basis_setup CHECKPOINT NEW_SELECTORS full
# OUTPUT must contain the retained candidate-input.bin before starting.
OMP_NUM_THREADS=8 build/test/ckks/ckks_narrow_ll14_basis_ntt_fbs CHECKPOINT SWITCH56 all OUTPUT NEW_SELECTORS NEW_SELECTORS/dd433-primary20.bin PACKED435 --cache-gib 90
```

These are fixed-fixture experimental runners, with checkpoint metadata and
origin guards, not general key import APIs. Reuse the generated selector bank
and DD key on subsequent runs. `--cache-gib 0` selects streamed operation; a
positive budget must fit the full spectrum cache. The budget covers retained
spectra, not total process memory. GCC builds disable IPA-ICF for the replay
because the existing multi-type S2C instantiations require it.

`--tensor-rns 1` uses exact RNS tensors in both the factor product trees and
frequency-basis reconstruction. Both paths reuse the retained relinearization
keys and preserve the existing rescale locations. `--basis-tensor-rns 0`
selects the FFT tensor for reconstruction; `--basis-tensor-pair 1` compares
FFT/RNS/RNS/FFT after a warmup for one requested function. Reports include
`basis_tensor_rns` and `basis_series_seconds`. The paired original-degree Tanh
run reduced reconstruction from 1.275s to 0.357s and total evaluation from
10.453s to 9.551s on the qualified 24-thread CPU configuration; all output
bytes match the independently validated reference. See
[narrow-basis-rns-tensor-validation.json](narrow-basis-rns-tensor-validation.json)
for the complete configuration, checks and timing limits.

For independent verification, set `PYTHONPATH` to the sibling
`Parameter-Selection/python/estimates`, then run
`python3 test/ckks/validate_narrow_ll14_basis.py CHECKPOINT/checkpoint.bin OUTPUT/gelu-outputs.bin OUTPUT/gelu-report.json validation.json`
(and repeat for Tanh/Sigmoid). No key contents are printed.

Lower-scale basis candidates at delta35 failed the fixed accuracy gate
(GeLU error approximately 1.2e-5 to 1.4e-5). The qualified profile retains
delta44 throughout the tree; those failed candidates are not performance wins.


The grouped DD product update tiles 256 NTT points and reduces after at most
three lazy Shoup products. The bound is below 7p and fits uint64 for p<2^61;
other moduli retain the prior calculation. A retained-key N16384 DD433/20/20
benchmark, pinned to eight physical P cores and run before/after in both
orders, measured approximately 6% and 8% faster complete DD applications.
All component outputs were byte-identical. This is a component improvement,
not a whole-refresh speedup estimate.

The repository frequency-basis runner then completed all three functions at
40.35/75.70/39.73 seconds (GeLU/Tanh/Sigmoid), unpinned OMP8, with all saved
ciphertext files byte-identical to the independently decoded reference.
Public setup took 236.52 seconds in this run; it was excluded from the online
timers and reused across functions. Evidence, compiler commands, sampling
scope and the rejected matrix-register experiment are recorded in
`narrow-dd-grouped-ntt-products-validation.json`. The 10x goal remains open.


The CPU basis replay now prepares the public S2C plaintext and Galois spectra
once, before the online timer. `--s2c-cache-gib 9` is the default spectrum
budget; `--s2c-cache-gib 0` uses the evaluation-local path. This budget is
separate from `--cache-gib`, which controls the selector cache. Neither budget
is a cap on total process memory. The retained LL14 profile uses
4,269,801,472 bytes of plaintext spectra plus 734,003,200 bytes of Galois spectra
(about 4.66 GiB combined), excluding small container overhead.

Only high plaintext digits proved zero from the public compact coefficients
are omitted. The prepared plan still supports the entire signed int64 range.
All input-dependent baby rotations, digit transforms, fused products, giant
rotations, sparse-key encapsulation, and rescaling are recomputed online.
The library's existing three-argument coefficient preparation API retains its
evaluation-local behavior; the explicit transform overload permits prepared
public backends without duplicating the lift/switch/rescale pipeline.

The S2C component measured 14.57 seconds on the baseline and
8.17/8.24/8.27 seconds with prepared spectra, after 3.01 seconds of one-time
public preparation. All three coefficient ciphertexts were byte-identical;
28 prepared Galois spectra were reused without rebuilding any online.
`ckks_prepared_compact_linear` covers direct-only, conjugate-only, combined,
full signed-int64 coefficients, omitted zero digits, all-zero plans, budget
rejection, and repeated cache reuse. See
`narrow-dd-prepared-s2c-component-validation.json` for measurement scope and
the pre-existing scalar linear-transform test compilation failure, reproduced
at the baseline commit.


The integrated full replay passed GeLU/Tanh/Sigmoid at
34.96/72.61/35.18 seconds, with S2C at 8.75/8.83/8.78 seconds. All three
ciphertext files, including the remaining-level square, are byte-identical
to the independently decoded reference. The one-time public preparation and
key loading remain separately reported. A live process observation during
Tanh recorded a peak RSS so far of about 103.4 GiB. Full timings are single
unpinned trials and include host variation; the 10x objective remains open.
See `narrow-dd-prepared-s2c-full-validation.json` for the command, complete
reports, output hashes and memory observation. The same saved private keys,
public evaluation keys and input were reused throughout.


The retained LL14 CPU basis runner now defaults to direct small-prime selector
convolutions (`--direct-rns 1`). This replaces nine 52-bit digit convolutions
with a single exact convolution across 17 roughly 30-bit primes. It uses the
same saved selector ciphertexts and evaluation keys. No encryption modulus,
noise distribution, output scale, or rounding rule changes. A compile-time
bound proves that twice the maximum complete signed convolution is smaller
than the CRT product; runtime table checks enforce the public coefficient
bound. The final centered integer is reconstructed exactly modulo the torus.

The selector spectra require about 39.5 GiB, compared with about 83.7 GiB for
the digit backend. `--cache-gib 45` fits the new retained cache; the separate
S2C cache remains controlled by `--s2c-cache-gib 9`. Budgets cover retained
spectra and documented container overhead, not all process memory.
`--direct-rns 0` selects the previous digit backend, and `--cache-gib 0` uses
streamed selectors with either backend. Public setup is reused across
functions in one process; saved key artifacts are reused across process runs.

`modular_rns32` checks the primes, root orders, complete selector encoding and
exact centered CRT against arbitrary-precision integers, including boundary
values and both SIMD and scalar paths. `ckks_narrow_factor_rns32` compares
full-width random-mask head/tail factors byte for byte, rejects inadequate
cache budgets and unsupported table bounds, and verifies zero online selector
reads/transforms. Retained N16384 head0/tail0 four/nine-column component trials
also match exactly. See `narrow-dd-direct-rns-arithmetic-validation.json` and
`narrow-dd-direct-rns-factor-validation.json` for proof and timing scope.


The integrated direct-RNS full replay completed GeLU/Tanh/Sigmoid at
29.59/58.68/29.69 seconds. All saved ciphertext files, including the subsequent
square, match the independently decoded reference byte for byte. The retained
selector allocation is 42,442,068,672 bytes (about 39.53 GiB including vector
objects). Public loading/preparation took 238.47 seconds once and was reused
across all three functions; it is excluded from online timings. A live Tanh
observation recorded VmHWM of 62,084,044 KiB (about 59.2 GiB), compared with
the prior profile's roughly 103.4 GiB observation. These are single unpinned
OMP8 full trials, with host variation. The remaining 10x goal is not achieved.
See `narrow-dd-direct-rns-full-validation.json` for exact reports, command,
output hashes and the memory-observation scope.


Prepared DD products now automatically use two roughly 30-bit NTT primes when
the complete signed digit-convolution bound fits their centered CRT interval
and both primes support the ring degree. Otherwise they retain the 61-bit
backend. The existing primary and secondary gadgets, auxiliary lift and
rounding are unchanged. The fifth `CKKSPreparedDDRelinKeyNTT` template argument
can be set to `false` to force the 61-bit backend for component controls.
Publicly derived automorphism caches and product-tree caches select the new
path through the same existing factories; no key regeneration is required.

Retained DD433/20/20 tests pinned to eight physical P cores, in both orders,
measured about 27% less time per complete DD application. DD spectrum storage
halved from 253,755,392 to 126,877,696 bytes. The independent coefficient tests
cover the small-prime path, forced 61-bit single-prime path, and wide-gadget
two-prime fallback; prepared FFT/NTT basis reconstruction remains byte-identical.
See `narrow-dd-small-prime-dd-component-validation.json`.

The full repository runner completed GeLU/Tanh/Sigmoid at
28.26/55.25/27.43 seconds, with all saved ciphertext files byte-identical to
the independently decoded reference, including the subsequent square.
Public loading/preparation took 97.61 seconds in this run and was reused across
functions. These are single unpinned OMP8 timings; host and filesystem-cache
variation prevent attributing the whole difference to the DD change.
`narrow-dd-small-prime-dd-full-validation.json` records exact timings, cache
selection, output hashes and scope. The remaining 10x goal is still open.


Prepared compact S2C plans now omit exactly zero public branches, giant groups,
and diagonals before allocating their spectra. Valid all-zero operators
explicitly overwrite the result with zero. The retained LL14 source contains
128 completely zero conjugate groups; their rotations and products previously
consumed about 3.4 seconds despite contributing zero. The prepared plan keeps
128 direct groups and no conjugate groups, with byte-identical coefficients.

Prepared FFT batches also use the maximum per-term plaintext digit L1 norm
in each public group, with the balanced top ciphertext-digit bound. A batch
coefficient is bounded by its term count times the input digit magnitude bound
times that public norm. Qualified N512/N16384 profiles use a 2^47 coefficient
budget, giving more headroom than the prior fixed-size worst case. Low-norm
digits may be fused more widely; large-norm digits receive smaller batches.
The optional third prepared-plan constructor argument `false` retains fixed
batching for diagnostics. A 512-term constant-polynomial stress test checks
all coefficients against independent uint128 multiplication and signed
rescaling: the bounded path passes, while the old fixed control has an observed
rounding difference. This synthetic control is separate from the retained
normal-input runs, whose coefficient ciphertexts match in every trial.

The integrated full GeLU/Tanh/Sigmoid replay measures
24.31/51.68/24.09 seconds, with S2C at 4.94/5.42/5.08 seconds. All saved output
ciphertext bytes, including the subsequent square, equal the independently
decoded reference. Public plaintext spectra drop by exactly 1 GiB to
3,196,059,648 bytes; Galois spectra remain 734,003,200 bytes. A 4 GiB public
S2C spectrum budget is demonstrated (`--s2c-cache-gib 4`), separate from the
45 GiB selector-cache budget. The runner releases the copied compact source
plan after prepared spectra own everything they need. A live Tanh observation
records VmHWM of about 54.9 GiB. Single unpinned OMP8 full timings include host
variation, and the 10x objective remains open. Exact reports and measurement
scope are in `narrow-dd-pruned-bounded-s2c-full-validation.json`; component
controls and the stress-reference correction are in
`narrow-dd-pruned-bounded-s2c-component-validation.json`.
