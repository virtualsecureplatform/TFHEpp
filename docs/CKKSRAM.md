# CKKS batch refresh for a TLWE RAM

This CPU prototype keeps the existing TFHE RAM's circuit bootstrapped
address, CMUX read tree, and conditional write. It replaces the refresh of
every updated cell with one CKKS Boolean bootstrap batch. It is an
experimental implementation, with no end-to-end security certification.

`RefreshableRAM<KS, CB, PKS, BR, AddressBits, WordBits>::access` takes an
array of address-bit TLWEs, a write-enable TLWE, and a word of data TLWEs.
It returns the **old word before the write**, also as TLWEs. All use the
same TFHE secret and standard `+/-mu` Boolean convention. There is no
cleartext address in the evaluator. The public refresh callback sees the
updated values of every cell, regardless of the encrypted write-enable.

The shared implementation supports two refresh callbacks:

1. TFHE: identity key switch, blind rotation, and sample extraction per
   cell, matching the algorithm in `test/tfhe/ram.cpp`.
2. CKKS: functional coefficient packing, one narrow Boolean bootstrap,
   slots-to-coefficients decoding, and direct LWE key switching back to
   TFHE. Every access refreshes the whole RAM; deferred refresh is not
   assumed in this comparison.

Only coefficient zero of a cell's TRLWE carries a meaningful value.
Inverse sample extraction can therefore install a refreshed TLWE directly;
it does not need another TFHE bootstrap. CMUX reads/writes do not rotate
the cell polynomial. Other coefficients are unspecified and must not be
read as additional memory bits.

## Bridge and bit restoration

The packing key encrypts each TFHE secret coefficient, multiplied by four
radix-256 weights, as a constant polynomial under the bootstrap's actual
structured sparse CKKS secret. Public digit polynomials place the TLWE
phases into consecutive coefficients. Exact two-prime negacyclic NTT/CRT
arithmetic computes this linear map. No padded TFHE ring secret is used.
The arithmetic regression compares both ciphertext components against a
direct integer negacyclic convolution, including wrap and unfilled slots.

The packing level is q=2^48. A 32-bit TLWE is lifted by 16 bits; its
Boolean phase `+/-1/8` becomes `x=+/-2` at CKKS scale 2^44. The restoring
function is

```
sin(4*pi*x/16).
```

It equals `+/-1` at the two Boolean values and has zero first derivative
there. A frequency-two narrow blind rotation with weight `-i/2`, added to
its conjugate, evaluates this function. This is a real encrypted CKKS
bootstrap, not an identity refresh or a client decrypt/re-encrypt step.
The existing retuned selector, DD, and native automorphism keys are reused.

The result is decoded into coefficients under the dense CKKS secret.
Rescaling by 15 bits and reducing the active modulus to 32 bits gives
scale 2^29, restoring the TFHE `+/-1/8` phase convention. The return key
uses ordinary TFHE TLWE encryptions (`BaseBit=2, Levels=16`) with the
target parameter's noise. Its payload is approximately 2 GiB. The forward
packing key is approximately 1 GiB; its prepared NTT cache is approximately
4 GiB. These are additional to the retained CKKS bootstrap caches.

The benchmark owner decrypts only for validation, outside the timed RAM
access. The refresh callback and RAM access API accept no private key.

## Run and reuse

Build from the TFHEpp root:

```sh
cmake -S . -B build
cmake --build build --target ckks_ram_setup ckks_ram_fixture_check ckks_tlwe_coefficient_packing ckks_narrow_ll14_retuned_ram -j2
build/test/ckks/ckks_tlwe_coefficient_packing
```

Generate the RAM fixture **once**, using the checkpoint matching the
retained CKKS bootstrap keys:

```sh
build/test/ckks/ckks_ram_setup CHECKPOINT_DIRECTORY NEW_RAM_FIXTURE_DIRECTORY
```

Setup refuses an existing destination. Repeated evaluations reuse
`owner.bin`, `tfhe-eval.bin`, `pack48.bin`, `to-tfhe32.bin`, and
`trace-4x2.bin`, and `trace-256x1.bin`. Private owner files are test fixtures and must stay private.
The evaluator generates no keys or fresh input ciphertexts.
Owner deserialization uses a constructor-free archive reader: TFHE's
ordinary `SecretKey` constructor would otherwise generate discarded
secrets before loading the retained ones. The separate
`ckks_ram_fixture_check RAM_FIXTURE_DIRECTORY` validates the reader
against both saved encrypted traces and their initial memories.

Use the retuned bootstrap runner's public-key arguments and performance
options, selecting the RAM executable and fixture directory:

```sh
OMP_NUM_THREADS=8 build/test/ckks/ckks_narrow_ll14_retuned_ram \
  CHECKPOINT_DIRECTORY SWITCH56_FILE tanh RAM_FIXTURE_DIRECTORY \
  RETUNED_SELECTOR_DIRECTORY RETUNED_DD433_FILE unused-with-native-auto \
  --native-auto-directory RETUNED_NATIVE32_DIRECTORY \
  --cache-gib 45 --s2c-cache-gib 4 \
  --direct-rns 1 --s2c-rns 1 --tensor-rns 1 --auto-rns 1 \
  --parallel-conjugation 1 --parallel-tail-alignment 1 \
  --batch-product-trees 1 --dd-wide-rns 1 --dd-primary-bits 40 \
  --dd-level-chunks 1 --matrix-joint 1 --parallel-table-ntt 1 \
  --s2c-baby-workers 8 --s2c-group-workers 8
```

`tanh` is the inherited runner's required profile selector; the RAM path
always executes the frequency-two Boolean function above. This prototype
currently requires fixed primary-4 S2C caches.

Set `CKKS_RAM_PROFILE=256x1` in the environment to run the 256-word,
one-bit profile matching the dimensions in the existing TFHE RAM test.
Its report is `comparison-256x1.json`. Without that variable the runner
uses four words of two bits and writes `comparison.json`.
To add missing traces to an older fixture without regenerating any keys,
run `ckks_ram_setup --add-traces RAM_FIXTURE_DIRECTORY`. Existing trace
archives are also reused.

The first trace has four words of two bits and 16 sequential accesses. It
visits every address, writes both bit values, disables writes, overwrites
the same address, and feeds a prior encrypted read output into a write.
It checks the read result and **every memory cell after every access**.
Both backends start with the same saved encrypted memory and requests.
`comparison.json` contains total access and individual stage timings.
Setup/cache preparation and validation are excluded from access timings.

## Paired crossover experiment

**No end-to-end crossover was observed at the tested sizes**, including a
full 8,192-bit batch. The earlier 600–700-bit extrapolation from the initial
timings is not supported by these controlled measurements.

| Stored bits (one-bit words) | TFHE mean access | CKKS mean access | CKKS / TFHE |
| ---: | ---: | ---: | ---: |
| 256 | 2.255 s | 6.659 s | 2.953x |
| 512 | 4.184 s | 8.500 s | 2.032x |
| 1,024 | 7.928 s | 11.948 s | 1.507x |
| 8,192 | 60.175 s | 61.235 s | 1.018x |

TFHE was faster in **every one of the four measured pairs at each size**.
At full capacity, TFHE access times ranged from 60.143 to 60.227 seconds;
CKKS ranged from 61.178 to 61.292 seconds. All 40 operations, including
warm-ups, passed: 99,840 post-access cell-bit checks and 40 read-bit checks.
See [the complete receipts and provenance](ckks-ram-crossover.json).

At 8,192 bits the CKKS refresh stages averaged 0.262 seconds for packing,
3.750 seconds for the Boolean bootstrap, 0.604 seconds for decoding, and
**54.866 seconds for conversion back to TLWE**. Complete CKKS refresh
averaged 59.486 seconds versus 58.427 seconds for TFHE. The CKKS core
amortizes well; its per-cell output conversion consumes almost all of
that benefit. Core-bootstrap timings alone do not establish a RAM crossover.

These results supersede the initial comparison below for scaling decisions.
The old measurements had different runtime/load conditions and much slower
TFHE timings. Their causes were not isolated. All sizes above were measured
under the same explicit affinity policy; the older and newer numbers must
not be combined into one fitted curve. This finite experiment does not
exclude a crossover at every untested size or under other implementations
and hardware. It does not support claiming one for the current RAM path.

Use `CKKS_RAM_PROFILE=crossover` for 256-, 512-, and 1,024-bit RAMs and
`CKKS_RAM_PROFILE=crossover-full` for 8,192 bits. Both modes use one-bit
words. Add missing encrypted traces once with
`ckks_ram_setup --add-traces RAM_FIXTURE_DIRECTORY`; all existing keys
and input archives are reused. The full case fills one CKKS batch.

These modes keep two independent RAM states, starting from identical
saved ciphertexts. They apply five identical encrypted requests to each
backend, alternate which backend runs first, and exclude the first access
of each backend as warm-up. The remaining four paired accesses include a
write to the final address and disabled writes reading the middle and
the two endpoints. Every read and every cell after every access is checked
outside the timed operation. Stage timings include both conversions.

The controlled runs use `OMP_NUM_THREADS=8`, `OMP_DYNAMIC=FALSE`,
`OMP_PROC_BIND=close`, and
`OMP_PLACES='{0},{2},{4},{6},{8},{10},{12},{14}'`, with process affinity
`taskset -c 0-15`. The TFHE refresh loop remains serial. Set the process
stack limit to 128 MiB (`ulimit -s 131072` in a subshell) before running
the full case: the inherited `RAMUX` allocates its large temporary tree
on the stack. This changes no cryptographic parameter or gate algorithm.

Each profile writes `crossover-SIZEx1.jsonl` with individual operation
receipts and `crossover-SIZEx1.json` with its summary. Validate and collect
them with:

```sh
python3 test/ckks/analyze_ram_crossover.py RAM_FIXTURE_DIRECTORY \
  --output crossover-results.json
```

The analyzer requires all ten successful receipts per profile and checks
the order, warm-up classification, and agreement with summary timings.
It reports mean ratios and individual paired differences, distinguishing
a mean win from a win in every measured pair. A finite experiment cannot
prove that a crossover exists or is absent for every machine and size.

## Initial measured comparison

The [machine-readable record](ckks-ram-comparison.json) contains complete
timings, commands, source hashes, and validation scope. On an i9-12900K,
with eight CKKS threads pinned to P cores:

| RAM | Accesses per backend | TFHE mean access | CKKS mean access | CKKS / TFHE |
| --- | ---: | ---: | ---: | ---: |
| 4 words x 2 bits | 16 | 1.006 s | 11.742 s | 11.67x |
| 256 words x 1 bit | 8 | 11.078 s | 17.353 s | 1.57x |

Both traces passed every read and every post-access memory-cell check:
48 complete operations and 4,352 cell-bit checks across the two backends.
The small trace also passed encrypted read-output feedback into a write.
The packing phase audit errors were 1.16e-8 and 2.09e-8 in CKKS message
units. The largest CKKS-backed memory error from the ideal `+/-1/8` TLWE
phase was 0.01020 in the 256-word trace; these are Boolean correctness
tests, not a statistical failure-probability qualification.

The mean CKKS refresh-stage times were:

| Stage | 8 cell bits | 256 cell bits |
| --- | ---: | ---: |
| TLWE coefficient packing | 0.509 s | 0.523 s |
| CKKS Boolean bootstrap | 8.934 s | 8.900 s |
| Slots-to-coefficients decoding | 1.416 s | 1.396 s |
| Direct switch back to TLWE | 0.145 s | 4.025 s |

The fixed bootstrap cost amortizes with batch size, while direct output
switching grows with the number of cells. **There is no measured speed
win yet.** The larger case closes much of the small-batch gap but does
not establish a crossover point. The combined benchmark process peaked
at about 55.8 GiB RSS, including retained CKKS caches and both schemes'
evaluation material. This is not an isolated TFHE memory measurement.

The TFHE callback follows the existing serial per-cell refresh loop;
this is not a comparison against a newly parallelized TFHE implementation.
The host had unrelated CPU load. The RAM adapter uses retained RNS
selectors, batched product trees, and independent native automorphisms;
the inherited replay runner's windowed-DD/native-auto-batch trial flags
are not applied by this adapter. These are prototype measurements rather
than best-tuned throughput claims.

The final owner-loading correction described above was independently
checked after these runs. It changes excluded setup work and leaves the
timed RAM and refresh arithmetic unchanged.

## Limits of the assessment

The CKKS core uses the
[retuned classical-128 target profile](narrow-ll14-retuned-profile.md),
including its stated heuristic models and circular/KDM assumptions. The
new bridge adds q=2^48 sparse-secret packing keys and ordinary destination
TFHE LWE switching keys. These additional public key families and their
cross-scheme circular dependencies have not received a separate complete
Parameter-Selection assessment. A successful RAM trace does not extend
the core's security qualification to the whole bridge.

The saved dense/sparse development secrets also have historical low-noise
evaluation keys in offline diagnostics. They are not deployment secrets.
This prototype does not use the rejected generic scheme-switching profile.

A four-word trace measures small-RAM latency. The additional 256-word
trace checks eight sequential accesses at addresses 0, 127, 128, and 255,
including writes and disabled writes; every cell is checked after each
access. These two profiles cannot establish the
crossover point for hundreds or thousands of batched bits, a failure
probability bound, or a safe number of accesses between refreshes. The
packing and direct output switching are straightforward prototype paths;
include their cost when comparing with TFHE.

## Experimental subring output switch

**Security status:** investigation found that the historical v1 generator used
noncryptographic mask expansion. It has been removed. Current v2 ring keys use
domain-separated BLAKE3 expansion with the full 256-bit seed; v1 ring archives
and legacy seeded CKKS rows are rejected. Generate replacement public keys in
a separate directory. Earlier v1 timings remain historical observations;
circular security is accepted for the current assessment, while the sigma128
family misses the strict classical screen. See the [sigma128 assessment](#sigma128-assessment-with-circular-security-assumed), the
[security investigation](#security-investigation-2026-10-25) and its
[reproducible receipt](ckks-ram-ring-switch-security.json).

The opt-in `ring` backend moves the output key switch before sample extraction.
It switches the q32 coefficient ciphertext from the dense degree-16384 CKKS
secret to the existing degree-1024 TFHE secret embedded as `t(X^16)`. Splitting
by coefficient residue modulo 16 produces sixteen ordinary small-ring TRLWEs;
coefficient `r + 16*j` becomes coefficient `j` in small ciphertext `r`.
The existing first-8192-coefficient layout uses 512 coefficients of each small
ciphertext. Ordinary sample extraction then returns TLWEs under the original
TFHE secret, at scale 2^29, without another per-bit key switch.

This is a separate construction from the older contiguous-zero-padding
`CKKSToTFHERingSwitchKey` API. The new public API is
`CKKSToTFHESubringSwitchKey32<P,L>` with
`PreparedCKKSToTFHESubringSwitch32<P,L>::switchAndSplit` and `apply`.
It requires single-polynomial secrets, a power-of-two destination degree
that divides the source degree, a multi-limb source through degree 32768,
and a 32-bit destination torus. `apply` extracts a nonempty prefix of at most
`P::n` coefficients in the original order. RAM packing remains limited to
`P::n/2` values. Prepared spectra are immutable; per-call scratch is local.

The bridge uses 16 complete balanced radix-4 rows, q=2^32, and continuous
integer Gaussian sigma 256 (normalized sigma 2^-24). No auxiliary ciphertext
modulus, gadget truncation, or floating FFT convolution is used. The existing
RNS32 gadget implementation provides exact convolution with two NTT primes.
Raw row payload is 2 MiB; prepared spectra occupy 4 MiB. These are payload
counts, not isolated process RSS measurements.

### Add public keys to a retained RAM fixture

Build the normal 128-bit topology (`USE_BLOCK_BINARY=OFF`, `ENABLE_TEST=ON`).
The new archive fingerprints use the default BLAKE3 backend.

```sh
cmake --build build --target ckks_ram_setup ckks_subring_switch \
  ckks_ram_output_switch_benchmark ckks_narrow_ll14_retuned_ram -j2
OMP_NUM_THREADS=8 build/test/ckks/ckks_ram_setup --add-ring-switch-keys \
  CHECKPOINT_DIRECTORY RAM_FIXTURE_DIRECTORY
```

This mode reads the retained checkpoint and constructor-free owner archive.
It generates no private keys, trace requests, or fresh input encryptions. It
adds `to-tfhe32-subring-sigma256-v2.bin` and
`to-tfhe32-direct-sigma256-v1.bin`, retaining the original sigma-128
`to-tfhe32.bin`. The second new file is a direct-switch baseline at matched
bridge-key noise. The original three-argument `CKKSToTFHEKeyGen` keeps its
sampler; its new explicit-noise overload accepts normalized Gaussian sigma.

Both new archives validate a versioned profile, source-checkpoint and owner
BLAKE3-256 fingerprints, degrees, embedding stride, modulus, scale, radix,
row count, and sigma before reading their fixed-size payloads. Setup reuses
matching files, rejects incompatible files, and publishes new files without
overwriting an existing destination. A stale `.lock` or `.partial` file from
an interrupted setup must be inspected and removed explicitly before retrying.
Key generation rejects source parameter types whose noise floor would silently
raise the fixed sigma-256 profile.

`ring-switch-key-audit.json` records an independent schoolbook small-ring phase
audit over all 262144 ring-key error coefficients, plus a distributed 4096-row
sample of the matched direct key. It publishes only aggregate residuals and
fingerprints. The ring audit checks the actual uniform added-error bound
`2*sum(abs(key_error))/2^32` against the assessed limit 0.025746934137548232.

### Select and compare the bridges

Apply these environment variables to the existing retained RAM command above:

- `CKKS_RAM_OUTPUT_SWITCH=direct` selects the original sigma-128 direct key
  and remains the default; `direct256` selects the matched direct key;
  `ring` selects the prepared subring bridge.
- `CKKS_RAM_REPORT_DIRECTORY=NEW_DIRECTORY` separates output receipts from
  the retained fixture. Existing report filenames remain unchanged for normal
  and legacy crossover profiles.
- `CKKS_RAM_PROFILE=ring-compare` runs 256-, 512-, 1024-, and 8192-bit profiles
  with four independent states: serial TFHE and all three CKKS output bridges.
  Each profile applies the same five retained requests, excludes the first
  access as warm-up, and rotates execution order so each backend occupies
  each measured position once. Every read and every post-access cell is
  checked outside timing.
- `CKKS_RAM_PROFILE=crossover-full` with `CKKS_RAM_CAPTURE_ONLY=1` derives
  one public q32 input from the retained encrypted initial memory. It checks
  all restored bits and saves `output-input32.bin` without encrypting a new
  input. This mode is for isolating the output-switch cost.

All new validation paths audit the pre-switch phase using exact integer
arithmetic and enforce the assumed 1/16 error budget for the used coefficients.
Input audits and key preparation are outside RAM access timings. Reports
separate packing, Boolean bootstrap, decoding, and output conversion;
ring kernel/split/extraction timers in bridge metadata are cumulative and
include warm-ups. The standalone benchmark records these timers per trial.
Combined RSS includes all loaded baselines and retained bootstrap caches;
it must not be described as the ring-only memory footprint.

A reproducible driver uses the retained paths/flags in the paired crossover
receipt, pins eight P cores, sets the 128 MiB stack limit, hashes retained
files before/after, and samples external CPU activity:

```sh
python3 test/ckks/tools/run_ram_ring_switch.py \
  --binary build/test/ckks/ckks_narrow_ll14_retuned_ram \
  --bridge-benchmark build/test/ckks/ckks_ram_output_switch_benchmark \
  --output NEW_EXPERIMENT_DIRECTORY
python3 test/ckks/analyze_ram_ring_switch.py NEW_EXPERIMENT_DIRECTORY \
  --output ring-switch-results.json
```

`--baseline-receipt` can select a receipt with paths valid on another machine.
The default driver reuses one warmed bootstrap cache across the capture, small
traces and four-way comparison. It has
`--phase pipeline|capture|small|comparison|bridge` for staged execution;
phase directories must be fresh. The small phase checks all 16 accesses of
the 4x2 trace, including encrypted read-output feedback, for each bridge.
The bridge phase compares the same saved q32 ciphertext at all four batch
sizes: one warm-up per backend, then all six permutations of execution order.
The analyzer requires complete successful journals, matching stage sums,
all cell/read checks, and retained-file integrity receipts. It reports
stage-only and whole-access gains separately; contended observations do not
establish an uncontended or universal speedup.

### Qualification boundary

[The parameter assessment](ckks-ram-ring-switch-assessment.json) screens the
new sigma-256 family at the effective degree 1024, not 16384. The weakest
classical generic LWE estimates are 135.19 bits (BDGL16) and 133.09 (MATZOV).
These use variance-matched proxies for truncated continuous Gaussian noise.
The candidate sigma-128 family reaches only 127.79 under MATZOV. Adding the
new key does not qualify or repair other historical lower-noise public material.

The analytical 0.02575 bridge-error bound is conditional on independent ideal
Gaussian key sampling and exact arithmetic. Combined with a separately bounded
1/16 input phase error, it stays below the 1/8 Boolean decision margin. A
successful trace verifies this budget for those observed inputs, not for every
possible RAM execution. Small-ring RLWE structure, complete upstream noise,
finite-precision sampling, and cross-scheme circular/KDM security remain
separate qualification tasks. The prototype stays opt-in and all reports keep
`security_validated: false`.

### Observed prototype run (2026-10-25)

The retained-fixture run used four independent RAM states, rotated execution
order, one warm-up, and four measured accesses per backend. Mean whole-access
times in seconds were:

| Bits | Serial TFHE | Direct σ128 | Direct σ256 | Ring σ256 |
| ---: | ---: | ---: | ---: | ---: |
| 256 | 2.315 | 6.876 | 6.864 | 5.215 |
| 512 | 4.274 | 8.677 | 8.682 | 5.333 |
| 1024 | 8.120 | 12.100 | 12.076 | 5.416 |
| 8192 | 61.487 | 60.176 | 60.134 | 6.643 |

At 8192 bits, ring switching reduced whole-access time by
9.05× against the matched-noise direct bridge and
9.26× against serial TFHE in this run.
The isolated output-conversion benchmark changed from
65.5354 s to 0.011118 s
(5895×). The remaining dominant CKKS stage was Boolean
bootstrapping, averaging 3.961 s per full-batch access;
output conversion averaged 0.011225 s in the RAM comparison.

The ring key payload was 2 MiB and its prepared cache 4 MiB; the direct key
payload was 2050 MiB. The evaluator releases the raw ring key after preparing
its cache. Combined process peak RSS includes the retained bootstrap cache,
all loaded baselines, RAM states and scratch; it is not ring-only RSS.

All 80 four-way RAM accesses, 84 isolated conversion trials (including
warm-ups), the full 8192-bit capture, and each bridge's 16-access small-RAM
trace passed. The small traces include encrypted read-output feedback.
Retained-file hashes were unchanged. The generated ring key's independently
audited uniform added-error bound was 0.024777563754469156, below the assessed
0.025746934137548232 limit. All observed pre-switch errors stayed below 1/16.

The strict whole-job external-CPU screen did not pass. These are observations
from this run and do not establish a universal or uncontended speedup.
The default remains the original direct bridge, and complete protocol security
remains unqualified. [The validation receipt](ckks-ram-ring-switch-validation.json)
contains raw timings, correctness receipts, source and binary provenance,
key audits, negative analyzer checks, and eight freshly reproduced
Parameter-Selection estimates for the implemented sigma-256 profile.

### Security investigation (2026-10-25)

The v1 prototype is **not qualified for deployment**. The investigation found
a concrete mask-generation defect, in addition to the existing parameter and
composition qualifications. No secret-key recovery was demonstrated. No private
keys or evaluation-key payloads were read, and no retained material was changed
during this investigation. The earlier performance and correctness receipts
remain historical observations.

The algebraic construction has a useful precise interpretation. Write
`F(X) = sum_r X^r F_r(X^16)`, and let `Y=X^16`. A public row
`B_j=A_j*t(X^16)+g_j*s+E_j` splits into sixteen equations
`B_(j,r)=A_(j,r)*t+g_j*s_r+E_(j,r)` in
`Z_(2^32)[Y]/(Y^1024+1)`. Deinterleaving merely permutes coefficients: it
does not divide noise by 16 or lower the secret dimension below 1024.
With independent source/target secrets and uniform masks, giving away `s`
and subtracting its messages leaves 256 ordinary small-ring RLWE polynomial
samples. This establishes the appropriate *conditional* assumption for the
one-way bridge. It does not establish the security of cyclic auxiliary keys.
The small-ring reduction in [Gentry, Halevi and Smart, Section 3.1](https://people.csail.mit.edu/shaih/pubs/ring-switching-old.pdf)
is a related precedent; that preliminary version's odd-modulus theorem is
not asserted verbatim for our even-modulus torus implementation.

The historical v1 implementation failed the uniform-mask premise. Its call to
`exactEncryptPolynomialAtLevel` uses `CKKSSeedStream`, which folds a 256-bit
seed into one 64-bit SplitMix state. `USE_BLAKE3` does not replace this expander.
The mixer is invertible, so two distinct seeds with identical streams can be
constructed directly. A C++ probe verified 16384 identical consecutive outputs
for such a pair. At q32, enumerating the missing upper 32 bits of the first
output and checking subsequent words gives a mask distinguisher requiring
at most 2^32 candidate states. The executable Python demonstration searches
a deliberately restricted 2^12 range, not all 2^32 candidates. These are
mask-distribution results: masks are already public, so they do **not** prove
a 2^32 secret-key attack or a 64-bit cap on encryption security. They do rule
out treating this expander as a 128-bit cryptographic PRG.

An extended Parameter-Selection run completed all 36 requested estimates:
seven attacks per classical cost model at each noise setting, plus four
finite-sample controls per model for sigma256. The finite control uses
262144 flattened scalar coordinates as a generic LWE heuristic, not a claim
that ring rotations are independent scalar samples.

| Conditional uniform-mask profile | BDGL16 minimum | MATZOV minimum |
| --- | ---: | ---: |
| n1024, q32, sigma256, seven attacks | 135.19 bits | 133.09 bits |
| n1024, q32, sigma128, seven attacks | 129.65 bits | 127.79 bits |
| sigma256, finite scalar-coordinate control | 135.19 bits | 133.09 bits |

The added attacks were MITM and two primal hybrid variants. None reduced the
previous minima. These are classical estimates conditional on the mask and
noise models; no quantum or complete ring-specific qualification was obtained.
Truncating a continuous Gaussian also differs from the variance-matched
discrete Gaussian used by the estimator. For sigma256 the ideal per-coefficient
total-variation distance is approximately 0.00169037, so equality of variance
does not supply a negligible-distribution-distance argument.

The public-key graph has further obligations: packing encrypts TFHE secret
coefficients under the sparse CKKS secret; selectors encrypt functions of that
sparse secret under the dense CKKS secret; output switching encrypts the dense
secret under the TFHE secret. This cycle already exists with direct output
switching. The dense-to-sparse switch, automorphism keys and relinearization
material add further secret-dependent relationships. Ordinary independent-key
IND-CPA security does not automatically justify the whole graph. A compatible
circular/KDM assumption or a separate composition proof is required.

The inspected full-mask packing, direct-switch, native automorphism, selector
and compact-DD setup paths do not use this particular legacy seeded expander.
However, sigma128 TFHE encryptions and evaluation material remain separate
128-bit-screen blockers. Raising only the new output key's noise cannot repair
them. If weaker material was previously disclosed, replacement public keys
cannot revoke that information; affected deployment secrets must be replaced.
The existence of an old private checkpoint alone does not imply disclosure.

The required repair is full masks from the existing CSPRNG with explicit
exact encryption, or a separately reviewed cryptographic expansion format.
The ring archive already stores full masks, so this repair need not enlarge
the 2 MiB payload or change evaluation arithmetic. A new archive version must
bind the mask-generation law and distinguish repaired keys from v1. Then
generate replacement material in a separate directory, repeat noise/correctness
and performance checks, and assess every exposure under the shared secrets.
Finite key layers can avoid some inter-epoch cycles at additional setup cost;
they do not by themselves prove repeated fixed-key RAM or self-key security.

The independently audited fixed-key bridge-error bound of 0.02477756 still
holds for the tested key irrespective of the Gaussian model. Universal
upstream RAM noise and finite-precision sampler guarantees remain open.
[The security receipt](ckks-ram-ring-switch-security.json) contains the complete
equations, source hashes, C++ probe, numerical sampler comparison, key graph,
all attack outputs, and remediation order. Reproduce the checks with:

```sh
python3 test/ckks/tools/assess_ram_ring_security.py --output algebra.json
sage -python test/ckks/tools/assess_ram_ring_security.py --mode estimator \
  --selection ../Parameter-Selection --output extended-estimates.json
```

### Removal of noncryptographic expansion

The current implementation removes both CKKS SplitMix expansion and BFV's
SplitMix fallback. A shared `SeededMaskStream` uses keyed BLAKE3 with a complete
256-bit seed, explicit little-endian serialization, separate CKKS/BFV domains,
and buffered XOF output. It throws if BLAKE3 support is disabled; there is no
noncryptographic fallback. `BLAKE3PRNG` also no longer accepts an implicit zero
seed or a single 64-bit seed. Supply OS-backed randomness or a full 256-bit seed.

Seeded CKKS rows carry a fixed magic, dimensions, modulus and algorithm tag
before their payload. Legacy untagged rows are rejected. Ring key magic and
filename advance to v2; the unaffected direct256 and captured-input formats
remain compatible. BFV seeded transition keys use format 3 and reject format 1;
its full-mask format 2 is unchanged. Retained historical files are not migrated
or overwritten: regenerate the affected keys rather than changing their tags.
The historical security receipt describes the defect before this removal.

Deterministic noncryptographic generators remain in tests and numerical probes
for public test vectors. They are not used by production encryption/key-generation
implementations. Circular security is accepted as an assumption for the current RAM assessment.
The remaining sigma128 public-material screen is assessed below.

### Sigma128 assessment with circular security assumed

For this assessment, circular security of the RAM evaluation-key graph is an
accepted assumption. It is not treated as a remaining parameter-screen blocker.
At the active level-1 parameters (`n=1024`, `q=2^32`, iid ternary secret,
`alpha=2^-25`), sigma128 is borderline approximately 128-bit classical strength,
but fails a strict requirement that both tested models return at least 128 bits.
This small model-dependent shortfall is not a demonstrated practical attack.

The complete seven-attack grid gives:

| Continuous sigma | BDGL16 minimum | MATZOV minimum |
| ---: | ---: | ---: |
| 128 | 129.65 bits | 127.79 bits |
| 132 | 129.94 bits | 128.06 bits |
| 136 | 130.23 bits | 128.33 bits |
| 144 | 130.52 bits | 128.64 bits |
| 160 | 131.39 bits | 129.33 bits |
| 192 | 132.85 bits | 130.74 bits |
| 256 | 135.19 bits | 133.09 bits |

These use Parameter-Selection's pinned estimator, BDGL16 and MATZOV classical
costs, max-beta 4096, and variance-matched discrete-Gaussian proxies. The actual
sampler truncates a continuous Gaussian: sigma128 gives ideal integer RMS
127.601738. Even using nominal discrete-Gaussian RMS128 instead gives a
MATZOV dual-hybrid estimate of 127.795 bits, so the shortfall is not an artifact
of the small RMS correction. The seven attacks comprise uSVP, BDD, dual,
dual-hybrid, MITM, and two primal hybrid variants. See the estimator's
[dual](https://lattice-estimator.readthedocs.io/en/latest/algorithms/lwe-dual.html)
and [primal](https://lattice-estimator.readthedocs.io/en/latest/algorithms/lwe-primal.html)
model documentation.

There is enough public material for the finite-sample estimate. In
`privksk4cb_1`, the function polynomial is `{1}`. Each public row has the form
`b(X)=a(X)*t(X)+c+e(X)`, where `c` is a constant even though its value depends
on the source secret. Extracting coefficient 1 gives an ordinary known-zero
LWE sample under `t`. The mask is a signed permutation of the original uniform
mask coefficients. Taking one fixed coefficient from each independent row
avoids treating correlated ring rotations as independent samples. The private
switch has `(2048+1)*8*4 = 65568` such independent rows. The finite-m65568
screen has exactly the same minima as the unbounded screen; the reported
optimal dual-hybrid instance uses only 1024 samples. The 8192 initial-memory
TRLWEs also have constant messages and supply independent known-zero
coefficient-1 samples. Restricting the application to few input requests, or
retiring only the old direct bridge, therefore does not remove this exposure.

Independent schoolbook integer audits of synthetic encryptions checked 65536
coefficients at each of sigma128 and sigma256. The observed RMS values were
127.50 and 256.29, consistent with their intended noise. No multiplication
rounding residual was observed on those 131072 coefficients of the tested FFT
backend; this empirical check does not prove a universal rounding bound.
A small-domain private-switch probe confirmed the known-zero extraction and
another important detail: setting `PKS::alpha` for sigma256 while leaving its
**target** at sigma128 still produces sigma128 noise. `privkskgen` calls
`trlweSymEncryptZero` with the target parameter's default noise.

I recommend sigma256 (`alpha=2^-24`) across the affected level-1 family. Its
5.09-bit MATZOV margin is preferable to sigma132's 0.057-bit margin, and it
matches the current ring and matched-direct bridge profiles. Regeneration
must cover ordinary level-1 encryption, the level-1 bootstrapping key,
private-switch keys, and any direct-switch material under that secret.
Changing only the new ring key does not fix the other public exposures.
Gate correctness, circuit bootstrapping, RAM feedback and upstream CKKS noise
must be checked again after a retune. No encryption parameters were changed
by this assessment.

A retuned fixture also needs explicit noise-profile/version binding.
`portablelvl1param` formerly left alpha uninitialized and omitted it from its
standalone serialization and equality. The aggregate `lweParams` serializer did
write the alpha field, but did not initialize or validate it. The original RAM
fixture header did not bind the noise profile either. A source-constant
change is therefore insufficient to reject old weaker archives. If weaker
material was disclosed for a deployment secret, regenerate under a fresh
secret as well: replacing public evaluation keys cannot revoke old exposure.
The retained diagnostic fixtures remain untouched.

For context, the other active TFHE levels pass the same seven-attack screens:
level 0 (`n=630`, q16, binary secret, integer RMS5.6783) has minima
130.68/128.51 bits; level 2 (`n=2048`, q64, ternary secret, integer RMS8191.60)
has minima 131.54/129.62 bits. Level 0 becomes the smallest screened margin
if level 1 is raised to sigma256. This is a classical parameter assessment
under the stated assumptions, rather than an exact-sampler or quantum proof.

[The assessment receipt](ckks-ram-sigma128-assessment.json) contains all 168
completed estimates (166 grid/control results plus two nominal-RMS controls),
finite-sample derivations, source hashes and the synthetic-noise probe. Reproduce
the combined grid with:

```sh
sage -python test/ckks/tools/assess_tfhe_sigma128.py \
  --selection ../Parameter-Selection --output tfhe-sigma128-estimates.json
```


### Opt-in sigma256 prototype

Build the default 128-bit parameter family with `-DUSE_BLOCK_BINARY=OFF
-DUSE_LVL1_SIGMA256=ON`. The option is off by default and rejects conflicting
parameter families. It raises level-1 fresh integer noise to sigma256
(`alpha=2^-24`), uses four base-32 gadget rows, and uses twelve base-4
level-2-to-level-1 private-switch rows. The smaller gadget digits control the
noise accumulated by the full RAM tree. Merely doubling sigma with the old
three base-64 / eight base-8 decompositions passed the synthetic and small RAM
checks but failed the full-size CKKS input budget; that attempt is retained in
the qualification receipt. The pre-switch budget remains 1/16.

New TFHE parameter archives bind dimensions, noise, gadget decompositions,
key-switch parameters, sampler and mask versions before loading evaluation
material. RAM fixtures additionally bind the same profile. The sigma256 build
rejects original unframed archives and incompatible framed profiles. The
default build can read its legacy archives as explicitly unbound material,
but refuses to re-save them under the new binding. `portablelvl1param` now
initializes, serializes and compares alpha.

For a diagnostic retune of the original unframed fixture, retaining its private
owners while regenerating all TFHE evaluation material and input ciphertexts:

```sh
build/test/ckks/ckks_ram_setup --retune-lvl1-sigma256 \
  CHECKPOINT ORIGINAL_FIXTURE NEW_FIXTURE
build/test/ckks/ckks_ram_fixture_check NEW_FIXTURE --ram-noise
python3 test/ckks/tools/run_ram_ring_switch.py \
  --binary build/test/ckks/ckks_narrow_ll14_retuned_ram \
  --fixture NEW_FIXTURE --phase sigma256 --output RUN_ROOT/results
```

The destination must be new. This private-owner import is specifically for the
original legacy fixture; it does not reuse public evaluation material. It
creates all five traces and both output bridges. `--ram-noise` checks one full
RAM write against ideal sine restoration before the expensive CKKS trial.
The qualification run checks capture, sixteen accesses per backend on small
RAM, and five alternating accesses per backend on the 8192-cell RAM, including
one warmup each. Every access checks every cell. The analyzer also requires
synthetic NAND, circuit-bootstrap and encrypted-feedback checks and measured
initial-memory/private-switch noise near sigma256.

This is a correctness and timing experiment under assumed circular security.
The existing classical screen gives level 1 minima of 135.19/133.09 bits under
BDGL16/MATZOV; the active level-0 minimum remains 128.51 bits. The new private
switch still supplies 49176 independent known-zero samples, exceeding the
1024 samples used by the screened optimum. Changing decomposition therefore
does not remove that exposure. The retained diagnostic owners are not a fresh
security deployment: old disclosed weaker material remains useful against the
same owners. Deployment requires fresh connected secrets and regenerated
material. Exact-sampler, quantum and composition qualifications remain open.


The [sigma256 qualification receipt](ckks-ram-sigma256-validation.json) records
the successful run and the failed noise-only attempt. All ten full-size
accesses passed, checking 81920 cell bits; four accesses per backend were
measured after warmup. Mean access times were 7.377 seconds with CKKS/ring
refresh and 81.851 seconds with serial TFHE (11.096x observed speedup).
The worst full-size pre-switch error was 0.042944, below 0.0625. Measured CKKS
refresh stages averaged 0.2684 seconds packing, 4.0057 seconds bootstrap,
0.5913 seconds decoding and 0.01185 seconds output switching. All 22 original
fixture files retained their SHA256 hashes. Initial-memory and private-switch
integer noise RMS were 252.43 and 256.76, respectively. The synthetic gate
checks also passed (64 NAND checks, 8192 circuit-bootstrap output bits and
16 encrypted-feedback accesses).

The host-wide external-load monitor peaked at 0.255 CPU cores and failed the
strict 0.2-core isolation threshold. These timings are observations from this
paired run. The receipt does not claim uncontended performance or complete
security validation. Default sigma128 parameters remain unchanged.
