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

## Measured comparison

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
