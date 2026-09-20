# Scalar TFHE → BFV multiplication → TFHE

`scheme_switch_comparison bfv` evaluates the same 256 deterministic unsigned
8-bit operand pairs as the existing `clpx` and `ckks` modes, and verifies every
16-bit product. BFV and CLPX both consume eight Boolean TFHE ciphertexts per
operand and return sixteen Boolean ciphertexts. CKKS retains its existing
single integer-phase TLWE interface. These are comparable arithmetic workloads,
but not identical TFHE interfaces.

The benchmark also has a [separately tuned 16-bit BFV profile](#tuned-16-bit-bfv-profile).

## Measured comparison

On an Intel Core i9-9900X, with eight OpenMP threads and the 24-GiB cap,
all three schemes passed all 256 products. Times below are seconds per product;
forward includes both operands. Key generation, encryption, decryption, and
diagnostic measurements are excluded.

| Scheme | Forward | Multiply | Reverse | Total | Peak RAM (GiB) |
|---|---:|---:|---:|---:|---:|
| BFV | 2.280294 | 0.029605 | 0.242755 | **2.552654** | 3.489 |
| CLPX | 0.120619 | 0.000129 | 3.117131 | 3.237879 | 1.076 |
| CKKS | 2.497689 | 0.674486 | 0.426352 | 3.598527 | 14.050 |

BFV took 21.2% less evaluation time than CLPX and 29.1% less than CKKS in this
run. BFV's forward conversion dominates its cost; CLPX spends most of its time
in reverse conversion. Key generation took 21.53 s, 5.58 s, and 18.59 s for
BFV, CLPX, and CKKS respectively. This is one same-host run per scheme, not a
confidence interval or a claim about other hardware.

The maximum measured BFV coefficient error was 4.72×10^-17 plaintext units
after forward conversion and 2.21×10^-10 after multiplication. Both are below
the conditional envelopes. The whole-run conditional bound is approximately
2^-84.99, subject to the assumptions below; the worst-case product envelope is
about 0.00196 plaintext units. Runtime parameter values and key storage were
checked against the analysis; source hashes identify the implementation used
for that original 8-bit run.

Evidence: [machine-readable comparison](bfv-switch-comparison-20260919/comparison.json),
[BFV log](bfv-switch-comparison-20260919/bfv.log),
[CLPX log](bfv-switch-comparison-20260919/clpx.log),
[CKKS log](bfv-switch-comparison-20260919/ckks.log),
[encrypted correctness tests](bfv-switch-comparison-20260919/correctness.txt),
[error-analysis tests](bfv-switch-comparison-20260919/error-analysis-tests.txt).

## Build and run

From the TFHEpp directory:

```sh
cmake -S . -B build-bfv-switch \
  -DENABLE_TEST=ON -DUSE_BLOCK_BINARY=OFF \
  -DUSE_SPQLIOS_INTL=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build-bfv-switch --target \
  bfv_scheme_switch cmux_double_decomposition scheme_switch_comparison simd_ops -j 4

OMP_NUM_THREADS=8 OMP_STACKSIZE=32M MALLOC_ARENA_MAX=2 \
  build-bfv-switch/test/bfv/bfv_scheme_switch --full

python3 benchmark/run_scheme_switch_comparison.py \
  --executable build-bfv-switch/test/clpx/scheme_switch_comparison \
  --output docs/my-bfv-comparison
```

The runner requires a new output directory. It launches BFV, CLPX, and CKKS in
separate sequential processes with eight OpenMP threads and a 24-GiB address
space limit. It saves raw logs, `/usr/bin/time` resource measurements, source
and executable hashes, compiler flags, and `comparison.json`. Encryption,
decryption, and diagnostic phase measurements are excluded from evaluation
timers. Key generation has a separate timer. Forward timing includes both
operands, multiplication includes relinearization, and reverse timing includes
all extraction and Boolean bootstraps. Per-product times are aggregate time
divided by the product count; they are not latency percentiles.

### Experimental width sweep

The benchmark accepts `scheme_switch_comparison SCHEME OPERAND_BITS PRODUCTS`.
The original one-argument command still selects 8-bit operands and 256 products.
CLPX accepts widths through 64, with a full 128-bit reference product at width
64. The scalar BFV benchmark instantiates widths 8, 9, 10, 12, 16, 24, and 32;
only width 16 selects the separately tuned larger-ring profile below. The
other widened profiles remain experimental, not a claim that all have enough precision.
CKKS's existing scalar benchmark remains restricted to 8-bit operands.

```sh
python3 benchmark/run_scheme_switch_comparison.py \
  --executable build-bfv-switch/test/clpx/scheme_switch_comparison \
  --output docs/my-width-16 --schemes bfv clpx --bits 16 --products 16
```

Inputs start with eight boundary pairs, including maximum times maximum and
values around the top bit, followed by deterministic pseudorandom pairs. The
reference product and decrypted output use unsigned 128-bit integers, avoiding
host-integer overflow. A scheme must pass every product to receive a successful
timing result. The runner preserves failures and continues to the next scheme.

`BFVSwitchParamFor<W>` and `TFHEBFVSwitchKeyFor<W>` expose the experimental
scalar profile. The original names remain aliases for `W=8`. The widened
plaintext modulus is `2^(2W)` and the scale is `2^(128-2W)`. Except for the
new `W=16` specialization, the original ring/noise/decomposition parameters
are retained. The reverse decomposition extracts all `2W` bits plus a guard
digit. The historical sweep below used the original parameters even at `W=16`.

The companion analysis accepts `python3 python/BFVswitch.py --operand-bits W`
from Parameter-Selection. It scales the forward noise, multiplication envelope,
reverse decomposition, and union-event counts with width. Under its stated
assumptions, the fixed profile still meets the whole-256-product `2^-40` target
at width 9 (bound approximately `2^-61.37`), but not at width 10 or above.
Use `--legacy-profile` to reproduce that old 16-bit analysis; without it,
`--operand-bits 16` now selects the tuned profile.
An empirical pass at those larger widths does **not** establish the requested
failure guarantee. A wider limb-based BFV implementation would be a different
algorithm and would need carry propagation included in both timing and analysis.

The sweep exposed and fixed a CLPX partial-block bug: the last block's weights
and guard bias used its shortened length, while HomDecomp still used the full
16-bit data window. Zero decoded incorrectly at output widths 18 and 24.
Weights and bias now retain the full window scale; only the number of processed
coefficients shrinks. Full 16-bit blocks are unchanged. The CLPX regression
test checks zero/max values at widths 2, 18, and 24, plus the existing 16-bit
cases, with both supported reverse decompositions (basebit 2 and 4).

#### Legacy width-sweep measurements (2026-09-20, before 16-bit tuning)

Times are seconds per full product, excluding setup and diagnostic decryption.
Forward includes both operands. Runs use the same eight-thread/24-GiB settings
as the original comparison. Widths 8–16 request 16 products (eight boundaries
and eight random pairs); widths 32 and 64 request eight boundary products.

| Operand bits | Full product bits | BFV | CLPX | BFV conditional target, 256 products |
|---|---:|---:|---:|---|
| 8 | 16 | 2.592957 | 3.233255 | Met |
| 9 | 18 | 4.504278 | 3.992835 | Met |
| 12 | 24 | 4.768454 | 5.148321 | Not met |
| 16 | 32 | Correctness failure | 6.464906 | Not met |
| 32 | 64 | Correctness failure | 13.200426 | Not met |
| 64 | 128 | Not implemented | 26.263508 | Not modeled |

Every CLPX run in this table passed all requested products. Successful BFV
rows also passed all requested products; failed BFV runs stop at the first
mismatch and do not receive an evaluation timing. The conditional target
column concerns BFV only, with the assumptions described below; it is not an
empirical error rate or a new wide-CLPX failure guarantee.

The eight-thread screening run found a crossover at **9-bit operands**:
CLPX took 3.992835 s/product versus BFV's 4.504278 s/product, approximately
11.4% less time. Both passed all 16 products. This is implementation- and
thread-count-specific: nine forward BFV bootstraps per operand require a second
wave on eight workers, increasing forward time from 2.315756 s at width 8 to
4.186142 s at width 9. It is **not a monotonic crossover**: at width 12, BFV
took 4.768454 s versus CLPX's 5.148321 s, both passing all 16 products, but
without a satisfactory BFV conditional failure bound at that width.

The 16-bit scalar BFV run failed on `57380 * 57482`, returning `3298317159`
instead of `3298317160`; its maximum observed product error was approximately
0.951189 plaintext units. At 32 bits, BFV already failed on `0 * 4294967295`.
A failed run is not a valid speed comparison. A fair BFV comparison at these
larger widths needs a new parameter set or limb arithmetic, not just a larger
plaintext modulus. CKKS is not part of this wider Boolean-interface sweep.

Raw evidence and source/executable fingerprints are in
[the width-sweep directory](bfv-width-sweep-20260920/), with
[the successful 9-bit comparison](bfv-width-sweep-20260920/9/comparison.json)
and [its checked conditional bound](bfv-width-sweep-20260920/bounds/9-measured.json).
The `pre-partial-block-fix` subdirectory preserves the initial CLPX failures;
those runs must not be used as successful CLPX timings. The width-8 pilot was
run before that fix, which leaves full 16-bit blocks unchanged. These are short
screening runs, not confidence intervals or empirical failure-rate estimates.

Verification passed: [CLPX encrypted partial-block regressions](bfv-width-sweep-20260920/clpx-partial-block-tests.txt),
[the original BFV encrypted regression](bfv-width-sweep-20260920/bfv-regression-tests.txt),
[seven BFV analysis tests](bfv-width-sweep-20260920/error-analysis-tests.txt),
and [nine CLPX analysis tests](bfv-width-sweep-20260920/clpx-analysis-tests.txt).
Source hashes identify the pre-tuning implementation used for these runs.

The default test checks isolated PBS and packing, forward boundary values,
reverse digit boundaries with signed injected error, and six complete products.
`--full` additionally checks all 256 forward values and all 4096 coefficients of
each resulting plaintext polynomial. Failures use runtime checks, including in
Release builds. `cmux_double_decomposition` tests encrypted selectors 0 and 1
at rotations 0, 1, N, and N+1.

## Tuned 16-bit BFV profile

`BFVSwitchParamFor<16>` and `TFHEBFVSwitchKeyFor<16>` now use a separate
parameter set. The 8-bit profile, Boolean interfaces, scalar constant-polynomial
encoding, and full 32-bit product for 16-bit operands are unchanged in meaning.
No limb arithmetic, truncation of the product, or secret-key evaluation is used.

| Parameter | 8-bit operands | 16-bit operands |
|---|---:|---:|
| Ring dimension | 4096 | 8192 |
| Ciphertext modulus | 2^128 | 2^128 |
| Plaintext modulus | 2^16 | 2^32 |
| BFV scale | 2^112 | 2^96 |
| Fresh relative noise | 2^-105 | 2^-124 |
| Primary gadget (levels × bits) | 4 × 21 | 12 × 10 |
| Auxiliary DD (levels × bits) | 8 × 16 | 8 × 16 |
| BFV→TFHE key switch (levels × bits) | 7 × 2 | 6 × 3 |
| Stored evaluation keys (GiB) | 3.399 | 18.913 |
| Conditional log2 failure, 256 products (original bitwise conversion) | -84.99 | -54.30 |

The original bitwise-forward [tuned comparison](bfv16-tuned-20260920/comparison/comparison.json) passed
all 16 full 32-bit products for both BFV and CLPX (eight boundary and eight
seeded-random pairs), including the previously failing `57380 * 57482` pair.
On the same eight-thread host, BFV took **24.852042 s/product**, versus
**6.452809 s/product** for CLPX: CLPX was **3.85× faster** in this run with
both schemes correct. This is a short benchmark, not a failure-rate estimate.

BFV forward/multiply/reverse costs were 24.299024 / 0.058802 / 0.494216 seconds
per product. Setup was separate at 236.914 seconds, and measured peak RSS was
19.054 GiB, below the 24-GiB address-space cap. Maximum observed forward and
product errors were 6.42e-21 and 2.96e-9 plaintext units; both are below the
[compiled-parameter-checked conditional envelopes](bfv16-tuned-20260920/bound-measured-16.json).
The 256-product bound is analytical; this measured run tested 16 products,
not 256. The 8-bit profile is unchanged.

The original bitwise-forward 16-bit search tests 636 memory-feasible candidates in the separately
screened `N=8192`, integer-unit noise `sigma=16` instance. It requires eight
extra failure-exponent bits beyond the requested `2^-40`, ranks by primary
gadget rows then stored key bytes, and reserves four GiB within the 24-GiB cap
for buffers/stacks. Three candidates pass this screen; the selected profile
has the smallest key size among those with the fewest rows. This is a bounded
cost-model search, not a globally optimal measured-runtime claim.

```sh
# From Parameter-Selection, reproduce parameter selection and its assumptions:
python3 python/BFVswitch.py --operand-bits 16 --tune-16bit \
  --products 256 --failure-bits 40 --output /tmp/bfv16-bound.json
sage -python python/estimates/bfv_switch_security.py \
  --operand-bits 16 --output /tmp/bfv16-security.json

# From TFHEpp, after building with USE_SPQLIOS_INTL=ON:
python3 benchmark/run_scheme_switch_comparison.py \
  --executable build-bfv-switch/test/clpx/scheme_switch_comparison \
  --output docs/my-tuned-16 --schemes bfv clpx --bits 16 --products 16
```

The [conditional report](bfv16-tuned-20260920/bound-256.json) uses the same
explicit sub-Gaussian/FFT assumptions as the original analysis, with 13 trace
stages and the new ring and key-switch parameters. It is not an unconditional
failure proof. The [four-attack LWE screen](bfv16-tuned-20260920/security.json)
is a heuristic and still excludes ring-structure and circular/KDM attacks.
The implementation truncates a continuous Gaussian to integer noise; the
security report therefore also screens a smaller `sigma=8` as a sensitivity
check, not as a proof that the sampler equals a discrete Gaussian.
The weakest ring estimates were 236.35 bits at sigma=16 and 234.29 bits at
sigma=8; the unchanged TFHE bridge remains the weakest screened instance
in the whole chain, at 129.46 bits.

Implementation support includes a thread-local 8192-point SPQLIOS digit FFT
processor and heap-backed ordinary TRGSW key-generation buffers, avoiding the
default main-stack limit. [Integer FFT and encrypted CMUX tests](bfv16-tuned-20260920/fft-cmux-tests.txt)
exercise both ring sizes; [analysis tests](bfv16-tuned-20260920/analysis-tests.txt)
verify the selected profile, memory bound, and signed product-error/carry
boundaries through 32-bit output values.
The [48 encrypted reverse-switch boundary/error cases](bfv16-tuned-20260920/reverse-tests.txt)
passed with signed injections of 1/16 plaintext unit, including all-one output
bits and the previously failing product. The [8-bit encrypted regression](bfv16-tuned-20260920/bfv8-regression.txt)
also passed after the shared FFT/key-generation changes. Reproduce those with
the `bfv16_reverse`, `cmux_double_decomposition`, and `bfv_scheme_switch` targets.

### Paired forward conversion

The current 16-bit forward path first uses small Boolean bootstraps to
re-encode each sign at `+/-Q/16`, then combines adjacent bits with weights
1 and 2. The four digit phases are `(-3,-1,1,3)/16`. A triangular,
anti-periodic large-ring LUT plus a public offset converts each pair into
its weighted radix-4 BFV digit. This reduces large-ring bootstraps from
32 to 16 per product, at the cost of 32 small Boolean bootstraps. It reuses
the existing keys; all parameters, the reverse path, and the 8-bit path
remain unchanged. Directly dividing a torus ciphertext by two is not used.

For the conditional error model, each packed address must satisfy
`3 B1 + K1h + A3_address < 1/16`. Initial normalization addresses have
radius `1/8`. Each operand accumulates eight rather than sixteen large-PBS
errors, followed by the same thirteen packing stages. Event counts include
the extra normalization bootstraps and packed key switches, without assuming
independence. `BFVswitch.py --operand-bits 16` models this circuit;
`--bitwise-forward` reproduces the previous tuned circuit, and
`--tune-16bit` retains the original bitwise parameter-selection search.
Benchmark validation checks `bfv_forward_digit_bits` as well as the
cryptographic parameters, rejecting a mismatched circuit report.

The [new 256-product conditional bound](bfv16-paired-20260920/bound-256.json)
uses the same explicit primitive-tail and FFT assumptions. The prior
[security screen](bfv16-tuned-20260920/security.json) still describes the
unchanged LWE instances; its source hash belongs to the earlier implementation,
not this optimization. It is not a new security proof for the circuit.

The [paired-forward benchmark](bfv16-paired-20260920/comparison/comparison.json)
passed all 16 full 32-bit products on the same eight-thread host and 24-GiB
address-space cap. BFV evaluation fell from **24.852042 to 12.628134 seconds
per product**, a **1.968× speedup** (49.2% less time). Forward conversion
fell from 24.299024 to 12.076074 seconds; multiplication and reverse conversion
took 0.058242 and 0.493819 seconds. Setup remained separate at 237.295 seconds.
This is a comparison of single short runs, not a statistical speed estimate.
The [measured-run bound check](bfv16-paired-20260920/bound-measured-16.json)
matches compiled parameters, digit width and key size, and verifies observed
errors lie within the envelopes. The analytical 256-product bound is
approximately `2^-83.91`, not an empirical failure-rate estimate.
Peak RSS remained **19.054 GiB**. The fresh CLPX run also passed all 16 cases
at **6.550864 seconds/product**, so CLPX remains **1.93× faster** for this
scalar 16-bit comparison. These timings exclude setup for both schemes.

Regression evidence: [analysis and radix-4 LUT tests](bfv16-paired-20260920/analysis-tests.txt),
[8-bit encrypted round trips](bfv16-paired-20260920/bfv8-regression.txt),
[16-bit reverse boundary/error cases](bfv16-paired-20260920/reverse-tests.txt),
and [integer FFT/CMUX checks](bfv16-paired-20260920/fft-cmux-tests.txt).

```sh
# Workspace root, after running the comparison above into a new directory:
python3 Parameter-Selection/python/BFVswitch.py --operand-bits 16 --products 16 \
  --benchmark-report TFHEpp/docs/my-tuned-16/comparison.json
python3 Parameter-Selection/python/BFVswitch.py --operand-bits 16 --products 256
```

## API and encoding

The public header is `include/bfv/scheme_switching.hpp`, also included by
`tfhe++.hpp`. It is enabled for the default 128-bit parameter family.

```cpp
using namespace TFHEpp;
Key<lvl1param> io_key;
Key<lvlhalfparam> half_key;
Key<BFVSwitchParam> bfv_key;
keyGen<lvl1param>(io_key);
keyGen<lvlhalfparam>(half_key);
keyGen<BFVSwitchParam>(bfv_key);
TFHEBFVSwitchKey evaluation_key;
TFHEBFVSwitchKeyGen(evaluation_key, io_key, half_key, bfv_key);

// a_bits and b_bits: BFVSwitchInput, encrypted under io_key, LSB first.
TRLWE<BFVSwitchParam> a, b, product;
TFHEToBFV(a, a_bits, evaluation_key);
TFHEToBFV(b, b_bits, evaluation_key);
TRLWEMultFullDD<BFVSwitchParam>(product, a, b, *evaluation_key.relin);
BFVSwitchOutput result;
BFVToTFHE(result, product, evaluation_key);
```

The evaluator stores no secret key. The explicit three-key generation API
avoids accidentally mixing the library's independent and subset key chains.
Inputs and outputs use standard `lvl1param` Boolean phases ±μ.

BFV uses N=4096, Q=2^128, t=2^16, and exact Δ=2^112. Ring keys are uniform
ternary, with α=2^-105; the primary gadget is 4×21 bits and the auxiliary
decomposition is 8×16 bits. This is a scalar constant-polynomial encoding,
not the prime-modulus SIMD encoding. No modular wrap occurs for 8-bit products.

Forward conversion performs eight weighted Boolean PBS operations per operand.
For bit i the LUT amplitude is 2^i·Δ/2; adding the same public offset produces
0 or 2^i·Δ without amplifying bootstrap error by the bit weight. The summed
TLWE is inverse-sample-extracted and packed with the iterative annihilation
trace. This avoids the recursive general packer's large 128-bit stack frames.

Reverse conversion extracts coefficient zero, adds Δ/8, and applies
`HomDecomp<..., 2, 9, 18>`. The window contains two low guard bits and sixteen
data bits. The guard digit is discarded. Each base-4 data digit is centered
by Q32/8, switched to the half-level key, and converted into two Boolean bits.
Directly switching the original 16-bit phase into the small bootstrap domain
would not preserve enough precision.

The implementation also fixes the DD branch of
`CMUXFFTwithPolynomialMulByXaiMinusOne`. The fused single-decomposition path
indexed DD keys incorrectly. DD now uses `ExternalProduct` with all auxiliary
rows and limb recombination. The single-decomposition path is unchanged.

## Conditional error envelope

The report [bfv-switch-bound.json](bfv-switch-bound.json) is generated by
`Parameter-Selection/python/BFVswitch.py`. It targets a union bound of 2^-40
over the entire 256-product run. The default parameters meet that target under
the report's assumptions; no gadget tuning is needed.

This is **not an unconditional failure-probability proof**. Primitive variance
formulas are used as assumed sub-Gaussian proxy scales, uniformly over the
circuit inputs. DD FFT digit convolutions, including those used during key
generation, are assumed to round to the exact integer result. Non-DD Boolean
PBS has a separate deterministic numerical-error allowance of 2^-20 in
normalized torus units. Measured agreement does not prove these assumptions.

Composition uses absolute error envelopes and a union bound, not independence
of ciphertext errors or evaluations with reused keys. For a uniform cutoff z,
each primitive event is assumed to have tail at most `2 exp(-z²/2)`. The report
counts every coefficient of every packing automorphism and relinearization,
as well as the scalar PBS, key-switch, address-rounding, and input events.
It evaluates at 99% of the largest cutoff satisfying all strict inequalities.

All quantities below are normalized by the appropriate ciphertext modulus.
Let B3 bound a forward PBS error, A3 an automorphism error, and R3 a
relinearization error. Each operand has coefficient error at most

```
E = 8 B3 + 12 (A3 + (N+1)/Q).
```

For BFV multiplication, ciphertext coefficient lifts satisfy
`v = Δm + e + Qk`. Balanced DD digits represent coefficients of magnitude less
than Q, so `|k| < N+2` for the supported inputs. Retaining all lift/error cross
terms gives the conservative product envelope

```
Eproduct = 510 E + 2 t N (N+2) E + N t E² + (1+N+N²)/Q + R3.
```

In particular, `t k e` terms must not be dropped. The final rounding term
bounds componentwise tensor rounding against 1, s, and s². The scalar
messages have magnitude at most 255; other plaintext coefficients are zero.

Let K31 and K1h bound the reverse identity switches, B1 a Boolean PBS output
error including its numerical allowance, and A1 a Boolean address-rounding
error. The initial guard correction requires

```
2^16 Eproduct + K31 + K1h + A1 < 1/4.
```

For data digit j=1..8, set
`D_j = 2^(16-2j) Eproduct + K31 + B1 + K1h`.
The correction and high-bit LUTs require `D_j + A1 < 1/16`; the low-bit LUT
requires `2 D_j + A1 < 1/8`. Finally B1 must be below 1/8 for Boolean
decryption. The exact scalar reference checks all 65,536 outputs; tests also
inject signed product errors and per-primitive envelope errors at carry edges.

From the workspace root:

```sh
python3 -m unittest discover -s Parameter-Selection/python -p test_bfv_switch.py -v
python3 Parameter-Selection/python/BFVswitch.py \
  --benchmark-report TFHEpp/docs/my-bfv-comparison/comparison.json \
  --output TFHEpp/docs/my-bfv-bound.json
```

Supplying the benchmark report checks the modeled parameters and evaluation-key
size against compiled values. `--search` retains the default if it passes;
otherwise it searches both identity-switch gadgets over 1–10 levels and
1–8 bits per level, retaining at most 32 bits. It minimizes key storage and
then total levels, reserving 4 GiB of the 24-GiB budget for runtime buffers.
Search results are analytical candidates; a changed candidate must be compiled
and measured before use in a performance comparison.

## Security screen

[bfv-switch-security.json](bfv-switch-security.json) and its
[raw log](bfv-switch-security.log) record fresh four-attack lattice-estimator
results using RC.BDGL16 and unbounded samples. The weakest estimated costs are
129.65 bits for the ternary TFHE I/O key, 129.46 bits for the binary half-level
key, and 130.21 bits for the ternary BFV ring key. All exceed 128 bits.
The recorded run used SageMath 10.6 and estimator revision
`3e48ef421ec256afddb3e7d2249a77eab6e9ba12`.

This is a generic LWE/RLWE-as-LWE heuristic screen. It does not certify ring
structure, circular/KDM security, or the conditional error assumptions.
Reproduce with SageMath and the checked-in estimator:

```sh
sage -python Parameter-Selection/python/estimates/bfv_switch_security.py \
  --output TFHEpp/docs/my-bfv-security.json
```

Scope is one scalar multiplication with Boolean TFHE I/O. SIMD batching,
larger integers, and BFV bootstrapping are not part of this implementation.
