# Scalar TFHE → BFV multiplication → TFHE

`scheme_switch_comparison bfv` evaluates the same 256 deterministic unsigned
8-bit operand pairs as the existing `clpx` and `ckks` modes, and verifies every
16-bit product. BFV and CLPX both consume eight Boolean TFHE ciphertexts per
operand and return sixteen Boolean ciphertexts. CKKS retains its existing
single integer-phase TLWE interface. These are comparable arithmetic workloads,
but not identical TFHE interfaces.

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
checked against the analysis, and all recorded source hashes match the final
implementation.

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
divided by 256; they are not latency percentiles.

The default test checks isolated PBS and packing, forward boundary values,
reverse digit boundaries with signed injected error, and six complete products.
`--full` additionally checks all 256 forward values and all 4096 coefficients of
each resulting plaintext polynomial. Failures use runtime checks, including in
Release builds. `cmux_double_decomposition` tests encrypted selectors 0 and 1
at rotations 0, 1, N, and N+1.

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
