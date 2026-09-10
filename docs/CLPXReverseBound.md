# End-to-end error-envelope bound for the 16-bit CLPX reverse circuit

The current measured reverse circuit is **not certified at 2^-40 by this bound**.
The new analysis covers the whole reverse computation, including its two ramp
bootstraps, digit rounding, accumulator, HomDecomp, output-bit bootstraps, and
final decryption. It is a deterministic correctness implication plus a
**conditional probability bound**. It does not turn Parameter-Selection variance
estimates into proved tail bounds or certify floating-point FFT arithmetic.

## Results for 256 reverse conversions

CLPX radix remains 2 in every row. `basebit` specifies the HomDecomp digit width,
so `basebit=4` means radix-16 decomposition and `basebit=2` means radix-4.
All rows return 16 Boolean TFHE ciphertexts from one CLPX ciphertext.

| Reverse settings | Largest uniform cutoff | Conditional whole-run bound | Clears 2^-40 |
|---|---:|---:|---|
| Measured: basebit=4, numdigit=5, half n=760, alpha=2^-17 | 3.869 sigma | Vacuous (<=1) | No |
| basebit=2, numdigit=9; same cryptographic parameters | 7.716 sigma | <=2^-25.17 | No |
| basebit=2, numdigit=9; half n=832, alpha=2^-18 | 8.997 sigma | <=2^-40.32 | Conditionally |
| basebit=2, numdigit=9; half n=832, alpha=2^-19 | 9.652 sigma | <=2^-48.95 | Conditionally |

The evaluated cutoff is 99% of the maximum, reserving strict error-radius slack.
The last candidate's intermediate LWE instance passes the four-attack classical
RC.BDGL16 security screen at **128.4086 bits**; the other chain instances already
screen above 128 bits. This candidate has **not been implemented or benchmarked
with encrypted ciphertexts**. Its performance cannot be inferred from the
previous measured CLPX row. The original measured configuration remains intact.

Reports: [measured](clpx-reverse-bound.json),
[smaller decomposition digits](clpx-reverse-base2-bound.json),
[alpha=2^-18 candidate](clpx-reverse-candidate-bound.json),
[alpha=2^-19 candidate](clpx-reverse-low-noise-bound.json), and
[candidate security](clpx-reverse-low-noise-security.json).

Failure of a sufficient bound does **not** prove that the actual failure
probability exceeds the target. In particular, the measured row's bound of 1
is not a prediction that the circuit fails routinely.

## Fixed-dimension alternatives found subsequently

The intermediate TLWE dimension does **not** need to change. A subsequent
search retains `n=760`, `alpha=2^-17`, the existing ring parameters, and
`basebit=2, numdigit=9` HomDecomp, while tuning only the level-2-to-half identity
key-switch decomposition. Security instances are unchanged (the existing
unbounded-sample screen has a minimum of 129.46 classical bits).

| KS20 decomposition | Conditional whole-run bound | KS20 key storage alone |
|---|---:|---:|
| Existing: 7 levels, 2 bits/level | <=2^-25.17 | 0.0813 GiB |
| 3 levels, 5 bits/level | **<=2^-43.47** | **0.2787 GiB** |
| 2 levels, 8 bits/level | **<=2^-53.21** | **1.4863 GiB** |

The last row improves on the earlier 832-dimensional candidate's <=2^-48.95
conditional bound. The middle row is the smallest KS20 key meeting the original
2^-40 target in the searched family (1..10 levels, 1..12 bits/level, at most
32 retained bits and a 3 GiB KS20-key cap). This is not a global optimum.

The identity switch selects one encrypted signed digit per level; it does not
multiply a single noisy encryption by the digit value. Its key-noise term is
therefore proportional to the number of levels. Fewer levels reduce noise,
while the wider digits keep truncation error small. The price is exponentially
more stored digit encryptions. Storage above follows the current balanced key
layout: `2048 * levels * 2^(basebit-1) * 761 * 4` bytes, excluding other keys and
process overhead. Blind-rotation rounding itself is unchanged.

These remain analytical candidates under the same unproved primitive-tail and
FFT assumptions. Four scalar-trace tests pass, including error injections for
both new candidates. The t=3, basebit=5 candidate subsequently passed
[the encrypted 256-product performance benchmark](CLPXT3Performance.md); the
t=2, basebit=8 candidate has not been encrypted-benchmarked.
Reports: [smaller key](clpx-reverse-fixed760-small-key-bound.json),
[stronger bound](clpx-reverse-fixed760-bound.json), and
[search results](clpx-reverse-fixed760-search.json).

From Parameter-Selection, reproduce with:

```sh
python3 python/CLPXreverse.py --basebit 2 --ks20-levels 3 --ks20-basebit 5 --output /tmp/fixed760-small.json
python3 python/CLPXreverse.py --basebit 2 --ks20-levels 2 --ks20-basebit 8 --output /tmp/fixed760-strong.json
```

## Input contract and scope

All errors and radii below use normalized torus units. For the intended integer
P modulo 2^16, the 16 input coefficient phases must lie close to

    u_i = -P / 2^(i+2)  (mod 1/2),  i=0,...,15.

This is the radix-2 CLPX phase relation. Let E_in bound each coefficient's error
around these valid residues. The probability that this input contract fails
must be included. The numerical adapter uses the previously calculated CLPX
product coefficient variance 2^109.20114828158572 in 64-bit integer units as an
*assumed sub-Gaussian tail scale*, not as a proof about the product distribution.

The scope is one 16-bit block, repeated at most 256 times with reused keys.
It does not cover inter-block carry bootstraps for wider outputs. The existing
basebit=2 implementation uses two carry digits and basebit=4 uses one, but no
inter-block carry bootstrap executes for this output width.
For an entire TFHE->CLPX multiplication->TFHE proof, add upstream semantic
failure probabilities to the reverse bound; a favorable multiplication variance
alone is insufficient.

## Deterministic correctness implication

Name the primitive absolute error envelopes:

- K20, K21, K1h: the three identity key-switch error envelopes, after normalizing
  by their output modulus; they exclude propagated input error.
- R2, R1: blind-rotation address approximation errors for the two ring sizes.
- B2, B1: bootstrap output errors around the actual lookup-table value. Any FFT
  numerical error must also fit these envelopes.

No independence assumption is needed for the following inequalities.

### Ramp composition and digit recovery

For normalized x in [0,1), the first ideal ramp is

    F(x) = x                  if x < 1/2,
           1/2 - x            otherwise.

The C++ second table scales this ramp by 1/8. Its essential identity is

    F(F(x))/8 = x/8  (mod 1/16).

Projecting F onto the circle modulo 1/2 is 1-Lipschitz, including at its branch
boundaries. Thus the ramp pair need not choose a unique real lift and does not
need a separate exact-address success event. Its output error, measured modulo
1/16, is at most

    E_F = (E_in + 2 K20 + 2 R2 + B2)/8 + B2.

For i>0, the digit-rounding input is `previous_scaled - 2*fid_stage2 + 1/32`.
Its sufficient condition is

    3 E_F + K20 + R2 < 1/32.

For the first coefficient, replace 3 by 2 because `previous_scaled` is zero.
We also require `3(E_F+B2) < 1/32` to retain the intended integer-digit range of
the triangular rounding table. Arbitrary representative changes by multiples
of 1/16 yield integer digit changes `k_(i-1)-2*k_i`. They telescope in the
weighted sum modulo 2^16. Therefore an inconsistent *real* lift of a ramp
output is not automatically a failure.

Once these rounding decisions succeed, the 16 outputs accumulate as an integer
Z congruent to P modulo 2^16, encoded as Z/2^20, with error bounded by

    A = 16 B2.

The C++ code adds the exact bias 2^-22 before HomDecomp.

### HomDecomp and final Boolean extraction

For basebit b, let d=16/b+1 and c be the number of carry digits. The actual
source window is W=b(d+c): **24 bits for b=4, 22 bits for b=2**, not 64 bits.
For zero-based digit j, the error in the returned `cres[j]` is bounded by

    C_j = 2^(W-b(j+1)) A + K21 + [j>0] B1.

A subtract-PBS decision that feeds requested output bits is correct if

    C_j + K1h + R1 < h_j.

For bit k within returned digit j>=1, the final bit-PBS condition is

    2^(b-k-1) (C_j + K1h) + R1 < o_(j,k).

Finally, B1 < 1/8 guarantees decryption of the Boolean output encoded at +/-1/8.
A wrong internal decision is included in the failure union even if later
operations could happen to correct it, making the implication conservative.

The h and o radii are computed by an exact integer oracle following the C++
shifts, offsets, subtract-PBS signs, and final bit ordering. It verifies all
65,536 residues and records witnesses attaining each minimum. A multiple of
2^16 added to Z changes only discarded high digits, so these radii also cover
representative-dependent carries within this one output block.

For the measured basebit=4 path, subtract-PBS radii feeding requested output
bits are 7/32, 1/32, 17/1024, and 257/16384. The final output bit has radius
**4097/262144 = 0.015628814697265625**, attained at P=0. Its envelope coefficient
is approximately 0.004039846589 per sigma. At the target cutoff 8.871 sigma,
the error envelope is 0.0358381, which exceeds this radius. This is the measured
configuration's bottleneck. Using the final TLWE's fresh noise alone misses it.

## Probability composition and assumptions

For every primitive event e, suppose an error X_e obeys

    Pr[ |X_e| > t*s_e + d_e ] <= 2 exp(-t^2/2).

The s_e and d_e are supplied tail scales and deterministic bias allowances.
Substitute these envelopes into the deterministic conditions above. If all
conditions hold, then, for K total primitive events,

    Pr[any output is incorrect] <= min(1, 2 K exp(-t^2/2)).

This is a union bound, not a product of success probabilities. Calls may reuse
evaluation keys, inputs may be correlated, and all primitive errors may be
correlated. Equivalently, separately proved primitive probabilities can be
summed without using this uniform Gaussian-shaped envelope at all.

For the measured configuration, the per-product counts are:

| Event | Count |
|---|---:|
| Input coefficient error | 16 |
| K20 | 48 |
| R2 | 48 |
| B2 | 48 |
| K21 | 6 |
| K1h | 10 |
| R1 | 21 |
| B1 | 21 |

That is 218 events per product and K=55,808 for 256 products. Counts conservatively
include unused executed HomDecomp work. With basebit=2, K=61,952.
For the stronger candidate, t=9.555576496 yields log2(bound)=-48.946687933.

**What is proved and what remains assumed:** the circuit's sufficient
conditions, finite noiseless semantics, and union-bound composition are explicit.
Parameter-Selection currently supplies variance estimates; these do not prove
sub-Gaussian tails for adaptive key-switch/bootstrap errors, and a product of
noisy ciphertexts need not have Gaussian noise. The adapter also lacks a
certified FFT-error envelope. Consequently all reports deliberately set
`implementation_certified=false`, even for candidates whose conditional bound
passes. This is not yet an unconditional 2^-40 certificate for the C++ program.

## Implementation and validation

New Parameter-Selection files:

- `python/noiseestimation/clpx_reverse_bound.py`: exact margin oracle,
  deterministic envelope composition, event counting, and probability adapter.
- `python/CLPXreverse.py`: parameterized CLI and source-hash recording.
- `python/tests/test_clpx_reverse_bound.py`: full-domain semantics, ramp branch
  identities, bounded-error injection, and failure-budget composition checks.

Four tests pass. They verify both decomposition layouts over all 65,536 residues,
exercise ramp wrap boundaries, and replay 520 edge/random integers per layout
with adversarial-sign primitive errors at 90% of the permitted envelope cutoff.
These are deterministic/model validation tests, not statistical evidence for
2^-40 failure probability. Cryptographic implementation and benchmark parameters
were not changed by this analysis.

Reproduce from Parameter-Selection:

```sh
PYTHONPATH=python python3 -m unittest discover -s python/tests -p test_clpx_reverse_bound.py -v
python3 python/CLPXreverse.py --tfhepp /home/ubuntu/sources/SS-workspace/TFHEpp --output /tmp/clpx-measured.json
python3 python/CLPXreverse.py --tfhepp /home/ubuntu/sources/SS-workspace/TFHEpp --basebit 2 --half-n 832 --half-alpha-bits 19 --output /tmp/clpx-candidate.json
```

The reports record the Parameter-Selection base revision, analysis-file hashes,
and analyzed C++ source hashes. The circuit is traced from this checkout's
`include/clpx/bfv-clpx.hpp`, `include/tfhe/homdecomp.hpp`, and
`include/tfhe/gatebootstrapping.hpp`. The local Nagai draft
`/home/ubuntu/sources/noiseproof-workspace/Nagai-InternationalConference/iacrtrans/sub/algo_B2T.tex`
provides the algorithm context; differences in scaling are resolved using the
executable C++ tables rather than assuming the draft and implementation agree.
