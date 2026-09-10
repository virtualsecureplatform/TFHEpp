# CLPX t=3 performance experiment

This experiment compares the measured original CLPX configuration with the
fixed-760 candidate from [the reverse-bound analysis](CLPXReverseBound.md).
Both multiply unsigned 8-bit operands and verify 16 Boolean output ciphertexts.

| Parameter | Baseline | Candidate |
|---|---|---|
| lvl2h identity switch | t=7, basebit=2 | t=3, basebit=5 |
| HomDecomp | numdigit=5, basebit=4 | numdigit=9, basebit=2 |
| Intermediate TLWE | n=760, alpha=2^-17 | unchanged |
| Ring dimensions, secret distributions, fresh ring noise | existing parameters | unchanged |

These are two coordinated changes, so the total runtime difference does not
isolate the identity-switch optimization alone. Smaller HomDecomp digits require
more decomposition operations. The balanced lvl2h key grows from 87,277,568 bytes
to 299,237,376 bytes (0.0813 to 0.2787 GiB), while each identity switch processes
fewer levels.

## Results

Both configurations passed all 256 products without OOM. Evaluation timings
are means per product; key generation is once per run.

| Measurement | Baseline | t=3 candidate |
|---|---:|---:|
| Forward (both operands) | 41.562 ms | 41.938 ms |
| Multiplication + relinearization | 0.046 ms | 0.045 ms |
| Reverse (16 output bits) | 996.929 ms | 1010.436 ms |
| Total evaluation | 1038.537 ms | 1052.419 ms |
| Key generation | 1.985 s | 2.093 s |
| Peak RSS | 0.879 GiB | 1.076 GiB |

The candidate took 1.34% longer overall in this run,
so performance is approximately unchanged. A difference this small is not
sufficient to establish a reliable regression from one run per configuration.
Peak RAM increased by 0.198 GiB,
consistent with the larger lvl2h key. This experiment validates the candidate's
functional execution at the calibrated dimensions, but not its conditional
rare-event correctness bound. [Raw logs and source hashes](clpx-t3-performance.json).

## Method

A fresh baseline and candidate process run sequentially, each with a fresh key
set and the same 256 operand pairs: eight boundary cases and 248 deterministic
random pairs. Both use eight OpenMP threads and a 24 GiB address-space cap.
Input encryption and output decryption/checking are excluded from evaluation
timers. Forward timing includes both operands; multiplication includes
relinearization; reverse timing includes all 16 output bits. Peak RSS covers
the entire process, including key generation. One key set is measured per
configuration; timings are not confidence intervals from repeated runs.

The candidate source is an isolated copy in `../experiments/clpx-t3`, with only
the changes in [the reproduction patch](clpx-t3-benchmark.patch). The exact tested candidate changes have now been applied to the main checkout.
The original baseline is commit `1451663`; both measured executables use the
same `-O3 -march=native`, LTO, OpenMP, BLAKE3, and SPQLIOS_INTL build settings.

## Reproduction

The main checkout now contains the tested candidate. Build and run it with:

```sh
cmake -S . -B build -DENABLE_TEST=ON -DUSE_BLOCK_BINARY=OFF -DUSE_BLAKE3=ON
cmake --build build --target scheme_switch_comparison -j2
OMP_NUM_THREADS=8 /usr/bin/time -v prlimit --as=25769803776 -- build/test/clpx/scheme_switch_comparison clpx
```

For the baseline, run the same command from a separate checkout of commit
`1451663`. The patch also reproduces the candidate when applied to that commit. Do not
run the two benchmarks simultaneously.

The cryptographic instances retain the prior minimum classical security screen
of 129.46 bits. The candidate's <=2^-43.47 whole-run reverse failure bound remains
conditional on the stated input/primitive tail and FFT-error assumptions.
Empirical correctness in this benchmark cannot establish that rare-event bound.
