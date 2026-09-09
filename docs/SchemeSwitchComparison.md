# Native 8-bit scheme-switch multiplication comparison

The benchmark compares this checkout's TFHE-to-CLPX-to-TFHE path with its
integer-only TFHE-to-CKKS-to-TFHE path. Both multiply unsigned 8-bit operands
and verify every bit of the exact 16-bit product. It measures the implementations
on this machine, not the asymptotic schemes or a SIMD CKKS implementation.

## Measured results

Both paths passed all 256 products without OOM. Times below are means per
product, except key generation. Both operands are included in forward time.

| Measurement | Integer CKKS | CLPX |
|---|---:|---:|
| Forward switch (both operands) | 1121.627 ms | 41.998 ms |
| Multiplication + relinearization | 373.453 ms | 0.045 ms |
| Reverse switch | 145.133 ms | 953.787 ms |
| Total evaluation | 1640.213 ms | 995.830 ms |
| Key generation (once) | 8.601 s | 1.979 s |
| Peak resident memory | 14.053 GiB | 0.879 GiB |

CLPX was **1.65 times faster overall** and
used **15.99 times less peak memory** in this run.
Its reverse switch accounts for 95.8% of evaluation time.
For CKKS, forward switching accounts for 68.4%; its reverse stage was
6.57 times faster than CLPX's, but returns one integer TLWE instead of 16 Boolean TLWEs.
The maximum measured CKKS integer error was 0.000213026, below 0.5.

These are native-interface performance results, **not a claim that both paths
have a verified <=2^-40 end-to-end failure bound**. See the scope below.
Raw logs, source hashes, compiler flags, and derived measurements are in
[the machine-readable run record](scheme-switch-comparison.json).

## Method

One shared executable, `test/clpx/scheme_switch_comparison.cpp`, generates the
same 256 pairs for both paths: eight boundary cases and 248 pairs from
`std::mt19937(20260905)`. Operand draws are sequenced explicitly. Each process
generates fresh cryptographic keys and fresh input encryptions. All evaluation
keys remain resident and are reused for the 256 products. The two benchmark
processes run sequentially, with `OMP_NUM_THREADS=8` and a 24 GiB address-space
limit. There is one measured key set per scheme; results are averages over the
256 evaluations, not repeated-process confidence intervals.

Forward time includes switching both operands. Multiplication time includes
relinearization. Reverse CKKS time includes its required public scaling,
rescaling, modulus reduction, and coefficient-to-TFHE key switch. Outer temporary allocation, input encryption, and output decryption/verification are
outside the three evaluation timers. Key generation includes secret keys and
all evaluation keys. Peak RSS covers the entire process including key setup.

The common build uses GCC with `-O3 -march=native`, OpenMP, LTO, BLAKE3, and
SPQLIOS_INTL on an AMD Ryzen 9 9950X3D (16 cores, 32 logical CPUs). The default
TFHEpp parameter family is active, with block-binary and key bundling disabled.
The corrected BLAKE3 stream implementation is used for both paths.

## Representation and parameters

| Property | Integer CKKS path | CLPX path |
|---|---|---|
| Input per operand | One lvl2 integer TLWE | Eight lvl1 Boolean TLWEs |
| Output | One lvl2 16-bit integer TLWE | Sixteen lvl1 Boolean TLWEs |
| Arithmetic ring dimension | 16,384 | 2,048 |
| Arithmetic ciphertext modulus | 256 -> 192 -> 128 active bits | 64 bits |
| Arithmetic backend | Exact integer/chunked NTT | Existing CLPX FFT backend |
| Forward parameters | Constant-phase switching | 4 LUTs, shift 4, w=20; gadget 3 x 13 bits |
| Reverse parameters | q=2^48, scale=2^40, basebit 2, 24 levels | numdigit 5, basebit 4; lvl2 gadget 4 x 10 bits |

This is a comparison with each implementation's native TFHE interface. It does
not include converting Boolean ciphertexts into the CKKS path's integer input,
or extracting Boolean output ciphertexts from that integer output. Such costs
must be included for an application that requires an identical Boolean interface.
The secret-key parameter in the CLPX reverse API is unused by its implementation;
secret keys are used only for setup, input encryption, and result checking.

## Security and failure scope

A subsequent [whole-reverse envelope analysis](CLPXReverseBound.md) now covers
the measured 16-bit reverse circuit. Its conservative conditional bound does
not certify the measured parameters at 2^-40; a retuned analytical candidate
passes under explicit primitive-tail assumptions. Neither changes the measured
benchmark results below into an implementation-level correctness certificate.


Parameter-Selection's lattice-estimator gives minimum classical RC.BDGL16
four-attack costs of 131.54 bits for the CKKS path and 129.46 bits for CLPX.
CLPX's individual estimates are 129.65 bits for lvl1, 129.46 bits for lvlhalf,
and 131.54 bits for lvl2. These are generic unbounded-sample LWE screens, not
complete protocol security certificates. See [CLPX security results](clpx-comparison-security.json)
and [CKKS security results](integer4-security.json).

The CKKS circuit has the existing ideal-normal-sampler coefficient bound:
conditional error 0.439649 integer units < 0.5 and whole-run failure <=2^-46.09
for one key generation and 256 fresh-input multiplications. See
[the bound](integer8-failure-bound.json).

The installed Parameter-Selection CLPX preset has an older 4 x 9-bit forward
gadget. `test/clpx/comparison_noise.py` overrides it with the measured 3 x 13-bit
gadget without modifying Parameter-Selection. Its forward-plus-multiplication
screen is favorable, but it does not certify the full reverse circuit. In
particular, the installed reverse estimator models a different digit-recovery
sequence and two carry digits; the measured circuit uses a rounded-digit PBS
and one base-16 carry digit. Intermediate PBS semantic failures must be bounded
before claiming <=2^-40 end-to-end. A tiny scalar Gaussian tail from the paper
threshold is not that bound. See [the noise screen](clpx-comparison-noise.json).

The upstream [Parameter-Selection documentation](https://github.com/virtualsecureplatform/Parameter-Selection)
describes the CLPX noise tools; the results here pin the local checkout revision
and explicitly account for its mismatch with the measured implementation.

## Reproduction

From the TFHEpp directory:

```sh
cmake -S . -B build -DENABLE_TEST=ON -DUSE_BLOCK_BINARY=OFF -DUSE_BLAKE3=ON
cmake --build build --target scheme_switch_comparison -j2
OMP_NUM_THREADS=8 /usr/bin/time -v prlimit --as=25769803776 -- build/test/clpx/scheme_switch_comparison clpx
OMP_NUM_THREADS=8 /usr/bin/time -v prlimit --as=25769803776 -- build/test/clpx/scheme_switch_comparison ckks
python3 test/clpx/comparison_noise.py --parameter-selection /home/ubuntu/sources/CKKS-workspace/Parameter-Selection --output docs/clpx-comparison-noise.json
```
