# Scheme-switching parameter calibration

## Decision

For the subsequently specified 4-bit multiplication objective, the
[integer-only path](Integer4Multiplication.md) supersedes the general-path
component targets below. It eliminates EvalMod, uses a 256-bit CKKS modulus
at N=16384, and retains ordinary direct switching to `lvl2param`. See its
separate security estimates and analytical failure bound.

The previously measured 16.61 GiB roundtrip is a **functional prototype that
fails the 128-bit security screen**. Successful decryption and bounded memory
do not establish security. Do not use its parameters for confidential data.

The component targets selected for further implementation are CKKS
`N=32768, logQ<=832, sigma=3.2`, dense ternary secret, and TFHEpp's standard
`lvl2param` (`n=2048`, ternary secret, 64-bit torus, `alpha=2^-51`) through
the direct LWE switching path. The weakest measured component is 131.54 bits.
These are **component targets, not an accepted end-to-end parameter set**:
the reverse precision schedule must be retuned and validated. The runnable
prototype remains explicitly marked as rejected; its correctness measurements
must not be attributed to these targets.

## Reproducible estimates

The calculation uses the local `Parameter-Selection` checkout and the same
four attacks and classical `RC.BDGL16` cost model as its
`python/estimates/CKKS_lvl6.py`: primal uSVP, primal BDD, dual, and dual hybrid.
All four completed for every case. The table reports their minimum, using
unbounded samples as a conservative generic LWE screen.

| Instance | Unknown dimension | log2(q) | Integer noise sigma | Secret | Minimum bits |
|----------|------------------:|--------:|--------------------:|--------|-------------:|
| Prototype CKKS | 32768 | 1000 | 3.2 | ternary | **111.87 — reject** |
| Prototype compact key, message eliminated | 760 | 40 | 3.2 | binary | **61.59 — reject** |
| Standard fresh lvlhalf TLWE | 760 | 32 | 32768 | binary | 129.46 |
| Shorter CKKS chain | 32768 | 832 | 3.2 | ternary | **136.38** |
| Less conservative CKKS chain | 32768 | 880 | 3.2 | ternary | 128.22 |
| Larger CKKS ring | 65536 | 1000 | 3.2 | ternary | 242.83 |
| Standard lvl2 TLWE | 2048 | 64 | 8192 | ternary | **131.54** |
| Hypothetical compact key, full-modulus generic screen | 2048 | 46 | 3.2 | binary | 147.38 |
| Alternative small-noise TLWE | 2048 | 32 | 3.2 | binary | 214.16 |
| Alternative small-noise TLWE | 4096 | 32 | 3.2 | binary | 459.88 |

Machine-readable per-attack costs and revision identifiers are in
[scheme-switching-calibration.json](scheme-switching-calibration.json).
The pinned revisions are:

- Parameter-Selection: `9dcf0acb9adb8765c5aa2d5a9b1292a612083543`
- Its lattice-estimator submodule: `27a581bb8e9d49f5e9e2db315bd48ac769d5f5f5`

The [lattice-estimator](https://github.com/malb/lattice-estimator) provides
heuristic attack-cost estimates. This check follows the local four-attack
policy; it is not an exhaustive cryptanalysis, a quantum-security claim, or
a proof for evaluation keys with circular/KDM dependencies. Ring structure,
known-support secrets, implementation leakage, and decryption-oracle attacks
need separate treatment. A passing generic screen cannot certify the compact
padded-secret construction.

From the TFHEpp root, with Sage installed:

```sh
sage -python test/ckks/calibrate_scheme_switching.py \
  --parameter-selection /path/to/Parameter-Selection \
  --output docs/scheme-switching-calibration.json --jobs 4
```

On this machine, the existing Sage container is usable with:

```sh
apptainer exec --no-mount bind-paths,home,cwd \
  --env DOT_SAGE=/tmp/ss-calibration-sage --pwd /work \
  --bind /home/ubuntu/sources/CKKS-workspace/Parameter-Selection:/selection \
  --bind /home/ubuntu/sources/SS-workspace/TFHEpp:/work \
  /home/ubuntu/sources/CKKS-workspace/Parameter-Selection/python/sagemath.sif \
  sage -python /work/test/ckks/calibrate_scheme_switching.py \
  --parameter-selection /selection \
  --output /work/docs/scheme-switching-calibration.json --jobs 4
```

The mount exclusions avoid a local Apptainer default-bind failure. The
container uses SageMath 10.9. `--case NAME` restricts a run to a candidate.
Rejected candidates are expected calibration results; estimator errors cause
a nonzero exit and an incomplete result instead of a passing estimate.

## Why the compact key fails

`CKKSToTFHERingSwitchKeyGen` pads the 760-coefficient binary TFHE secret with
known zeros to the CKKS ring size. Its public seeded row has the form

```text
b = a * padded(s_TFHE) + 2^shift * s_CKKS + e  (mod 2^46).
```

With radix 256 and the current active level, the highest active row has
`shift=40`. Reducing that row modulo `2^40` eliminates its entire unknown
message, leaving a small-noise instance involving only the 760 unknown
TFHE coefficients. Seed compression does not conceal the uniformly generated
mask. Counting all 32768 padded positions as unknown would overestimate
security. The 61.59-bit result is a generic LWE attack-cost screen of the
exposed instance, not a demonstrated key-recovery attack implementation;
its ring-derived sample structure can require additional analysis.

The compact key uses the supplied CKKS noise (`sigma=3.2` here), not
`lvlhalfparam::alpha`. Consequently the 129.46-bit estimate for ordinary
fresh TLWE encryption cannot justify this switching key. Adding noise only
to the final TLWE outputs would not protect an already exposed public key.
Likewise, fixing only the CKKS modulus leaves this second failure intact.

The 2048-dimensional compact candidate passes the generic dimension/noise
screen, but its known-support ring structure remains uncalibrated. Prefer
the existing direct `CKKSToTFHEKeyGen` / `CKKSCoeffToTFHE` construction for
the next implementation: its rows are ordinary target TLWE encryptions and
use the target parameter's actual noise distribution.

## Precision and memory constraints

The existing reverse circuit consumes 60 bits for its linear transform,
five 58-bit basis multiplication levels, a 58-bit coefficient multiplication,
and nine 58-bit double-angle levels. Its output modulus is therefore

```text
1000 - 60 - (5 + 1 + 9)*58 = 70 bits.
```

Simply replacing 1000 with 832 would give **-98 bits**, so that change cannot
produce a working circuit. A modulus of 880 with 52-bit ciphertext and
coefficient scales and a 40-bit linear scale leaves 60 bits for the original
760-dimensional schedule, but has only a 0.22-bit security-screen margin and
has not been encrypted-validated. Increasing the TFHE dimension also changes
the unwrapped phase bound and hence the required approximation degree or
double-angle count. Neither modulus change is a drop-in calibrated fix.

For the preferred standard lvl2 destination, a direct key covering the
46-bit extraction level with `BaseBit=2, Levels=23` has size

```text
32768 * 23 * 2 * (2048 + 1) * sizeof(uint64_t)
= 24,708,120,576 bytes = 23.0112 GiB.
```

This is an exact payload calculation for the existing unseeded key type,
not measured peak RSS. Generate/use/release it before reverse-key generation.
The streamed transform still avoids the old 112 GiB plaintext plans. The
machine's approximately 91 GiB RAM makes this key feasible, but a new memory
cap must include ciphertexts, rotations, EvalMod keys, and scratch.

As an independent-noise variance screen, at most `32768*23` selected direct
key rows contribute fresh error per output. With `alpha=2^-51`, their
standard deviation in original message units (phase multiplied by 64) is
at most `64*sqrt(32768*23)*2^-51 = 2.47e-11`. This accounts only for key-row
noise; CKKS decoding, rescaling, modular reduction, numeric rounding, and
cross-output correlations are not covered. Full encrypted accuracy and
failure-probability validation remain necessary.

The larger-ring alternative has ample generic security margin while keeping
the 1000-bit modulus. It is not currently runnable through the exact compact
backend, which statically limits ring sizes to `2^15`. Its arithmetic support
and precision schedule would require changes, and the analogous direct lvl2
key alone would occupy 46.02 GiB. It is not the first choice for this
memory-sensitive implementation.
