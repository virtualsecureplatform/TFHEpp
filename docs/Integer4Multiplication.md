# TFHE → CKKS multiplication → TFHE for 4-bit integers

This path takes unsigned inputs `x,y` in `[0,15]`, each encoded as `x/256`
in a 64-bit TFHE torus, and returns the exact 8-bit product `x*y` in the
same encoding. It uses one value per CKKS ciphertext. It is an integer-only
construction, not a general conversion to canonical approximate real slots.

The default components pass the local classical 128-bit security screen.
A conservative analytical bound gives failure probability at most
`2^-77.8296`, below the requested `2^-40`, for one key generation and up to
256 fresh-input multiplications. This probability uses the standard ideal
Gaussian sampler model; it is not inferred from a small empirical failure
count. The arithmetic evaluation is exact modular integer/NTT arithmetic.

Implementation: [integer_multiplication.hpp](../include/ckks/integer_multiplication.hpp).

## Why no EvalMod is needed

Set `M=256`. Homomorphic partial decryption produces constant CKKS plaintexts
whose ideal values are lifts `u=x/M+k` and `v=y/M+l`, where `k,l` are integers.
Then

```text
M*u*v = x*y/M + x*l + y*k + M*k*l.
```

The last three terms are integers. Final conversion to the TFHE torus
removes them modulo one. Since `x*y<=225<256`, nearest-integer decoding
recovers the full product, not merely a truncated four-bit result.
For approximate real inputs the cross terms do not disappear; this shortcut
must not be used for arbitrary real arithmetic.

Each of the 2048 TFHE secret coefficients is encrypted as a constant
polynomial under the **full dense CKKS secret**. Public integer scalar
combinations compute `b-a*s` without slot encoding, rotations, or dense
linear transforms. The CKKS plaintext has the same constant ideal value in
every slot; coefficient noise is retained and bounded, not assumed zero.
The output uses ordinary direct LWE switching keys, not the rejected
zero-padded TFHE ring-secret construction.

All new keys store their complete uniformly sampled masks. The older
`CKKSSeedStream` compresses its state to 64 bits and uses SplitMix, not a
cryptographic PRG; this path deliberately does not use it. The tested build
uses TFHEpp's BLAKE3-based generator seeded from `std::random_device`.
This implementation also corrects that generator's refill offset: BLAKE3
seek offsets are bytes, so a 4096-byte buffer must advance by 4096 bytes,
not 512 64-bit words. The regression `ckks_prng_stream` reproduces the old
failure at the first refill and checks four buffers against the reference
XOF for both 32-bit and 64-bit output words. Keys generated with the old
overlapping stream must be regenerated; the calibrated validation uses
fresh keys after this fix.

## Parameters and security

| Component | Parameters | Weakest estimated attack |
|-----------|------------|-------------------------:|
| CKKS | N=16384, q=2^256, dense ternary secret, sigma=3.2 | 233.54 bits |
| TFHE `lvl2param` | n=2048, q=2^64, ternary secret, sigma=8192 (`alpha=2^-51`) | 131.54 bits |

Parameter-Selection's pinned lattice-estimator evaluated primal uSVP,
primal BDD, dual, and dual hybrid with the classical `RC.BDGL16` model and
unbounded samples. All eight estimates completed. The minimum is **131.54
bits**. See [integer4-security.json](integer4-security.json) for per-attack
costs and exact repository revisions. Lower-modulus relinearization rows
use the same full CKKS secret and the 3.2-unit noise floor.

As with conventional FHE parameter estimates, this is a heuristic classical
LWE/RLWE security screen, with the usual circular/KDM assumption for
evaluation keys. It is not a quantum-security claim or a proof against all
possible attacks. The [lattice-estimator](https://github.com/malb/lattice-estimator)
documents the attack-cost estimation model.

The active CKKS levels are:

```text
secret-coefficient encryptions: q=256, scale=64
public 64-bit scalar accumulation, rescale 64: q=192, scale=64
exact tensor multiplication, rescale 64, relinearize: q=128, scale=64
multiply by M^2=65536; rescale 24: q=104, scale=40
intentional level reduction: q=48, scale=40
direct TLWE extraction/switch: BaseBit=2, Levels=24, target q=2^64
```

Extraction divides the represented coefficient by `2^(48-40)=256`.
Consequently the public `M^2` multiplication implements the required `M*u*v`
torus phase. The analytical check verifies headroom before the intentional
final modulus reduction.

The exact NTT product uses 16-bit digits and prime
`p=4611686018427322369`. For the 192-bit tensor inputs it combines at most
12 digit products per inverse transform. Its signed coefficient bound is
`12*N*65535^2 < p/2`, enforced by a static assertion. Thus modular NTT
coefficients reconstruct the integer digit products without floating-point
roundoff or prime wraparound. Integer rescaling rounds with at most a
half-unit error per ciphertext component coefficient.

## Failure bound

The executable bound is
[integer4_failure_bound.py](../test/ckks/integer4_failure_bound.py), with
results in [integer4-failure-bound.json](integer4-failure-bound.json).
It uses exact rational arithmetic for the conditional error calculation.
Let `N=16384`, `n=2048`, `Delta=2^64`, and `U=n+1`. For each sampler output,
take a 12-standard-deviation cutoff and allow one extra integer unit for
rounding. This gives coefficient bounds `E_C=40` for CKKS and `E_T=98305`
for TFHE. The total number of Gaussian samples is at most:

| Source | Samples |
|--------|--------:|
| TFHE-secret coefficient encryptions | 33,554,432 |
| Relinearization rows | 262,144 |
| Direct output switching rows | 786,432 |
| Fresh inputs for 256 multiplications | 512 |

For a standard normal `Z`, `Pr(|Z|>12) <= 2*exp(-72)`. A union bound over
these 34,603,520 samples gives `log2(p) <= -77.8296`. The conditional
correctness argument does not require independence between ciphertext
errors, coefficients, selected key rows, or reused evaluation keys.

Here are the conditional error bounds, in decoded integer units:

| Contribution | Bound |
|--------------|------:|
| Original TFHE input noise after multiplication | 1.432e-6 |
| CKKS phase conversion errors and their product | 1.312e-6 |
| Tensor rescaling and relinearization | 5.246e-6 |
| Final rescaling | 7.452e-9 |
| Direct output LWE switching | 5.365e-7 |
| **Total** | **8.533e-6 < 0.5** |

More explicitly, put `eta=E_T/2^64` and
`A=(n*E_C+(1+N)/2)/Delta`. The ideal unwrapped TFHE phase has absolute value
at most `U`, and every coefficient error after CKKS phase conversion is at
most `A`. Multiplication contributes at most
`2*U*eta+eta^2 + 2*U*A+N*A^2` before the output factor `M^2`.
The `N*A^2` term bounds the entire noise-polynomial convolution, including
nonconstant coefficients.

The three rounded tensor components have secret weights `1,s,s^2`.
Their error is bounded by `(1+N)^2/(2*Delta)`. Relinearization contributes
at most `N*16*128*E_C/Delta`: 16 radix-256 rows, digit magnitude at most 128,
and a negacyclic convolution of N bounded terms. Multiply both bounds by
`M^2`. Final rescaling adds `(1+N)/(2*2^40)`. Direct output switching adds
`M*N*24*E_T/2^64`. These bounds allow every error to have the worst sign;
they do not extrapolate an empirical standard deviation.

Under the good sampler event, the total error is below half an integer for
every valid input pair, so the only modeled failure event is a sampler
cutoff violation. Reusing keys is included. Change `--multiplications` to
recompute the union bound for a different number of fresh input pairs.
The bound does not cover incorrectly encoded inputs, subsequent operations
on the output, or implementation faults. The C++ standard library normal
sampler and pseudorandom generator have not been formally verified here.

## Build, reproduce, and use

From the TFHEpp root:

```sh
cmake -S . -B build -DENABLE_TEST=ON -DUSE_BLOCK_BINARY=OFF -DUSE_BLAKE3=ON
cmake --build build --target ckks_integer_multiplication ckks_prng_stream -j2
build/test/ckks/ckks_prng_stream

# Small-ring exhaustive arithmetic regression (not a security-sized test).
OMP_NUM_THREADS=4 build/test/ckks/ckks_integer_multiplication

# Default calibrated dimensions: eight edge cases, then all 256 pairs.
OMP_NUM_THREADS=8 /usr/bin/time -v prlimit --as=25769803776 --core=0 -- \
  build/test/ckks/ckks_integer_multiplication --full
OMP_NUM_THREADS=8 /usr/bin/time -v prlimit --as=25769803776 --core=0 -- \
  build/test/ckks/ckks_integer_multiplication --full-exhaustive

python3 test/ckks/integer4_failure_bound.py --output docs/integer4-failure-bound.json
sage -python test/ckks/calibrate_scheme_switching.py \
  --parameter-selection /path/to/Parameter-Selection \
  --case ckks-integer4 --case tfhe2048-default \
  --output docs/integer4-security.json --jobs 4
```

The Sage container command in
[SchemeSwitchingCalibration.md](SchemeSwitchingCalibration.md) also works;
use the two case names and output path above.

```cpp
using C = TFHEpp::CKKSInteger4Param;
using L = TFHEpp::lvl2param;
auto ckks_secret = std::make_unique<TFHEpp::Key<C>>();
TFHEpp::Key<L> tfhe_secret{};
TFHEpp::keyGen<C>(*ckks_secret);
TFHEpp::keyGen<L>(tfhe_secret);
auto evaluation_key = std::make_unique<TFHEpp::TFHECKKSInteger4Key<C, L>>();
TFHEpp::TFHECKKSInteger4KeyGen(*evaluation_key, *ckks_secret, tfhe_secret);

TFHEpp::TLWE<L> a{}, b{}, product{};
TFHEpp::TFHEEncryptInteger4(a, 15, tfhe_secret);
TFHEpp::TFHEEncryptInteger4(b, 15, tfhe_secret);
TFHEpp::TFHECKKSInteger4Multiply(product, a, b, *evaluation_key);
auto result = TFHEpp::TFHEDecryptInteger8(product, tfhe_secret); // 225
```

Only setup, encryption, and final decryption use secret keys. The evaluator
needs the evaluation key and input ciphertexts. Its keys remain reusable.
The default TFHE bit/gate plaintext convention is not used: these helper
functions explicitly encode the 256-bin integer torus convention.

The direct key payload is 12,891,193,344 bytes (12.0059 GiB). The full-mask
phase key adds 2 GiB, and relinearization keys add 16 MiB. Arithmetic scratch
is additional. The `2^-40` target is supported by the analytical bound above;
empirical tests check implementation behavior rather than measuring such a
small failure probability.

The full-size exhaustive run with fresh keys after the PRNG correction
passed all **256 input pairs**, under the 24 GiB address-space cap:

| Measurement | Result |
|-------------|-------:|
| Peak RSS | 14.05 GiB (14,732,460 KiB) |
| Key generation | 8.66 s |
| 256 complete evaluations | 426.78 s |
| Total wall time | 7 min 15.72 s |
| Maximum integer-unit error | 2.285e-9 |

This first implementation evaluates one product per CKKS ciphertext, about
1.67 seconds per product in this run; it does not yet provide SIMD packing.
The [validation record](integer4-validation.json) contains the complete
progress output and measurements. Small-ring exhaustive tests, the exact
NTT/schoolbook product comparison, malformed-input checks, the BLAKE3 stream
regression, and the existing CKKS/TFHE switching regression also pass.
