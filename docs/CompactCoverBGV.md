# N=65536 compact-cover BGV bootstrap

`include/bfv/compact-cover-bgv.hpp` implements the bounded-memory carrier and
transition layer. `include/bfv/compact-cover-bgv-bootstrap.hpp` implements a
complete scalar-message BGV refresh at the certified N=65536 parameters.

## Concrete profile

The checked-in profile is the thin, non-cyclic power-of-two schedule

```text
ring degree                 65536
cyclotomic index            131072
plaintext prime             65537
Frobenius order             2
plaintext slots             32768
hypercube dimensions        (16384, 2)
hypercube generators        (5, -1)
baby dimensions             (2, 91)
giant dimensions            (1, 181)
distinct switch exponents   362
peak live ciphertexts       368
```

The C++ manifest is generated from the mixed-radix formulas.  Its test checks
both the count of distinct automorphisms and that they generate all 65536 odd
cyclotomic exponents; the 362 exponents are not maintained as an unrelated
literal table.

The RNS basis contains fifteen 61-bit primes supporting degree-65536
negacyclic NTTs. Every prime is also one modulo `65537^2`, permitting exact
plaintext-preserving BGV operations. A full width-368 ciphertext occupies
about 5.39 GiB, while the scalar specialization uses width one and a 75 MiB
five-row bootstrap key.

## Implemented operations

- heap-backed `[limb][frontier label][NTT slot]` storage;
- checked memory and disk resource estimates;
- HEXL pointwise addition and multiplication, with scalar fallbacks;
- exact Galois relabeling in the negacyclic NTT representation;
- forward and inverse degree-65536 transforms;
- an exact generator for the selected thin-BGV switch-key manifest;
- a versioned, disk-backed transition-key format;
- slice-at-a-time seeded mask expansion and body loading; and
- an optional seeded per-limb transition for integration tests;
- a stored-uniform transition format adding no PRG assumption; and
- balanced coefficient decomposition over the complete 915-bit CRT modulus.

For gadget base `B`, transition row `r` uses the convention

```text
key_body - key_mask * target_secret = -B^r * source_secret + error.
```

The certified implementation reconstructs each centered coefficient modulo
the complete RNS product before taking its five balanced digits. It transforms
each shared digit polynomial and consumes one key row at a time. The test
checks the resulting phase identity in every NTT slot and across multiple RNS
limbs.

Every mask slice has an independent 256-bit seed.  Normal TFHEpp builds expand
it with BLAKE3 and use rejection sampling to obtain unbiased RNS residues.
Builds explicitly configured with `USE_BLAKE3=OFF` retain a deterministic
SplitMix fallback for tests; that fallback must not be used to generate
cryptographic evaluation keys.  Storing truly uniform masks instead remains
the information-theoretic alternative to seeded compression.

## Scalar direct bootstrap

For an input BGV phase

```text
b - a*s = m + p*e,
```

the scalar circuit first multiplies both ciphertext components by `p`. Its
stored-uniform transition rows have phase

```text
-B^r*s + p^2*e_r.
```

Full-modulus key switching therefore produces phase

```text
p*m + p^2*(e + sum_r d_r*e_r).
```

Multiplication of both output components by `p^-1 mod Q` gives an ordinary
plaintext-`p` BGV ciphertext of `m`. CBD(20) is bounded, and the checked
worst-case absolute error is below `2^221`, far below the 915-bit modulus.
There is no rounding failure or digit-extraction failure in this scalar path.
The unlimited-sample attack-cost proxy for the actual `p^2*CBD(20)`
Binary-NTT source is 243.53 bits; the checked certificate retains 242.53 bits
after its explicit reduction reserve.

Key generation samples a weight-32 ternary operational secret, a shared
Binary-NTT witness, unit pivots, and the public quadratic hint used for
two-component ciphertext multiplication. The same evaluation key is reused by
repeated bootstraps.

## Build and run

```bash
cmake -S . -B build-compact-cover \
  -DENABLE_TEST=ON -DUSE_HEXL=ON -DUSE_AVX512=ON \
  -DCMAKE_POLICY_VERSION_MINIMUM=3.5
cmake --build build-compact-cover \
  --target compact_cover_bgv_65536 compact_cover_bgv_bootstrap_65536 -j4
./build-compact-cover/test/bfv/compact_cover_bgv_65536
./build-compact-cover/test/bfv/compact_cover_bgv_bootstrap_65536
```

The full-size bootstrap test checks encryption, quadratic-hint multiplication,
resumable filesystem key generation, one refresh, refresh of the refresh
output, and refresh after multiplication.

## Scope boundary

The executable API encrypts one integer as a constant polynomial. Arbitrary
packed values in all 32768 degree-two slots still require the general
width-368 coefficient-to-slot/digit-extraction circuit described by the thin
schedule. Public-key ring-Regev encryption is also outside this milestone;
the implemented interface is secret-key FHE with public evaluation material.

The security statement is under ordinary Binary-NTT RLWE and the checked
quadratic/affine compiler. A reduction from conventional RLWE to Binary-NTT
RLWE remains a separate research premise.
