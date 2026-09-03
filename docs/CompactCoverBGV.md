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

The RNS basis contains twenty-three approximately 61-bit primes supporting degree-65536
negacyclic NTTs. Every prime is also one modulo `65537^2`, permitting exact
plaintext-preserving BGV operations. The scalar specialization uses width one;
its complete stored-uniform evaluation directory is approximately 8.18 GiB.

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
- balanced coefficient decomposition over the active CRT modulus.

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

## Scalar bootstrap

For an input BGV phase

```text
b - a*s = m + p*e,
```

the evaluator first modulus-switches the tired ciphertext to one RNS limb. A
two-row cross-modulus transition centers and scales its components by `p` and
lifts the result to the full 23-limb modulus. Its rows have phase

```text
-B^r*s + p^2*e_r.
```

The resulting plaintext modulo `p^2` is

```text
p*m + carry(X),
```

where every coefficient of the carry is bounded by 23 at the accepted input
level. A 16-stage Galois trace projects the constant coefficient, dropping one
RNS limb halfway and one afterward. The odd degree-93 removal polynomial maps
`p*m+carry[0]` to `p*m`; exact division by `p` returns `m`.

The deterministic recurrence bounds the 13-limb output error by `2^641.98`
against capacity `2^775.61`. Before multiplication, refreshed ciphertexts are
dropped to two limbs. One quadratic-hint multiplication and level drop then
produces a one-limb error below `2^4.17`, against capacity near `2^44`, ready
for the next bootstrap. Dropping two refreshed operands to one limb and adding
them gives the exact certified bound 36, also inside the next bootstrap's
accepted input set. The multiplication cycle contracts by about 27.4 bits.

The unlimited-sample comparison proxy for the actual `p^2*CBD(20)`
Binary-NTT source is 133.44 bits; reserving one bit for reduction accounting
leaves 132.44 bits. This number is an attack-cost proxy, not a proof reduction
to coefficient-binary LWE.

Key generation samples a weight-32 ternary operational secret, a temporary
Binary-NTT witness, unit pivots, and the public quadratic hint. The public
directory contains the phase-lift rows, sixteen level-specific Galois keys,
the hint, parameter/certificate hashes, and checksums. The master witness is
erased before evaluation, and the same directory is reused indefinitely.
The production key-generation and encryption overloads use TFHEpp's
cryptographically seeded generator. Engine-taking overloads are retained for
deterministic conformance tests and callers that explicitly manage their own
cryptographic generator.

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

The full-size bootstrap test checks encryption, nested BGV modulus reduction,
resumable filesystem key generation, one refresh, refresh of the refresh
output, addition and multiplication of refreshed outputs, and refresh of each
result. It also reconstructs the secret phase after every bootstrap stage and
checks the measured error against the formal certificate. The test-only phase
audit does not expose the secret to the evaluator API.

The N=65536 substrate test also evaluates a sparse polynomial at selected odd
root powers and checks that `NegacyclicNTTPlan::forward` uses the same natural
slot ordering as the mathematical Lean backend. This closes the public NTT
ordering boundary for both HEXL and scalar transform implementations.

## Scope boundary

The executable API encrypts one integer as a constant polynomial. Arbitrary
packed values in all 32768 degree-two slots still require the general
width-368 coefficient-to-slot/digit-extraction circuit described by the thin
schedule. Public-key ring-Regev encryption is also outside this milestone;
the implemented interface is secret-key FHE with public evaluation material.

The security statement is under ordinary Binary-NTT RLWE and the checked
quadratic/affine compiler. A reduction from conventional RLWE to Binary-NTT
RLWE remains a separate research premise.

FormalProof4FHE computes the correctness recurrence over exact natural numbers
and proves the concrete digit-removal polynomial for every scalar message and
supported carry. The reported 132.44-bit value remains a computational
lattice-attack estimate and is not part of the Lean correctness predicate.
