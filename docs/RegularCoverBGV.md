# Regular-cover same-key BGV prototype

TFHEpp now contains a small executable prototype of the regular-cover BGV
construction formalized in FormalProof4FHE.

The implementation is in `include/bfv/regular-cover-bgv.hpp`, with an
end-to-end test in `test/bfv/regular_cover_bgv.cpp`.

## Implemented construction

The prototype uses a split base ring and a cyclic regular cover. It implements:

- the regular lifted automorphism action;
- the full-dimensional invariant embedding;
- exact Binary-NTT affine completion;
- a coefficient-small operational secret and public quadratic hint;
- ordinary noisy encryption under one fixed operational secret;
- exact assembly and compilation of ordinary base-RLWE rows into cover rows;
- two-component quadratic-hint multiplication;
- direct uniform-mask gadget rows for `g*sigma(z)`;
- lifted automorphism key switching;
- a direct noisy gadget boot key for `g*z`;
- homomorphic phase evaluation by identity key switching;
- bounded digit removal by an interpolated polynomial; and
- rerandomization with a fresh noisy zero encryption.

The same noisy boot key is reused for repeated refreshes under the same
operational secret.

## Test profile

The executable test deliberately uses a small auditable profile:

```text
base cyclotomic     = F_q[X]/(X^2+1)
base split degree N = 2
cover group size    = 2
ciphertext modulus  = 65537
plaintext modulus   = 2
secret coefficients = {-1,0,1}, transformed to the two NTT roots
encryption error     = [-2,2]
boot-row error       = [-1,1]
gadget base          = 2
```

The nontrivial Galois automorphism swaps the two roots of `X^2+1`, so this is
an actual cyclotomic action rather than an arbitrary split-coordinate
permutation. The digit-removal polynomial covers the conservative phase interval
created by the input error and all boot-row errors. The test checks both
plaintexts over multiple trials and performs two consecutive refreshes using
the same key material.

It additionally checks:

- the public quadratic equation;
- binary completion of every witness coordinate;
- nonzero boot-row errors with the intended plaintext scaling;
- quadratic-hint multiplication;
- the lifted automorphism action;
- automorphism switching back to the same key; and
- the exact ordinary-source regular-cover security compiler identity.

## Build and run

With tests enabled, CMake discovers the test automatically:

```bash
cmake -S . -B build -DENABLE_TEST=ON
cmake --build build --target regular_cover_bgv
./build/test/bfv/regular_cover_bgv
```

## Scope and limitations

This is a functional confirmation of the proved construction, not a practical
parameter set. The implementation stores split/NTT coordinates directly and
uses a cyclic action so that every algebraic step remains visible.

The proof-derived secure candidate is extraordinarily oversized. The current
Parameter-Selection screen indicates that the original `N=512`, 900-bit
example is insecure for the ordinary source assumption, while a conservative
`N=65536` source passes the LWE attack-cost proxy but makes a full regular
cover and its dense-binary public key far too large for execution.

Production integration would require:

- replacing the split-coordinate carrier with TFHEpp's RNS polynomial type;
- generating the exact automorphism subgroup used by the selected bootstrap;
- instantiating native centered rounding and modulus switching componentwise;
- importing the native bootstrap noise recurrence; and
- selecting a smaller automorphism subgroup or a more compact statistically
  hiding public-key encoding if performance becomes relevant.

These limitations do not affect the functional result: the same-key regular
cover supports noisy boot rows, automorphism switching, multiplication, and
repeated refresh without an independent key chain.
