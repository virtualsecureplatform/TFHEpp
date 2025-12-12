# Core API

This page summarizes the main public types and functions. Most users only need
`include/tfhe++.hpp`.

## Namespaces and headers

All symbols are in namespace `TFHEpp`. The umbrella header is:

```cpp
#include <tfhe++.hpp>
```

## Types

TFHEpp uses `std::array`‑based fixed‑size types for strong compile‑time checking.
Definitions live in `include/params.hpp`.

- `TLWE<P>` — vector ciphertext of length `P::k * P::n + 1`.
- `TRLWE<P>` — ring ciphertext, `(P::k + 1)` polynomials of degree `P::n`.
- `TRGSW<P>` — GSW ciphertext over `TRLWE<P>`.
- `Polynomial<P>` and related FFT/NTT views.

`P` is a parameter struct (e.g. `lvl1param`) selected by build options.

## Keys

### SecretKey

```cpp
TFHEpp::SecretKey sk;
```

`SecretKey` holds the LWE secrets for each level. It is cheap to construct and
is required for encryption/decryption and evaluation‑key generation.

### EvalKey (evaluation / cloud key)

```cpp
TFHEpp::EvalKey ek;
ek.emplacebkfft<lvl01param>(sk);
ek.emplaceiksk<lvl10param>(sk);
```

`EvalKey` stores bootstrapping keys, key‑switching keys, circuit‑bootstrapping
keys, and other auxiliary keys. You only need to generate the pieces required
by the operations you plan to use.

Common generators (see `include/cloudkey.hpp`):

- `emplacebk<P>(sk)` / `emplacebkfft<P>(sk)` / `emplacebkntt<P>(sk)`
- `emplaceiksk<P>(sk)` — identity key switching.
- `emplacesubiksk<P>(sk)` — subset/partial key switching.
- `emplaceprivksk<P>(name, func_poly, sk)` — private key switching for circuit
  bootstrapping.
- `emplaceahk<P>(sk)` — annihilate key switching support.
- `emplacecbsk<P>(sk)` — circuit‑bootstrapping scheme‑switch key.

Low‑level accessors: `ek.getbkfft<P>()`, `ek.getiksk<P>()`, etc.

## Encryption and decryption

### Boolean TLWE

Batch encrypt/decrypt bits at a given level:

```cpp
auto c = bootsSymEncrypt<lvl1param>(plain_bits, sk);
auto p = bootsSymDecrypt<lvl1param>(c, sk);
```

Single‑ciphertext primitives:

- `tlweSymEncrypt<P>(torus_message, sk)`
- `tlweSymDecrypt<P>(ciphertext, sk)`
- `tlweSymIntEncrypt<P, plain_modulus>(integer, sk)`
- `tlweSymIntDecrypt<P, plain_modulus>(ciphertext, sk)`

### TRLWE / TRGSW

Ring encryption lives in:

- `include/trlwe.hpp` — `trlweSymEncrypt`, `trlweSymIntEncrypt`, decryptors.
- `include/trgsw.hpp` — `trgswSymEncrypt` and helpers.

These are used for LUT bootstrapping and circuit bootstrapping.

## Boolean gates

High‑level gates are in `include/gate.hpp`. They consume TLWE ciphertexts and
an `EvalKey`:

```cpp
HomNAND(res, ca, cb, ek);
HomAND(res, ca, cb, ek);
HomXOR(res, ca, cb, ek);
HomNOT(res, ca);
HomMUX(res, cs, c1, c0, ek);
```

Supported 2‑input gates include:

`AND`, `OR`, `XOR`, `NAND`, `NOR`, `XNOR`,
and negated‑input variants `ANDYN`, `ANDNY`, `ORYN`, `ORNY`.

Most gates have two overloads differing in whether identity key switching is
applied before or after blind rotation. If unsure, use the overload that takes
`TLWE<lvl1param>` inputs and returns `TLWE<lvl1param>`; it tends to have lower
noise growth.

## Gate bootstrapping / LUTs

Low‑level bootstrapping is in `include/gatebootstrapping.hpp`. The most common
entry points are:

- `GateBootstrappingTLWE2TLWE<P>(out, in, ek.getbkfft<P>(), testvector)`
- `GateBootstrappingManyLUT<P, num_out>(outs, in, ek.getbkfft<P>(), testvector)`

Here `testvector` encodes a LUT over the torus. For concrete examples, see
`include/aes.hpp` and tests like `test/gatebootstrapping*.cpp`.

