# Advanced Features

This page collects TFHEpp features beyond basic Boolean gates.

## Programmable bootstrapping with multiple outputs

`GateBootstrappingManyLUT<P, num_out>` evaluates a LUT and returns `num_out`
TLWE ciphertexts in one bootstrapping. It is useful for SIMD‑style bit‑slicing,
S‑boxes, and packed evaluation.

See:

- API: `include/gatebootstrapping.hpp`
- Example use: `include/aes.hpp`

## Circuit bootstrapping (TLWE → TRGSW)

Circuit bootstrapping turns a TLWE ciphertext into a TRGSW ciphertext to enable
faster arithmetic over `TRLWE`/`TRGSW` or to refresh noise before complex
circuits.

Entry points (see `include/circuitbootstrapping.hpp`):

- `CircuitBootstrapping<brP, privksP>(trgsw, tlwe, ek)` (outputs TRGSW)
- `CircuitBootstrapping<brP, privksP>(trgswfft, tlwe, ek)` (outputs TRGSWFFT)
- `CircuitBootstrappingSub<iksP, bkP, privksP>(...)` for subset keys.

You must generate:

- a bootstrapping key: `ek.emplacebkfft<brP>(sk)`
- private key‑switching keys for CB:
  `ek.emplaceprivksk4cb<privksP>(sk)` (or subset variant)

## Subset / partial keys

With `-DUSE_SUBSET_KEY=ON`, TFHEpp can generate subset key‑switching keys
(`emplacesubiksk`, `emplacesubprivksk`) that allow evaluation with only part of
the secret key. This is based on Klemsa/Zama partial key ideas.

Look at:

- Types: `SubsetKeySwitchingKey<P>` in `include/params.hpp`
- Keygen: `subikskgen` / `subprivkskgen` in `include/evalkeygens.hpp`

## Annihilate key switching / Modified Chen packing

`AnnihilateCircuitBootstrapping` and `AnnihilateKeySwitching` implement the
packing technique described in the Modified Chen’s packing paper. Enable with
the default 128‑bit parameters; generate:

```cpp
ek.emplaceahk<AHlvl2param>(sk);
ek.emplacecbsk<AHlvl2param>(sk);
```

See `include/circuitbootstrapping.hpp` and `test/annihilate*.cpp`.

## Double Decomposition (Bivariate Representation)

Double Decomposition is a technique from the [TFHE Double Decomposition paper](https://eprint.iacr.org/2023/771) that uses bivariate polynomial representation for more flexible parameter choices in TFHE external products.

### Overview

In standard TFHE, the external product decomposes input coefficients into `l` digits in base `Bg`. Double Decomposition adds an auxiliary decomposition level, decomposing each digit further into `l̅` sub-digits in base `B̅g`. This creates a bivariate structure:

```
a ≈ Σᵢ Σⱼ aᵢⱼ * Bg^(l-i) * B̅g^(l̅-j)
```

When `l̅=1`, this reduces to standard decomposition.

### Parameters

Double Decomposition introduces additional parameters in the `*param` structs:

| Parameter | Description |
|-----------|-------------|
| `l` | Primary decomposition depth |
| `Bgbit` | Primary base bits (Bg = 2^Bgbit) |
| `l̅` | Auxiliary decomposition depth |
| `B̅gbit` | Auxiliary base bits (B̅g = 2^B̅gbit) |
| `lₐ` | Primary decomposition depth (nonce/key part) |
| `Bgₐbit` | Primary base bits (nonce part) |
| `l̅ₐ` | Auxiliary decomposition depth (nonce part) |
| `B̅gₐbit` | Auxiliary base bits (nonce part) |

The total decomposition bits should satisfy: `l * Bgbit + l̅ * B̅gbit ≤ width` where `width` is the Torus bit width (32, 64, or 128).

### 128-bit Torus Support

TFHEpp supports 128-bit Torus (`__uint128_t`) for Double Decomposition via `lvl3param`. The FFT backend uses 64-bit operations internally, which works because DD always operates on decomposition digits (small integers), not raw Torus values.

See `include/params/128bit.hpp` for the `lvl3param` definition.

### Usage

Double Decomposition is automatically used when `l̅ > 1` or `l̅ₐ > 1` in your parameter set. The `ExternalProduct` function detects this and uses the DD path:

```cpp
// Standard parameters (l̅=1): uses standard decomposition
// DD parameters (l̅>1): automatically uses Double Decomposition
TFHEpp::ExternalProduct<P>(result, trlwe, trgswfft);
```

### TRGSW Encryption with DD

When using DD, TRGSW encryption follows this process:
1. Create an ordinary TRGSW with `k*lₐ + l` rows
2. Apply auxiliary decomposition to each row, expanding to `k*lₐ*l̅ₐ + l*l̅` rows

This is handled automatically by `trgswSymEncrypt` when DD parameters are detected.

### Key Functions

Located in `include/trgsw.hpp`:

- `DecompositionImpl<P, IsNonce, AuxOnly>` — Unified decomposition (standard and DD)
- `TRLWEBaseBbarDecompose` — Decompose TRLWE to base B̅g
- `RecombineTRLWEFromDD` — Recombine l̅ TRLWEs back to single TRLWE
- `ExternalProduct` — Automatically uses DD when `l̅ > 1`

### Tests

- `test/externalproductdoubledecomposition.cpp` — Tests both standard (l̅=1) and DD (l̅=2) paths
- `test/gatebootstrappingtlwe2tlwedoubledecomposition.cpp` — 128-bit gate bootstrapping with DD

## BFV++ (experimental)

TFHEpp includes partial BFV/BFV++ support in `include/bfv++.hpp` for users who
want ring arithmetic alongside TFHE. This is not as feature‑complete as the TFHE
stack; treat it as experimental.

## Performance tips

- Prefer FFT backends matching your hardware (SPQLIOS for AVX2, AVX512 SPQLIOS,
  AArch64 backend, or MKL/FFTW3).
- Generate only the evaluation keys you need; some (BK/CB) are large.
- Use batch encrypt/decrypt APIs which are OpenMP‑parallelized.

