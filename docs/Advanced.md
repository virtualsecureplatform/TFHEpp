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

## BFV++ (experimental)

TFHEpp includes partial BFV/BFV++ support in `include/bfv++.hpp` for users who
want ring arithmetic alongside TFHE. This is not as feature‑complete as the TFHE
stack; treat it as experimental.

## Performance tips

- Prefer FFT backends matching your hardware (SPQLIOS for AVX2, AVX512 SPQLIOS,
  AArch64 backend, or MKL/FFTW3).
- Generate only the evaluation keys you need; some (BK/CB) are large.
- Use batch encrypt/decrypt APIs which are OpenMP‑parallelized.

