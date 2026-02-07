# Parameters and Security Levels

TFHEpp chooses a parameter set at compile time. The selection happens in
`include/params.hpp` by checking build‑time macros set by CMake.

## Built‑in parameter sets

Located under `include/params/`:

- `128bit.hpp` — default 128‑bit security parameters.
- `concrete.hpp` — fast 128‑bit parameters inspired by CONCRETE.
- `CGGI19.hpp` — parameters from CGGI19.
- `CGGI16.hpp` — legacy 80‑bit parameters (selected by `USE_80BIT_SECURITY`).
- `compress.hpp` — compressed ciphertext parameters.
- `tfhe-rs.hpp` — TFHE‑rs‑style parameters.
- `ternary.hpp` — ternary‑key parameters.

Select one with CMake, for example:

```bash
cmake .. -DUSE_CONCRETE=ON
```

Only one of these options should be enabled at a time.

## Key distribution

Many level‑1+ parameter sets use ternary secrets (`{-1,0,1}`) for better noise
growth. This is reflected in each `*param` struct in `include/params/*.hpp`.

## Double Decomposition parameters

For Double Decomposition (bivariate representation), parameter sets include
auxiliary decomposition parameters:

| Parameter | Description |
|-----------|-------------|
| `l̅` | Auxiliary decomposition depth (default: 1 = standard decomposition) |
| `B̅gbit` | Auxiliary base bits |
| `l̅ₐ` | Auxiliary depth for nonce/key part |
| `B̅gₐbit` | Auxiliary base bits for nonce part |



## Notation correspondence

TFHEpp parameter structs (in `include/params/*.hpp`) follow the TFHE notation with a few
implementation-oriented splits:

| Concept | TFHEpp field | Meaning | Common paper notation |
|--------|--------------|---------|------------------------|
| Ring degree | `nbit`, `n` | `n = 2^nbit` polynomial degree | `N` (often `N=2^nbit`) |
| TLWE/TRLWE dimension | `k` | number of polynomials minus 1 (`k=1` for TFHE-style) | `k` |
| Primary (TRGSW) decomposition levels | `l`, `lₐ` | #digits for gadget decomposition (main / nonce part) | `ℓ` / `ℓ_a` |
| Primary base (bits) | `Bgbit`, `Bgₐbit` | base is `Bg = 2^{Bgbit}` (main / nonce) | `B_g` bits / `B_g` |
| Primary base (value) | `Bg`, `Bgₐ` | `2^{Bgbit}` (main / nonce) | `B_g` |
| DD limb count | `l̅`, `l̅ₐ` | #limbs for bivariate/DD row representation (main / nonce) | `\\bar{ℓ}` / `\\bar{ℓ}_a` |
| DD limb size (bits) | `B̅gbit`, `B̅gₐbit` | limb base is `2^{B̅gbit}` (main / nonce) | `K` (limb size in bits) |
| Fresh noise | `α` (or `η`) | error distribution parameter (depends on `errordist`) | `α` / `σ` |
| Torus type | `T` | torus word type (e.g. `uint64_t`, `__uint128_t`) | `\mathbb{T}` discretization |

Notes:
- In TFHEpp's current DD external product, the *input* TRLWE is still decomposed using the
  primary `(l, Bgbit)`; DD parameters `(l̅, B̅gbit)` control how TRGSW rows are stored/recombined.
- For FFT-based polynomial multiplication with `double`, a practical safety constraint is
  `Bgbit + B̅gbit + nbit < 53` (often keep a few extra margin bits for repeated accumulations).

When `l̅=1`, standard decomposition is used. When `l̅>1`, Double Decomposition
is automatically enabled. See `docs/Advanced.md` for details.

The 128-bit parameter set (`lvl3param` in `include/params/128bit.hpp`) uses
`__uint128_t` Torus and is designed for Double Decomposition.

## Custom parameters

To add your own parameter set:

1. Create `include/params/myparams.hpp` defining `lvl0param`, `lvl1param`, …
2. Add a macro guard and include in `include/params.hpp` (mirroring the
   existing pattern).
3. Build with `-DUSE_MY_PARAMS=ON` (or add a new CMake option).

When changing parameters, remember that evaluation keys (bootstrapping,
key‑switching, circuit‑bootstrapping keys) must be regenerated.

