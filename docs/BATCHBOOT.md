# BatchBoot

TFHEpp implements the core algorithms from “BatchBoot: Fast Batched
Bootstrapping for TFHE Scheme and Practical Applications” (USENIX Security
2026):

- radix-4 efficient monomial–polynomial multiplication (EMPmul), including
  parametrized external products and hoisted Fourier transforms;
- full and sparse functional BatchBoot, including RLWE-to-MLWE extraction,
  repacking, homomorphic trace, and optional RLWE key switching;
- the Half.BatchCBS conversion and the external-product tree used by batched
  circuit bootstrapping.

The implementation is in `include/tfhe/batchboot.hpp`. It is independently
written against TFHEpp's native ciphertext and FFT interfaces from the
published algorithms. No source file from the reference artifact is included
or copied into TFHEpp. The paper and artifact are listed for provenance:

- <https://www.usenix.org/conference/usenixsecurity26/presentation/li-zhihao>
- <https://doi.org/10.5281/zenodo.17936945>

## Functional BatchBoot

For sparse packing, first extract the active coefficient coset:

```cpp
constexpr std::uint32_t slots = 8;
using RingP = TFHEpp::lvl1param;
using ModuleP = TFHEpp::BatchRingSwitchP<RingP, slots>;

auto module_key =
    TFHEpp::BatchRingSwitchSecret<RingP, slots>(sparse_ring_key);
TFHEpp::TRLWE<ModuleP> module_ciphertext;
TFHEpp::BatchRingSwitch<RingP, slots>(module_ciphertext, ring_ciphertext);
```

Generate the BatchBoot and packing keys under the accumulator key, build a
test polynomial, and bootstrap:

```cpp
TFHEpp::BatchBootKey<ModuleP, RingP> bsk;
TFHEpp::BatchBootKeyGen(bsk, module_key, accumulator_key);

TFHEpp::AnnihilateKey<RingP> automorphism_keys;
TFHEpp::annihilatekeygen(automorphism_keys, accumulator_key);

auto lut = TFHEpp::MakeBatchBootTestVector<RingP, 3>(
    [](std::uint32_t x) { return x; });
TFHEpp::TRLWE<RingP> output;
TFHEpp::BatchBootstrap(output, module_ciphertext, lut.polynomial, bsk,
                       automorphism_keys, lut.exponent_bias);
```

As in the paper, one input bit is reserved for anti-cyclic padding. A test
vector with `InputBits` therefore accepts messages in
`[0, 2^(InputBits-1))`. `MakeBatchBootTestVector` places those messages at the
centers of the LUT intervals so modulus-switch rounding does not put an exact
encoding on a table boundary. Repacked results occupy coefficients
`i * RingP::n / slots`.

## Circuit BatchBoot

`BatchHalfCircuitBootstrap` and `BatchHalfCircuitBootstrapFFT` implement the
first half of Algorithm 4 and return one RGSW monomial control per active
coefficient. `BatchExternalProductTree` implements Algorithm 8. For a
single-output LUT, `BatchCircuitBootstrap` composes both stages; callers with
multiple output functions can generate the controls once and reuse them in
one product tree per output.

## Constraints

- The input secret for `BatchBootKeyGen` must be sparse binary.
- The serialized/evaluator-visible key shape includes each module component's
  Hamming weight. Use a fixed public weight profile if that distribution must
  not depend on secret key generation.
- The active slot count must be a power of two and no greater than the target
  polynomial degree.
- The Algorithm 8 radix must be a power of two in `[2, N]`.
- EMPmul currently uses TFHEpp parameter sets with standard decomposition
  (`l̅ = l̅ₐ = 1`), including the default level-1 and level-2 parameters.
- `BatchHalfCircuitBootstrap` currently supports `k = 1`, equal nonce/body
  gadget lengths, and a power-of-two gadget length. The default level-2
  parameter has the required four-row gadget.
- Keys are intentionally move-only and should normally be heap allocated;
  realistic BatchBoot keys are large.

The implementation validates structural preconditions but does not select
cryptographic parameters. Production deployments must use a parameter set and
failure analysis appropriate to their message precision, sparse weight, and
number of EMP invocations.
