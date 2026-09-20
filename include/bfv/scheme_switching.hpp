#pragma once

#include <array>
#include <memory>

#include "bfv++.hpp"
#include "bfv-slots.hpp"
#include "tfhe/evalkeygens.hpp"
#include "tfhe/gatebootstrapping.hpp"
#include "tfhe/homdecomp.hpp"

#ifdef TFHEPP_DEFAULT_128BIT_PARAMS
namespace TFHEpp {

// Scalar BFV: an operand occupies only the constant coefficient. Wider
// profiles are experimental: representability does not imply enough noise
// margin. The checked-in 8-bit bound must not be applied to them unchanged.
template <unsigned InputBits>
struct BFVSwitchParamFor : lvl3param {
    static_assert(InputBits >= 1 && InputBits <= 32);
    static constexpr unsigned input_bits = InputBits;
    static constexpr uint32_t plain_modulusbit = 2 * InputBits;
    static constexpr T plain_modulus = T{1} << plain_modulusbit;
    static constexpr T delta_int = T{1} << (128 - plain_modulusbit);
    static constexpr uint64_t Q_mod_t = 0;
    static constexpr double Δ = static_cast<double>(delta_int);
};

// Separate 16x16->32 profile. The larger ring permits smaller relative noise
// without weakening the generic LWE screen. The primary gadget keeps 120 bits;
// full 8x16-bit DD coverage and the 128-bit ciphertext modulus are unchanged.
// See docs/BFVSchemeSwitching.md for the conditional bound and measured cost.
template <>
struct BFVSwitchParamFor<16> : lvl3param {
    static constexpr uint32_t nbit = 13;
    static constexpr uint32_t n = 1U << nbit;
    static constexpr uint32_t l = 12, lₐ = 12;
    static constexpr uint32_t Bgbit = 10, Bgₐbit = 10;
    static constexpr uint32_t Bg = 1U << Bgbit, Bgₐ = 1U << Bgₐbit;
    static inline const double α = std::ldexp(1.0, -124);
    static constexpr unsigned input_bits = 16;
    static constexpr uint32_t plain_modulusbit = 32;
    static constexpr T plain_modulus = T{1} << plain_modulusbit;
    static constexpr T delta_int = T{1} << 96;
    static constexpr uint64_t Q_mod_t = 0;
    static constexpr double Δ = static_cast<double>(delta_int);
};

template <unsigned InputBits>
struct TFHEToBFVBootstrapParamFor {
    using domainP = lvlhalfparam;
    using targetP = BFVSwitchParamFor<InputBits>;
#ifdef USE_KEY_BUNDLE
    static constexpr uint32_t Addends = 2;
#else
    static constexpr uint32_t Addends = 1;
#endif
};

template <unsigned InputBits>
struct BFVToTFHEKeySwitchParamFor : lvl31param {
    using domainP = BFVSwitchParamFor<InputBits>;
    static constexpr uint32_t t = InputBits == 16 ? 6 : lvl31param::t;
    static constexpr uint32_t basebit = InputBits == 16 ? 3 : lvl31param::basebit;
};

using BFVSwitchParam = BFVSwitchParamFor<8>;
using TFHEToBFVBootstrapParam = TFHEToBFVBootstrapParamFor<8>;
using BFVToTFHEKeySwitchParam = BFVToTFHEKeySwitchParamFor<8>;
using BFVSwitchInput = std::array<TLWE<lvl1param>, 8>;
using BFVSwitchOutput = std::array<TLWE<lvl1param>, 16>;

// No secret key is retained by the evaluator. All keys use the explicitly
// supplied independent secret chain, including the half-level bridge.
template <unsigned InputBits>
struct TFHEBFVSwitchKeyFor {
    std::unique_ptr<KeySwitchingKey<lvl1hparam>> to_half;
    std::unique_ptr<BootstrappingKeyFFT<TFHEToBFVBootstrapParamFor<InputBits>>> to_bfv;
    std::unique_ptr<AnnihilateKey<BFVSwitchParamFor<InputBits>>> packing;
    std::unique_ptr<relinKeyFFT<BFVSwitchParamFor<InputBits>>> relin;
    std::unique_ptr<KeySwitchingKey<BFVToTFHEKeySwitchParamFor<InputBits>>> to_tfhe;
    std::unique_ptr<BootstrappingKeyFFT<lvlh1param>> boolean_bootstrap;
};
using TFHEBFVSwitchKey = TFHEBFVSwitchKeyFor<8>;

template <unsigned InputBits>
inline void TFHEBFVSwitchKeyGen(TFHEBFVSwitchKeyFor<InputBits> &out,
                              const Key<lvl1param> &io_key,
                              const Key<lvlhalfparam> &half_key,
                              const Key<BFVSwitchParamFor<InputBits>> &bfv_key)
{
    using BFVSwitchParam = BFVSwitchParamFor<InputBits>;
    using TFHEToBFVBootstrapParam = TFHEToBFVBootstrapParamFor<InputBits>;
    using BFVToTFHEKeySwitchParam = BFVToTFHEKeySwitchParamFor<InputBits>;
    out.to_half = std::make_unique<KeySwitchingKey<lvl1hparam>>();
    ikskgen<lvl1hparam>(*out.to_half, io_key, half_key);
    out.to_bfv =
        std::make_unique<BootstrappingKeyFFT<TFHEToBFVBootstrapParam>>();
#ifdef USE_KEY_BUNDLE
    bkfftgen<TFHEToBFVBootstrapParam>(*out.to_bfv, half_key, bfv_key);
#else
    // The half-level secret is binary; key rows are independent, and the
    // library's CSPRNG and FFT workspaces are thread-local.
#pragma omp parallel for schedule(static)
    for (int i = 0; i < lvlhalfparam::n; ++i) {
        Polynomial<BFVSwitchParam> plain{};
        plain[0] = half_key[i];
        trgswSymEncrypt<BFVSwitchParam>((*out.to_bfv)[i][0], plain, bfv_key);
    }
#endif
    out.packing = std::make_unique<AnnihilateKey<BFVSwitchParam>>();
    annihilatekeygen<BFVSwitchParam>(*out.packing, bfv_key);
    out.relin = makeRelinKeyFFT<BFVSwitchParam>(bfv_key);
    out.to_tfhe =
        std::make_unique<KeySwitchingKey<BFVToTFHEKeySwitchParam>>();
    ikskgen<BFVToTFHEKeySwitchParam>(*out.to_tfhe, bfv_key, io_key);
    out.boolean_bootstrap = std::make_unique<BootstrappingKeyFFT<lvlh1param>>();
    bkfftgen<lvlh1param>(*out.boolean_bootstrap, half_key, io_key);
}

template <unsigned InputBits>
inline void TFHEToBFV(TRLWE<BFVSwitchParamFor<InputBits>> &out,
                      const std::array<TLWE<lvl1param>, BFVSwitchParamFor<InputBits>::input_bits> &bits,
                      const TFHEBFVSwitchKeyFor<InputBits> &key)
{
    using P = BFVSwitchParamFor<InputBits>;
    constexpr unsigned digit_bits = InputBits == 16 ? 2 : 1;
    auto weighted = std::make_unique<std::array<TLWE<P>, InputBits / digit_bits>>();
    if constexpr (digit_bits == 2) {
        // Re-encode signs at Q/16 before combining. Dividing an existing
        // torus ciphertext by two would not preserve its encrypted phase.
        std::array<TLWE<lvl1param>, InputBits> normalized{};
        Polynomial<lvl1param> normalize_lut;
        normalize_lut.fill(uint32_t{1} << 28);
#pragma omp parallel for schedule(static)
        for (int bit = 0; bit < InputBits; ++bit) {
            TLWE<lvlhalfparam> small{};
            IdentityKeySwitch<lvl1hparam>(small, bits[bit], *key.to_half);
            GateBootstrappingTLWE2TLWE<lvlh1param>(
                normalized[bit], small, *key.boolean_bootstrap, normalize_lut);
        }
#pragma omp parallel for schedule(static)
        for (int digit = 0; digit < InputBits / 2; ++digit) {
            TLWE<lvl1param> packed{};
            for (size_t i = 0; i < packed.size(); ++i)
                packed[i] = normalized[2 * digit][i] +
                            2 * normalized[2 * digit + 1][i];
            TLWE<lvlhalfparam> small{};
            IdentityKeySwitch<lvl1hparam>(small, packed, *key.to_half);
            const typename P::T half_weight = (P::delta_int / 2) << (2 * digit);
            Polynomial<P> lut;
            // Anti-periodic extension maps phases (-3,-1,1,3)/16 to
            // (-3,-1,1,3)*half_weight, with address radius 1/16.
            for (size_t i = 0; i < P::n; ++i)
                lut[i] = (i < P::n / 4 || i >= 3 * P::n / 4)
                             ? half_weight : 3 * half_weight;
            GateBootstrappingTLWE2TLWE<TFHEToBFVBootstrapParamFor<InputBits>>(
                (*weighted)[digit], small, *key.to_bfv, lut);
            (*weighted)[digit][P::n] += 3 * half_weight;
        }
    }
    else {
#pragma omp parallel for schedule(static)
    for (int bit = 0; bit < InputBits; ++bit) {
        TLWE<lvlhalfparam> small{};
        IdentityKeySwitch<lvl1hparam>(small, bits[bit], *key.to_half);
        // PBS returns +/- half the desired weight; the public offset makes
        // the two outputs 0 and 2^bit * Delta without scaling PBS noise.
        const typename P::T half_weight = (P::delta_int / 2) << bit;
        Polynomial<P> lut;
        lut.fill(half_weight);
        GateBootstrappingTLWE2TLWE<TFHEToBFVBootstrapParamFor<InputBits>>(
            (*weighted)[bit], small, *key.to_bfv, lut);
        (*weighted)[bit][P::n] += half_weight;
    }
    }
    TLWE<P> sum{};
    for (const auto &bit : *weighted)
        for (size_t i = 0; i < bit.size(); ++i) sum[i] += bit[i];
    // A single occupied coefficient needs only the iterative trace, avoiding
    // the large recursive stack frames of the general multi-LWE packer.
    InvSampleExtractIndex<P>(out, sum, 0);
    AnnihilateKeySwitching<P>(out, out, *key.packing);
}

template <unsigned InputBits>
inline void BFVToTFHE(std::array<TLWE<lvl1param>, 2 * InputBits> &out,
                      const TRLWE<BFVSwitchParamFor<InputBits>> &in,
                      const TFHEBFVSwitchKeyFor<InputBits> &key)
{
    using P = BFVSwitchParamFor<InputBits>;
    TLWE<P> extracted{};
    SampleExtractIndex<P>(extracted, in, 0);
    // The window has two low guard bits, followed by 2*InputBits data bits.
    // Bias by half a guard unit. HomDecomp removes the guard digit before
    // Boolean extraction; directly key-switching a 16-bit phase would lose
    // precision in the small TFHE bootstrap domain.
    extracted[P::n] += P::delta_int / 8;
    std::array<TLWE<lvl1param>, InputBits + 1> digits{};
    HomDecomp<BFVToTFHEKeySwitchParamFor<InputBits>, lvl1hparam, lvlh1param,
              2, InputBits + 1, 2 * InputBits + 2>(
        digits, extracted, *key.to_tfhe, *key.to_half,
        *key.boolean_bootstrap);
    Polynomial<lvl1param> boolean_lut;
    boolean_lut.fill(-lvl1param::μ);
#pragma omp parallel for schedule(static)
    for (int digit = 1; digit <= InputBits; ++digit) {
        digits[digit][lvl1param::n] += uint32_t{1} << 29;
        TLWE<lvlhalfparam> small{}, low{};
        IdentityKeySwitch<lvl1hparam>(small, digits[digit], *key.to_half);
        for (size_t i = 0; i < small.size(); ++i) low[i] = small[i] << 1;
        GateBootstrappingTLWE2TLWE<lvlh1param>(
            out[2 * (digit - 1)], low, *key.boolean_bootstrap, boolean_lut);
        GateBootstrappingTLWE2TLWE<lvlh1param>(
            out[2 * (digit - 1) + 1], small, *key.boolean_bootstrap,
            boolean_lut);
    }
}

}  // namespace TFHEpp
#endif  // TFHEPP_DEFAULT_128BIT_PARAMS
