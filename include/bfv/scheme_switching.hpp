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

// Scalar BFV over Z_(2^16)[X]/(X^4096+1). This is deliberately not a
// SIMD parameter: an operand occupies only the constant coefficient.
struct BFVSwitchParam : lvl3param {
    static constexpr uint32_t plain_modulusbit = 16;
    static constexpr T plain_modulus = T{1} << plain_modulusbit;
    static constexpr T delta_int = T{1} << 112;
    static constexpr uint64_t Q_mod_t = 0;
    static constexpr double Δ = static_cast<double>(delta_int);
};

struct TFHEToBFVBootstrapParam {
    using domainP = lvlhalfparam;
    using targetP = BFVSwitchParam;
#ifdef USE_KEY_BUNDLE
    static constexpr uint32_t Addends = 2;
#else
    static constexpr uint32_t Addends = 1;
#endif
};

struct BFVToTFHEKeySwitchParam : lvl31param {
    using domainP = BFVSwitchParam;
};

using BFVSwitchInput = std::array<TLWE<lvl1param>, 8>;
using BFVSwitchOutput = std::array<TLWE<lvl1param>, 16>;

// No secret key is retained by the evaluator. All keys use the explicitly
// supplied independent secret chain, including the half-level bridge.
struct TFHEBFVSwitchKey {
    std::unique_ptr<KeySwitchingKey<lvl1hparam>> to_half;
    std::unique_ptr<BootstrappingKeyFFT<TFHEToBFVBootstrapParam>> to_bfv;
    std::unique_ptr<AnnihilateKey<BFVSwitchParam>> packing;
    std::unique_ptr<relinKeyFFT<BFVSwitchParam>> relin;
    std::unique_ptr<KeySwitchingKey<BFVToTFHEKeySwitchParam>> to_tfhe;
    std::unique_ptr<BootstrappingKeyFFT<lvlh1param>> boolean_bootstrap;
};

inline void TFHEBFVSwitchKeyGen(TFHEBFVSwitchKey &out,
                              const Key<lvl1param> &io_key,
                              const Key<lvlhalfparam> &half_key,
                              const Key<BFVSwitchParam> &bfv_key)
{
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

inline void TFHEToBFV(TRLWE<BFVSwitchParam> &out,
                      const BFVSwitchInput &bits,
                      const TFHEBFVSwitchKey &key)
{
    using P = BFVSwitchParam;
    auto weighted = std::make_unique<std::array<TLWE<P>, 8>>();
#pragma omp parallel for schedule(static)
    for (int bit = 0; bit < 8; ++bit) {
        TLWE<lvlhalfparam> small{};
        IdentityKeySwitch<lvl1hparam>(small, bits[bit], *key.to_half);
        // PBS returns +/- half the desired weight; the public offset makes
        // the two outputs 0 and 2^bit * Delta without scaling PBS noise.
        const P::T half_weight = (P::delta_int / 2) << bit;
        Polynomial<P> lut;
        lut.fill(half_weight);
        GateBootstrappingTLWE2TLWE<TFHEToBFVBootstrapParam>(
            (*weighted)[bit], small, *key.to_bfv, lut);
        (*weighted)[bit][P::n] += half_weight;
    }
    TLWE<P> sum{};
    for (const auto &bit : *weighted)
        for (size_t i = 0; i < bit.size(); ++i) sum[i] += bit[i];
    // A single occupied coefficient needs only the iterative trace, avoiding
    // the large recursive stack frames of the general multi-LWE packer.
    InvSampleExtractIndex<P>(out, sum, 0);
    AnnihilateKeySwitching<P>(out, out, *key.packing);
}

inline void BFVToTFHE(BFVSwitchOutput &out,
                      const TRLWE<BFVSwitchParam> &in,
                      const TFHEBFVSwitchKey &key)
{
    using P = BFVSwitchParam;
    TLWE<P> extracted{};
    SampleExtractIndex<P>(extracted, in, 0);
    // The 18-bit window has two low guard bits, followed by 16 data bits.
    // Bias by half a guard unit. HomDecomp removes the guard digit before
    // Boolean extraction; directly key-switching a 16-bit phase would lose
    // precision in the small TFHE bootstrap domain.
    extracted[P::n] += P::delta_int / 8;
    std::array<TLWE<lvl1param>, 9> digits{};
    HomDecomp<BFVToTFHEKeySwitchParam, lvl1hparam, lvlh1param, 2, 9, 18>(
        digits, extracted, *key.to_tfhe, *key.to_half,
        *key.boolean_bootstrap);
    Polynomial<lvl1param> boolean_lut;
    boolean_lut.fill(-lvl1param::μ);
#pragma omp parallel for schedule(static)
    for (int digit = 1; digit <= 8; ++digit) {
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
