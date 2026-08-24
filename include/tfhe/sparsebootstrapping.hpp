#pragma once

#include <bit>
#include <cstdint>
#include <limits>
#include <random>
#include <span>
#include <stdexcept>

#include "gatebootstrapping.hpp"

namespace TFHEpp {

// Encrypt a Boolean TLWE whose mask and phase are explicitly represented in
// Z_q and embedded in the surrounding power-of-two torus.  The regular
// tlweSymEncrypt samples a full torus mask, which is correct for TFHEpp's
// default parameters but does not instantiate the q = 512 LWE layer in the
// shallow-bootstrap paper.
template <class P>
void tlweSymEncryptModQ(TLWE<P> &res, const bool message,
                        const Key<P> &key, const double alpha = P::α)
{
    static_assert(requires { P::q; P::qbit; },
                  "explicit-q encryption requires q and qbit parameters");
    static_assert(std::has_single_bit(P::q),
                  "explicit-q encryption requires a power-of-two modulus");
    static_assert(P::qbit <= std::numeric_limits<typename P::T>::digits,
                  "explicit-q modulus does not fit the torus type");

    using T = typename P::T;
    constexpr uint32_t shift = std::numeric_limits<T>::digits - P::qbit;
    std::uniform_int_distribution<T> maskdist(0, static_cast<T>(P::q - 1));
    res = {};
    res[P::k * P::n] =
        ModularGaussian<P>(message ? static_cast<T>(P::μ)
                           : static_cast<T>(-P::μ),
                           alpha);
    for (uint32_t i = 0; i < P::k * P::n; i++) {
        res[i] = static_cast<T>(maskdist(generator) << shift);
        res[P::k * P::n] += res[i] * key[i];
    }
}

// LUT for a binary value embedded as +/- q/8 when the blind-rotation modulus
// q is smaller than 2N.  The ordinary constant TFHE test vector places its
// sign boundary at N and therefore cannot distinguish values in Z_512 inside
// an N = 1024 accumulator.
template <class P, typename P::targetP::T mu>
Polynomial<typename P::targetP> ShallowBootBinaryIdentityTestVector()
{
    constexpr uint32_t modulus = BlindRotationModulus<P>;
    static_assert(modulus <= P::targetP::n,
                  "small-modulus shallow LUT requires q <= N");
    constexpr uint32_t exponent_scale = 2 * P::targetP::n / modulus;
    Polynomial<typename P::targetP> result = {};
    for (uint32_t exponent = 0; exponent < modulus; exponent++) {
        const uint32_t ring_exponent = exponent * exponent_scale;
        const auto desired = exponent < modulus / 2
                                 ? static_cast<typename P::targetP::T>(-mu)
                                 : static_cast<typename P::targetP::T>(mu);
        if (ring_exponent == 0)
            result[0] = desired;
        else if (ring_exponent < P::targetP::n)
            // X^e * X^(N-e) = -1 in Z[X]/(X^N + 1).
            result[P::targetP::n - ring_exponent] = -desired;
        else if (ring_exponent == P::targetP::n)
            result[0] = -desired;
        else
            result[2 * P::targetP::n - ring_exponent] = desired;
    }
    return result;
}

// Algorithm 2 of Jain--Lin--Liu--Saha (ePrint 2026/1730), specialized to a
// structured sparse binary LWE secret.  `block_offsets` partitions the LWE
// secret coordinates.  The secret must contain exactly one 1 in every block;
// this lets all the RGSWs in a block be added before a single external product.
inline void ValidateStructuredSparseBlocks(
    const std::span<const std::uint32_t> block_offsets,
    const std::uint32_t dimension)
{
    if (block_offsets.size() < 2 || block_offsets.front() != 0 ||
        block_offsets.back() != dimension)
        throw std::invalid_argument(
            "sparse bootstrap blocks must partition the LWE coordinates");
    for (std::size_t block = 0; block + 1 < block_offsets.size(); block++)
        if (block_offsets[block] >= block_offsets[block + 1])
            throw std::invalid_argument(
                "sparse bootstrap blocks must be nonempty and ordered");
}

// Sample the structured sparse binary secret required by
// StructuredSparseBlindRotate.  Keeping this separate from SecretKey makes
// the altered LWE secret distribution an explicit opt-in.
template <class P, class URBG>
void structuredSparseBinaryKeyGen(
    Key<P> &key, const std::span<const std::uint32_t> block_offsets,
    URBG &engine)
{
    static_assert(P::k == 1,
                  "structured sparse key generation requires an LWE key");
    static_assert(P::key_value_min == 0 && P::key_value_max == 1,
                  "structured sparse key generation requires binary keys");
    ValidateStructuredSparseBlocks(block_offsets, P::n);

    key.fill(typename P::T{0});
    for (std::size_t block = 0; block + 1 < block_offsets.size(); block++) {
        std::uniform_int_distribution<std::uint32_t> choose(
            block_offsets[block], block_offsets[block + 1] - 1);
        key[choose(engine)] = typename P::T{1};
    }
}

template <class P>
void StructuredSparseBlindRotate(
    TRLWE<typename P::targetP> &res, const TLWE<typename P::domainP> &tlwe,
    const BootstrappingKeyFFT<P> &bkfft,
    const Polynomial<typename P::targetP> &testvector,
    const std::span<const std::uint32_t> block_offsets)
{
    static_assert(P::Addends == 1,
                  "structured sparse blind rotation is incompatible with key "
                  "bundling");
    static_assert(P::domainP::k == 1,
                  "structured sparse blind rotation requires an LWE domain");
    static_assert(P::domainP::key_value_min == 0 &&
                      P::domainP::key_value_max == 1,
                  "structured sparse blind rotation requires binary keys");

    ValidateStructuredSparseBlocks(block_offsets, P::domainP::n);

    ModswitchTLWE<typename P::domainP> moded;
    BRModSwitch<P, 1>(moded, tlwe);
    res = {};
    PolynomialMulByXai<typename P::targetP>(
        res[P::targetP::k], testvector, moded[P::domainP::n]);

    for (std::size_t block = 0; block + 1 < block_offsets.size(); block++) {
        const std::uint32_t first = block_offsets[block];
        const std::uint32_t last = block_offsets[block + 1];
        CMUXFFTwithScheduledPolynomialMulByXaiMinusOne<P>(
            res, std::span<const BootstrappingKeyElementFFT<P>>(
                     bkfft.data() + first, last - first),
            std::span<const typename P::domainP::T>(moded.data() + first,
                                                    last - first));
    }
}

template <class P>
void StructuredSparseGateBootstrappingTLWE2TLWE(
    TLWE<typename P::targetP> &res, const TLWE<typename P::domainP> &tlwe,
    const BootstrappingKeyFFT<P> &bkfft,
    const Polynomial<typename P::targetP> &testvector,
    const std::span<const std::uint32_t> block_offsets)
{
    alignas(64) TRLWE<typename P::targetP> acc;
    StructuredSparseBlindRotate<P>(acc, tlwe, bkfft, testvector,
                                   block_offsets);
    SampleExtractIndex<typename P::targetP>(res, acc, 0);
}

}  // namespace TFHEpp
