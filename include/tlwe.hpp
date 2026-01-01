#pragma once

#include <array>
#include <cstdint>
#include <vector>

#include "key.hpp"
#include "params.hpp"
#include "utils.hpp"

namespace TFHEpp {
template <class P>
void tlweSymEncrypt(TLWE<P> &res, const typename P::T p, const double α,
                    const Key<P> &key)
{
    res = {};
    res[P::k * P::n] = ModularGaussian<P>(p, α);
    for (int k = 0; k < P::k; k++)
        for (int i = 0; i < P::n; i++) {
            res[k * P::n + i] = UniformTorusRandom<P>();
            res[P::k * P::n] += res[k * P::n + i] * key[k * P::n + i];
        }
}

template <class P>
void tlweSymEncrypt(TLWE<P> &res, const typename P::T p, const uint η,
                    const Key<P> &key)
{
    std::uniform_int_distribution<typename P::T> Torusdist(0, P::q - 1);
    res = {};
    res[P::k * P::n] =
        p + CenteredBinomial<P>(η)
        << (std::numeric_limits<typename P::T>::digits - P::qbit);
    for (int k = 0; k < P::k; k++)
        for (int i = 0; i < P::n; i++) {
            res[k * P::n + i] =
                Torusdist(generator)
                << (std::numeric_limits<typename P::T>::digits - P::qbit);
            res[P::k * P::n] += res[k * P::n + i] * key[k * P::n + i];
        }
}

template <class P>
void tlweSymEncrypt(TLWE<P> &res, const typename P::T p, const Key<P> &key)
{
    if constexpr (P::errordist == ErrorDistribution::ModularGaussian)
        tlweSymEncrypt<P>(res, p, P::α, key);
    else
        tlweSymEncrypt<P>(res, p, P::η, key);
}

template <class P>
void tlweSymEncrypt(TLWE<P> &res, const typename P::T p, const SecretKey &sk)
{
    tlweSymEncrypt<P>(res, p, sk.key.get<P>());
}

template <class P, uint plain_modulus = P::plain_modulus>
void tlweSymIntEncrypt(TLWE<P> &res, const typename P::T p, const double α,
                       const Key<P> &key)
{
    const double Δ = std::pow(2.0, std::numeric_limits<typename P::T>::digits) /
                     plain_modulus;
    tlweSymEncrypt<P>(res, static_cast<typename P::T>(p * Δ), α, key);
}

template <class P, uint plain_modulus = P::plain_modulus>
void tlweSymIntEncrypt(TLWE<P> &res, const typename P::T p, const uint η,
                       const Key<P> &key)
{
    constexpr double Δ =
        std::pow(2.0, std::numeric_limits<typename P::T>::digits) /
        plain_modulus;
    tlweSymEncrypt<P>(res, static_cast<typename P::T>(p * Δ), η, key);
}

template <class P, uint plain_modulus = P::plain_modulus>
void tlweSymIntEncrypt(TLWE<P> &res, const typename P::T p, const Key<P> &key)
{
    if constexpr (P::errordist == ErrorDistribution::ModularGaussian)
        tlweSymIntEncrypt<P, plain_modulus>(res, p, P::α, key);
    else
        tlweSymIntEncrypt<P, plain_modulus>(res, p, P::η, key);
}

template <class P, uint plain_modulus = P::plain_modulus>
void tlweSymIntEncrypt(TLWE<P> &res, const typename P::T p, const SecretKey &sk)
{
    tlweSymIntEncrypt<P, plain_modulus>(res, p, sk.key.get<P>());
}

template <class P>
typename P::T tlweSymPhase(const TLWE<P> &c, const Key<P> &key)
{
    typename P::T phase = c[P::k * P::n];
    for (int k = 0; k < P::k; k++)
        for (int i = 0; i < P::n; i++)
            phase -= c[k * P::n + i] * key[k * P::n + i];
    return phase;
}

template <class P>
bool tlweSymDecrypt(const TLWE<P> &c, const Key<P> &key)
{
    typename P::T phase = tlweSymPhase<P>(c, key);
    bool res =
        static_cast<typename std::make_signed<typename P::T>::type>(phase) > 0;
    return res;
}

template <class P>
bool tlweSymDecrypt(const TLWE<P> &c, const SecretKey &sk)
{
    return tlweSymDecrypt<P>(c, sk.key.get<P>());
}

template <class P, uint plain_modulus = P::plain_modulus>
typename P::T tlweSymIntDecrypt(const TLWE<P> &c, const Key<P> &key)
{
    constexpr double Δ =
        2 *
        static_cast<double>(
            static_cast<typename P::T>(1) << (std::numeric_limits<typename P::T>::digits - 1)) /
        plain_modulus;
    const typename P::T phase = tlweSymPhase<P>(c, key);
    typename P::T res = static_cast<typename P::T>(std::round(phase / Δ));
    return res >= plain_modulus / 2 ? res - plain_modulus : res;
}

template <class P, uint plain_modulus = P::plain_modulus>
typename P::T tlweSymIntDecrypt(const TLWE<P> &c, const SecretKey &sk)
{
    return tlweSymIntDecrypt<P, plain_modulus>(c, sk.key.get<P>());
}

template <class P, std::make_signed_t<typename P::T> μ>
void bootsSymEncrypt(std::vector<TLWE<P>> &c, const std::vector<uint8_t> &p,
                     const Key<P> &key)
{
    c.resize(p.size());
#pragma omp parallel for
    for (int i = 0; i < p.size(); i++)
        tlweSymEncrypt<P>(c[i], p[i] ? μ : -μ, key);
}

template <class P>
void bootsSymEncrypt(std::vector<TLWE<P>> &c, const std::vector<uint8_t> &p,
                     const Key<P> &key)
{
    bootsSymEncrypt<P, P::μ>(c, p, key);
}

template <class P = lvl1param>
void bootsSymEncrypt(std::vector<TLWE<P>> &c, const std::vector<uint8_t> &p,
                     const SecretKey &sk)
{
    bootsSymEncrypt<P>(c, p, sk.key.get<P>());
}

template <class P = lvl1param, std::make_signed_t<typename P::T> μ>
void bootsSymEncrypt(std::vector<TLWE<P>> &c, const std::vector<uint8_t> &p,
                     const SecretKey &sk)
{
    bootsSymEncrypt<P, μ>(c, p, sk.key.get<P>());
}

template <class P>
std::vector<uint8_t> bootsSymDecrypt(const std::vector<TLWE<P>> &c,
                                     const Key<P> &key)
{
    vector<uint8_t> p(c.size());
#pragma omp parallel for
    for (int i = 0; i < c.size(); i++) p[i] = tlweSymDecrypt<P>(c[i], key);
    return p;
}

template <class P = lvl1param>
std::vector<uint8_t> bootsSymDecrypt(const std::vector<TLWE<P>> &c,
                                     const SecretKey &sk)
{
    return bootsSymDecrypt<P>(c, sk.key.get<P>());
}

/**
 * @brief Adds an arbitrary number of TLWE ciphertexts element-wise.
 *
 * This function calculates res = c_1 + c_2 + ... + c_n.
 * It uses a C++17 fold expression within a loop to sum the elements
 * of all provided ciphertexts directly into the result.
 *
 * @tparam P The TLWE parameter type.
 * @tparam Args A parameter pack of TLWE<P> types.
 * @param res The output ciphertext where the sum is stored.
 * @param first The first ciphertext in the sum.
 * @param rest The remaining ciphertexts in the sum.
 */
template <class P, class... Args>
void TLWEAdd(TLWE<P> &res, const TLWE<P> &first, const Args &...rest)
{
    for (int i = 0; i <= P::k * P::n; i++) {
        // A binary fold expression sums all corresponding elements at once.
        res[i] = (first[i] + ... + rest[i]);
    }
}

/**
 * @brief Subtracts multiple TLWE ciphertexts element-wise.
 *
 * Calculates res = c1 - c2 - c3 - ...
 * NOTE: This implementation requires at least two input ciphertexts.
 */
template <class P, class... Args>
void TLWESub(TLWE<P> &res, const TLWE<P> &first, const Args &...rest)
{
    // A binary fold requires the parameter pack 'rest' to be non-empty.
    static_assert(
        sizeof...(Args) > 0,
        "This TLWESub implementation requires at least two arguments.");

    for (int i = 0; i <= P::k * P::n; i++) {
        // Binary fold over the '-' operator.
        // Expands to (((first[i] - c2[i]) - c3[i]) - ...)
        res[i] = (first[i] - ... - rest[i]);
    }
}

}  // namespace TFHEpp
