/**
 * @author Kotaro Matsuoka
 */

#pragma once

#include "gatebootstrapping.hpp"

namespace TFHEpp {
/*!
 * @brief Generates a Polynomial with each coefficient subtracted by a base
 * value
 * @tparam P The parameter set for the Polynomial
 * @tparam basebit The base to be subtracted
 * @return A Polynomial of the parameter type with coefficients subtracted by
 * the base value
 */
template <class P, uint basebit>
constexpr Polynomial<P> subtractpolygen()
{
    Polynomial<P> poly;
    for (int i = 0; i < P::n; i++)
        poly[i] =
            1ULL << (std::numeric_limits<typename P::T>::digits - basebit - 2);
    return poly;
}

// https://eprint.iacr.org/2023/645
/*!
 * @brief Homomorphically decomposes an input ciphertext into an array of level
 * 1 ciphertexts
 * @tparam high2midP The parameter set for the transition from high to mid level
 * @tparam mid2lowP The parameter set for the transition from mid to low level
 * @tparam brP The bootstrapping parameter set
 * @tparam basebit The base value
 * @tparam numdigit The number of digits
 * @param cres Array of output ciphertexts
 * @param cin Input ciphertext
 * @param kskh2m The key switching key for high to mid level
 * @param kskm2l The key switching key for mid to low level
 * @param bkfft The bootstrapping key FFT
 */
template <class high2midP, class mid2lowP, class brP, uint basebit,
          uint numdigit>
void HomDecomp(std::array<TLWE<typename high2midP::targetP>, numdigit> &cres,
               const TLWE<typename high2midP::domainP> &cin,
               const KeySwitchingKey<high2midP> &kskh2m,
               const KeySwitchingKey<mid2lowP> &kskm2l,
               const BootstrappingKeyFFT<brP> &bkfft)
{
    TFHEpp::TLWE<typename mid2lowP::targetP> tlwelvlhalf;
    TFHEpp::TLWE<typename high2midP::targetP> subtlwe;

    // cres will be used as a reusable buffer
    constexpr uint32_t plain_modulusbit = basebit * numdigit;
#pragma omp parallel for default(none) shared(cin, cres, kskh2m)
    for (int digit = 1; digit <= numdigit; digit++) {
        TFHEpp::TLWE<typename high2midP::domainP> switchedtlwe;
        for (int i = 0; i <= high2midP::domainP::k * high2midP::domainP::n; i++)
            switchedtlwe[i] = cin[i] << (plain_modulusbit - basebit * digit);
        IdentityKeySwitch<high2midP>(cres[digit - 1], switchedtlwe, kskh2m);
    }
    for (int digit = 1; digit <= numdigit; digit++) {
        if (digit != 1) {
            for (int i = 0; i <= high2midP::targetP::k * high2midP::targetP::n;
                 i++)
                cres[digit - 1][i] += subtlwe[i];
            cres[digit - 1][high2midP::targetP::k * high2midP::targetP::n] -=
                1ULL << (std::numeric_limits<
                             typename high2midP::targetP::T>::digits -
                         basebit - 1);
        }
        IdentityKeySwitch<mid2lowP>(tlwelvlhalf, cres[digit - 1], kskm2l);
        tlwelvlhalf[mid2lowP::targetP::k * mid2lowP::targetP::n] +=
            1ULL
            << (std::numeric_limits<typename mid2lowP::targetP::T>::digits -
                basebit - 1);
        if (digit != numdigit)
            GateBootstrappingTLWE2TLWEFFT<brP>(
                subtlwe, tlwelvlhalf, bkfft,
                subtractpolygen<typename high2midP::targetP, basebit>());
    }
}

}  // namespace TFHEpp