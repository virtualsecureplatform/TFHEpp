#include "keyswitch.hpp"

#include <cstdint>
#include <limits>

#include "utils.hpp"

namespace TFHEpp {
using namespace std;

template <class P>
void IdentityKeySwitch(TLWE<typename P::targetP> &res,
                       const TLWE<typename P::domainP> &tlwe,
                       const KeySwitchingKey<P> &ksk)
{
    constexpr typename P::domainP::T prec_offset =
        1ULL << (numeric_limits<typename P::domainP::T>::digits -
                 (1 + P::basebit * P::t));
    constexpr uint32_t mask = (1U << P::basebit) - 1;
    res = {};
    constexpr uint32_t domain_digit =
        std::numeric_limits<typename P::domainP::T>::digits;
    constexpr uint32_t target_digit =
        std::numeric_limits<typename P::targetP::T>::digits;
    if constexpr (domain_digit == target_digit)
        res[P::targetP::n] = tlwe[P::domainP::n];
    else if constexpr (domain_digit > target_digit)
        res[P::targetP::n] = (tlwe[P::domainP::n] +
                              (1ULL << (domain_digit - target_digit - 1))) >>
                             (domain_digit - target_digit);
    else if constexpr (domain_digit < target_digit)
        res[P::targetP::n] = tlwe[P::domainP::n]
                             << (target_digit - domain_digit);
    for (int i = 0; i < P::domainP::n; i++) {
        const typename P::domainP::T aibar = tlwe[i] + prec_offset;
        for (int j = 0; j < P::t; j++) {
            const uint32_t aij =
                (aibar >> (numeric_limits<typename P::domainP::T>::digits -
                           (j + 1) * P::basebit)) &
                mask;
            if (aij != 0)
                for (int k = 0; k <= P::targetP::n; k++)
                    res[k] -= ksk[i][j][aij - 1][k];
        }
    }
}
#define INST(P)                                                               \
    template void IdentityKeySwitch<P>(TLWE<typename P::targetP> & res,       \
                                       const TLWE<typename P::domainP> &tlwe, \
                                       const KeySwitchingKey<P> &ksk)
TFHEPP_EXPLICIT_INSTANTIATION_KEY_SWITCH(INST)
#undef INST

template <class P>
void PrivKeySwitch(TRLWE<typename P::targetP> &res,
                   const TLWE<typename P::domainP> &tlwe,
                   const PrivKeySwitchKey<P> &privksk)
{
    constexpr uint32_t mask = (1 << P::basebit) - 1;
    constexpr uint64_t prec_offset =
        1ULL << (std::numeric_limits<typename P::domainP::T>::digits -
                 (1 + P::basebit * P::t));

    res = {};
    for (int i = 0; i <= P::domainP::n; i++) {
        const typename P::domainP::T aibar = tlwe[i] + prec_offset;

        for (int j = 0; j < P::t; j++) {
            const typename P::domainP::T aij =
                (aibar >> (std::numeric_limits<typename P::domainP::T>::digits -
                           (j + 1) * P::basebit)) &
                mask;

            if (aij != 0) {
                for (int p = 0; p < P::targetP::n; p++) {
                    res[0][p] -= privksk[i][j][aij - 1][0][p];
                    res[1][p] -= privksk[i][j][aij - 1][1][p];
                }
            }
        }
    }
}
#define INST(P)                                                           \
    template void PrivKeySwitch<P>(TRLWE<typename P::targetP> & res,      \
                                   const TLWE<typename P::domainP> &tlwe, \
                                   const PrivKeySwitchKey<P> &privksk)
TFHEPP_EXPLICIT_INSTANTIATION_KEY_SWITCH(INST)
#undef INST

}  // namespace TFHEpp
