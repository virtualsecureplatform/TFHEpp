#pragma once
#include <AES.h>
// Transciphering by AES
//  Based on Hippogryph
namespace TFHEpp {

template <class P>
inline Polynomial<P> AESInvSboxPoly(const uint8_t upperindex)
{
    Polynomial<P> poly;
    constexpr uint segment = P::n / 16;
    for (int i = 0; i < 16; i++)
        for (int j = 0; j < segment / 2; j++) {
            poly[i * segment + 2 * j] =
                (inv_sbox[upperindex][i] & 0xF) *
                (1ULL << (std::numeric_limits<typename P::T>::digits - 5));
            poly[i * segment + 2 * j + 1] =
                (inv_sbox[upperindex][i] >> 4) *
                (1ULL << (std::numeric_limits<typename P::T>::digits - 5));
        }
    return poly;
}

template <class iksP, class brP>
void AESInvSbox(std::array<TLWE<typename brP::targetP>, 2> &res,
                const std::array<TLWE<typename iksP::domainP>, 2> &tlwe,
                const EvalKey &ek)
{
    std::array<std::array<TLWE<typename brP::targetP>, 2>, 16> midtlwes;
    TLWE<typename iksP::targetP> shifted;
    IdentityKeySwitch<iksP>(shifted, tlwe[0], ek.getiksk<iksP>());
    shifted[iksP::targetP::k * iksP::targetP::n] +=
        1ULL << (std::numeric_limits<typename iksP::targetP::T>::digits - 6);
    for (int i = 0; i < 16; i++)
        GateBootstrappingManyLUT<brP, 2>(
            midtlwes[i], shifted, ek.getbkfft<brP>(),
            AESInvSboxPoly<typename brP::targetP>(i));
    IdentityKeySwitch<iksP>(shifted, tlwe[1], ek.getiksk<iksP>());
    shifted[iksP::targetP::k * iksP::targetP::n] +=
        1ULL << (std::numeric_limits<typename iksP::targetP::T>::digits - 6);
    for (int i = 0; i < 2; i++) {
        TRLWE<typename brP::targetP> trlwe;
        std::array<TLWE<typename iksP::domainP>, 16> tabletlwe;
        for (int j = 0; j < 16; j++) tabletlwe[j] = midtlwes[j][i];
        TLWE2TablePacking<typename brP::targetP, 16>(
            trlwe, tabletlwe, ek.getahk<typename brP::targetP>());
        GateBootstrappingTLWE2TLWEFFT<brP>(res[i], shifted, ek.getbkfft<brP>(),
                                           trlwe);
    }
}
}  // namespace TFHEpp