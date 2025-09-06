#pragma once
#include <AES.h>

#include <ranges>
// Transciphering by AES
//  Based on Hippogryph
namespace TFHEpp {

constexpr uint aesNb = 4;

inline void SubWord(std::array<uint8_t, 4> &a)
{
    int i;
    for (i = 0; i < 4; i++) {
        a[i] = sbox[a[i] / 16][a[i] % 16];
    }
}

inline void RotWord(std::array<uint8_t, 4> &a)
{
    uint8_t c = a[0];
    a[0] = a[1];
    a[1] = a[2];
    a[2] = a[3];
    a[3] = c;
}

inline uint8_t xtime(uint8_t b)  // multiply on x
{
    return (b << 1) ^ (((b >> 7) & 1) * 0x1b);
}

inline uint8_t Rcon(unsigned int n)
{
    uint8_t c = 1;
    for (uint i = 0; i < n - 1; i++) {
        c = xtime(c);
    }
    return c;
}

constexpr uint Nk = 4;  // Number of 32-bit words in the key
constexpr uint Nb = 4;  // Number of columns (32-bit words) comprising the state
constexpr uint Nr = 10;  // Number of rounds, which is a function of Nk and Nb

inline void KeyExpansion(std::array<uint8_t, 4 * Nb *(Nr + 1)> &w,
                  const std::array<uint8_t, 16> &key)
{
    std::array<uint8_t, 4> temp;

    unsigned int i = 0;
    while (i < 4 * Nk) {
        w[i] = key[i];
        i++;
    }

    i = 4 * Nk;
    while (i < 4 * Nb * (Nr + 1)) {
        temp[0] = w[i - 4 + 0];
        temp[1] = w[i - 4 + 1];
        temp[2] = w[i - 4 + 2];
        temp[3] = w[i - 4 + 3];

        if (i / 4 % Nk == 0) {
            RotWord(temp);
            SubWord(temp);
            temp[0] ^= Rcon(i / (Nk * 4));
        }
        else if (Nk > 6 && i / 4 % Nk == 4) {
            SubWord(temp);
        }

        w[i + 0] = w[i - 4 * Nk] ^ temp[0];
        w[i + 1] = w[i + 1 - 4 * Nk] ^ temp[1];
        w[i + 2] = w[i + 2 - 4 * Nk] ^ temp[2];
        w[i + 3] = w[i + 3 - 4 * Nk] ^ temp[3];
        i += 4;
    }
}

template <class P>
inline std::array<Polynomial<P>, (1 << 8) / (P::n / 8)> AESSboxROMPoly()
{
    std::array<Polynomial<P>, (1 << 8) / (P::n / 8)> polys;
    for (int i = 0; i < (1 << 8) / (P::n / 8); i++)
        for (int j = 0; j < P::n / 8; j++)
            for (int k = 0; k < 8; k++) {
                const uint index = i * (P::n / 8) + j;
                polys[i][j * 8 + k] =
                    ((sbox[index >> 4][index & 0xf] >> k) & 0x1)
                        ? 1ULL << (std::numeric_limits<typename P::T>::digits -
                                   2)
                        : -(1ULL
                            << (std::numeric_limits<typename P::T>::digits -
                                2));
            }
    return polys;
}

template <class iksP, class brP, class ahP>
void AESSboxROM(const std::span<TLWE<typename brP::targetP>, 8> res,
                const std::span<const TLWE<typename iksP::domainP>, 8> tlwe,
                const EvalKey &ek)
{
    constexpr uint32_t address_bit = 8;  // Address by words.
    constexpr uint32_t words_bit = 3;
    constexpr uint32_t width_bit =
        brP::targetP::nbit -
        words_bit;  // log_2 of how many words are in one TRLWE message.
    static_assert(address_bit >= width_bit);
    alignas(64) std::array<TRGSWFFT<typename brP::targetP>, 8> trgsw;
    for (int i = 0; i < 8; i++) {
        TLWE<typename iksP::domainP> shifted = tlwe[i];
        // shifted[iksP::targetP::k * iksP::targetP::n] -=
        // 1ULL << (std::numeric_limits<typename iksP::domainP::T>::digits - 2);
        if (i >= width_bit)
            for (int j = 0; j <= iksP::domainP::k * iksP::domainP::n; j++)
                shifted[j] *= -1;
        AnnihilateCircuitBootstrappingFFT<iksP, brP, ahP>(trgsw[i], shifted, ek);
    }
    std::array<TRLWE<typename brP::targetP>, (1 << 8) / (brP::targetP::n / 8)>
        rom = {};
    for (int i = 0; i < (1 << 8) / (brP::targetP::n / 8); i++)
        rom[i][brP::targetP::k] = AESSboxROMPoly<typename brP::targetP>()[i];
    if constexpr (address_bit != width_bit) {
        TRLWE<typename brP::targetP> trlwe;
        UROMUX<typename brP::targetP, address_bit, width_bit>(trlwe, trgsw,
                                                              rom);
        LROMUX<typename brP::targetP, address_bit, width_bit>(res, trgsw,
                                                              trlwe);
    }
    else {
        LROMUX<typename brP::targetP, address_bit, width_bit>(res, trgsw,
                                                              rom[0]);
    }
}

template <class iksP, class brP, class ahP>
void AESSboxROM(std::array<TLWE<typename brP::targetP>, 8> &res,
                const std::array<TLWE<typename iksP::domainP>, 8> &tlwe,
                const EvalKey &ek)
{
    AESSboxROM<iksP, brP, ahP>(std::span(res), std::span(tlwe), ek);
}

template <class iksP, class brP, class ahP>
void SubBytes(std::array<TLWE<typename brP::targetP>, 128> &state,
              const EvalKey &ek)
{
    for (int i = 0; i < 16; i++)
        AESSboxROM<iksP, brP, ahP>(
            std::span(state).subspan(i * 8).template first<8>(),
            std::span(state).subspan(i * 8).template first<8>(), ek);
}

template <class iksP, class brP, class cbiksP, class cbbrP, class ahP>
void SubBytes(std::array<TLWE<typename brP::targetP>, 128> &state,
              const EvalKey &ek)
{
    for (int i = 0; i < 16; i++) {
        std::array<TLWE<typename cbbrP::targetP>, 8> temp;
        AESSboxROM<iksP, cbbrP, ahP>(
            std::span(temp),
            std::span(state).subspan(i * 8).template first<8>(), ek);
        for (int j = 0; j < 8; j++)
            GateBootstrapping<
                cbiksP, brP,
                1ULL << (std::numeric_limits<typename brP::targetP::T>::digits -
                         2)>(state[i * 8 + j], temp[j], ek);
    }
}

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

template <class iksP, class brP, class ahP>
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
    {
        TRLWE<typename brP::targetP> trlwe;
        std::array<std::array<TLWE<typename iksP::domainP>, 16>, 2> tabletlwe;
        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 16; j++) tabletlwe[i][j] = midtlwes[j][i];
        TLWE2TablePackingManyLUT<ahP, 16, 2>(
            trlwe, tabletlwe, ek.getahk<ahP>());
        GateBootstrappingManyLUT<brP, 2>(res, shifted, ek.getbkfft<brP>(),
                                         trlwe);
    }
}

template <class P>
inline std::array<Polynomial<P>, (1 << 8) / (P::n / 8)> AESInvSboxROMPoly()
{
    std::array<Polynomial<P>, (1 << 8) / (P::n / 8)> polys;
    for (int i = 0; i < (1 << 8) / (P::n / 8); i++)
        for (int j = 0; j < P::n / 8; j++)
            for (int k = 0; k < 8; k++) {
                const uint index = i * (P::n / 8) + j;
                polys[i][j * 8 + k] =
                    ((inv_sbox[index >> 4][index & 0xf] >> k) & 0x1)
                        ? 1ULL << (std::numeric_limits<typename P::T>::digits -
                                   2)
                        : -(1ULL
                            << (std::numeric_limits<typename P::T>::digits -
                                2));
            }
    return polys;
}

template <class iksP, class brP, class ahP>
void AESInvSboxROM(std::array<TLWE<typename brP::targetP>, 8> &res,
                   const std::array<TLWE<typename iksP::domainP>, 8> &tlwe,
                   const EvalKey &ek)
{
    constexpr uint32_t address_bit = 8;  // Address by words.
    constexpr uint32_t words_bit = 3;
    constexpr uint32_t width_bit =
        brP::targetP::nbit -
        words_bit;  // log_2 of how many words are in one TRLWE message.
    static_assert(address_bit >= width_bit);
    alignas(64) std::array<TRGSWFFT<typename brP::targetP>, 8> trgsw;
    for (int i = 0; i < 8; i++) {
        TLWE<typename iksP::domainP> shifted = tlwe[i];
        // shifted[iksP::targetP::k * iksP::targetP::n] -=
        // 1ULL << (std::numeric_limits<typename iksP::domainP::T>::digits - 2);
        if (i >= width_bit)
            for (int j = 0; j <= iksP::domainP::k * iksP::domainP::n; j++)
                shifted[j] *= -1;
        AnnihilateCircuitBootstrappingFFT<iksP, brP, ahP>(trgsw[i], shifted, ek);
    }
    std::array<TRLWE<typename brP::targetP>, (1 << 8) / (brP::targetP::n / 8)>
        rom = {};
    for (int i = 0; i < (1 << 8) / (brP::targetP::n / 8); i++)
        rom[i][brP::targetP::k] = AESInvSboxROMPoly<typename brP::targetP>()[i];
    if constexpr (address_bit != width_bit) {
        TRLWE<typename brP::targetP> trlwe;
        UROMUX<typename brP::targetP, address_bit, width_bit>(trlwe, trgsw,
                                                              rom);
        LROMUX<typename brP::targetP, address_bit, width_bit>(res, trgsw,
                                                              trlwe);
    }
    else {
        LROMUX<typename brP::targetP, address_bit, width_bit>(res, trgsw,
                                                              rom[0]);
    }
}

template <class P, uint index, uint n>
void ShiftRow(std::array<TLWE<P>, 128> &res)
{
    std::array<TLWE<P>, Nb * 8> tmp;
    for (int i = 0; i < Nb; i++)
        for (int j = 0; j < 8; j++)
            tmp[i * 8 + j] = res[index * Nb * 8 + ((i + n) % Nb) * 8 + j];
    for (int i = 0; i < Nb * 8; i++) res[index * Nb * 8 + i] = tmp[i];
}

template <class P>
void ShiftRows(std::array<TLWE<P>, 128> &res)
{
    ShiftRow<P, 1, 1>(res);
    ShiftRow<P, 2, 2>(res);
    ShiftRow<P, 3, 3>(res);
}

template <class P>
void InvShiftRows(std::array<TLWE<P>, 128> &res)
{
    ShiftRow<P, 1, aesNb - 1>(res);
    ShiftRow<P, 2, aesNb - 2>(res);
    ShiftRow<P, 3, aesNb - 3>(res);
}

template <class P>
void MixColumn(std::array<TLWE<P>, 32> &y_out, const std::array<TLWE<P>, 32> &x)
{
    // Temporary variables for intermediate XOR results, based on the 92-gate
    // circuit Naming corresponds to t0...t59 and y0...y31 as per the circuit.
    // We will use 't' for all intermediate values and 'y_out' for the final 32
    // output bits.
    std::array<TLWE<P>, 60> t;  // For t0...t59

    // Implement the 92 XOR gates from Listing 1 in the provided PDF [cite: 28,
    // 29] Note: The paper uses various symbols for XOR (e.g., ^, ~, -, ´). We
    // interpret all as XOR. The indices for x will be 0-31.

    // Implement the 92 XOR gates based on the user-provided circuit
    TLWEAdd<P>(t[0], x[0], x[8]);
    TLWEAdd<P>(t[1], x[16], x[24]);
    TLWEAdd<P>(t[2], x[1], x[9]);
    TLWEAdd<P>(t[3], x[17], x[25]);
    TLWEAdd<P>(t[4], x[2], x[10]);
    TLWEAdd<P>(t[5], x[18], x[26]);
    TLWEAdd<P>(t[6], x[3], x[11]);
    TLWEAdd<P>(t[7], x[19], x[27]);
    TLWEAdd<P>(t[8], x[4], x[12]);
    TLWEAdd<P>(t[9], x[20], x[28]);
    TLWEAdd<P>(t[10], x[5], x[13]);
    TLWEAdd<P>(t[11], x[21], x[29]);
    TLWEAdd<P>(t[12], x[6], x[14]);
    TLWEAdd<P>(t[13], x[22], x[30]);
    TLWEAdd<P>(t[14], x[23], x[31]);
    TLWEAdd<P>(t[15], x[7], x[15]);
    TLWEAdd<P>(t[16], x[8], t[1]);
    TLWEAdd<P>(y_out[0], t[15], t[16]);
    TLWEAdd<P>(t[17], x[7], x[23]);
    TLWEAdd<P>(t[18], x[24], t[0]);
    TLWEAdd<P>(y_out[16], t[14], t[18]);
    TLWEAdd<P>(t[19], t[1], y_out[16]);
    TLWEAdd<P>(y_out[24], t[17], t[19]);
    TLWEAdd<P>(t[20], x[27], t[14]);
    TLWEAdd<P>(t[21], t[0], y_out[0]);
    TLWEAdd<P>(y_out[8], t[17], t[21]);
    TLWEAdd<P>(t[22], t[5], t[20]);
    TLWEAdd<P>(y_out[19], t[6], t[22]);
    TLWEAdd<P>(t[23], x[11], t[15]);
    TLWEAdd<P>(t[24], t[7], t[23]);
    TLWEAdd<P>(y_out[3], t[4], t[24]);
    TLWEAdd<P>(t[25], x[2], x[18]);
    TLWEAdd<P>(t[26], t[17], t[25]);
    TLWEAdd<P>(t[27], t[9], t[23]);
    TLWEAdd<P>(t[28], t[8], t[20]);
    TLWEAdd<P>(t[29], x[10], t[2]);
    TLWEAdd<P>(y_out[2], t[5], t[29]);
    TLWEAdd<P>(t[30], x[26], t[3]);
    TLWEAdd<P>(y_out[18], t[4], t[30]);
    TLWEAdd<P>(t[31], x[9], x[25]);
    TLWEAdd<P>(t[32], t[25], t[31]);
    TLWEAdd<P>(y_out[10], t[30], t[32]);
    TLWEAdd<P>(y_out[26], t[29], t[32]);
    TLWEAdd<P>(t[33], x[1], t[18]);
    TLWEAdd<P>(t[34], x[30], t[11]);
    TLWEAdd<P>(y_out[22], t[12], t[34]);
    TLWEAdd<P>(t[35], x[14], t[13]);
    TLWEAdd<P>(y_out[6], t[10], t[35]);
    TLWEAdd<P>(t[36], x[5], x[21]);
    TLWEAdd<P>(t[37], x[30], t[17]);
    TLWEAdd<P>(t[38], x[17], t[16]);
    TLWEAdd<P>(t[39], x[13], t[8]);
    TLWEAdd<P>(y_out[5], t[11], t[39]);
    TLWEAdd<P>(t[40], x[12], t[36]);
    TLWEAdd<P>(t[41], x[29], t[9]);
    TLWEAdd<P>(y_out[21], t[10], t[41]);
    TLWEAdd<P>(t[42], x[28], t[40]);
    TLWEAdd<P>(y_out[13], t[41], t[42]);
    TLWEAdd<P>(y_out[29], t[39], t[42]);
    TLWEAdd<P>(t[43], x[15], t[12]);
    TLWEAdd<P>(y_out[7], t[14], t[43]);
    TLWEAdd<P>(t[44], x[14], t[37]);
    TLWEAdd<P>(y_out[31], t[43], t[44]);
    TLWEAdd<P>(t[45], x[31], t[13]);
    TLWEAdd<P>(y_out[15], t[44], t[45]);
    TLWEAdd<P>(y_out[23], t[15], t[45]);
    TLWEAdd<P>(t[46], t[12], t[36]);
    TLWEAdd<P>(y_out[14], y_out[6], t[46]);
    TLWEAdd<P>(t[47], t[31], t[33]);
    TLWEAdd<P>(y_out[17], t[19], t[47]);
    TLWEAdd<P>(t[48], t[6], y_out[3]);
    TLWEAdd<P>(y_out[11], t[26], t[48]);
    TLWEAdd<P>(t[49], t[2], t[38]);
    TLWEAdd<P>(y_out[25], y_out[24], t[49]);
    TLWEAdd<P>(t[50], t[7], y_out[19]);
    TLWEAdd<P>(y_out[27], t[26], t[50]);
    TLWEAdd<P>(t[51], x[22], t[46]);
    TLWEAdd<P>(y_out[30], t[11], t[51]);
    TLWEAdd<P>(t[52], x[19], t[28]);
    TLWEAdd<P>(y_out[20], x[28], t[52]);
    TLWEAdd<P>(t[53], x[3], t[27]);
    TLWEAdd<P>(y_out[4], x[12], t[53]);
    TLWEAdd<P>(t[54], t[3], t[33]);
    TLWEAdd<P>(y_out[9], y_out[8], t[54]);
    TLWEAdd<P>(t[55], t[21], t[31]);
    TLWEAdd<P>(y_out[1], t[38], t[55]);
    TLWEAdd<P>(t[56], x[4], t[17]);
    TLWEAdd<P>(t[57], x[19], t[56]);
    TLWEAdd<P>(y_out[12], t[27], t[57]);
    TLWEAdd<P>(t[58], x[3], t[28]);
    TLWEAdd<P>(t[59], t[17], t[58]);
    TLWEAdd<P>(y_out[28], x[20], t[59]);
}

// https://eprint.iacr.org/2024/1076
template <class P>
void MixColumnDepth4(std::array<TLWE<P>, 32> &y,
                     const std::array<TLWE<P>, 32> &x)
{
    // r0 … r64
    std::array<TLWE<P>, 65> r;

    // --- first stage --------------------------------------------------------
    TLWEAdd<P>(r[0], x[23], x[31]);
    TLWEAdd<P>(r[1], x[21], x[29]);
    TLWEAdd<P>(r[2], x[17], x[25]);
    TLWEAdd<P>(r[3], x[16], x[24]);
    TLWEAdd<P>(r[4], x[15], x[23]);
    TLWEAdd<P>(r[5], x[14], x[22]);
    TLWEAdd<P>(r[6], x[12], x[20]);
    TLWEAdd<P>(r[7], x[12], x[13]);
    TLWEAdd<P>(r[8], x[11], x[20]);
    TLWEAdd<P>(r[9], x[10], x[25]);
    TLWEAdd<P>(r[10], x[10], x[18]);
    TLWEAdd<P>(r[11], x[9], x[18]);
    TLWEAdd<P>(r[12], x[7], x[31]);
    TLWEAdd<P>(r[13], x[7], x[15]);
    TLWEAdd<P>(r[14], x[6], x[31]);
    TLWEAdd<P>(r[15], x[6], x[30]);
    TLWEAdd<P>(r[16], x[5], x[13]);
    TLWEAdd<P>(r[17], x[5], x[6]);
    TLWEAdd<P>(r[18], x[4], x[28]);
    TLWEAdd<P>(r[19], x[3], x[27]);
    TLWEAdd<P>(r[20], x[3], x[11]);
    TLWEAdd<P>(r[21], x[2], x[26]);
    TLWEAdd<P>(r[22], x[1], x[9]);
    TLWEAdd<P>(r[23], x[0], x[8]);

    // --- second stage -------------------------------------------------------
    TLWEAdd<P>(r[24], r[0], x[27]);
    TLWEAdd<P>(r[25], r[0], x[7]);
    TLWEAdd<P>(r[26], r[1], x[5]);
    TLWEAdd<P>(r[27], r[1], x[4]);
    TLWEAdd<P>(r[28], r[2], x[1]);
    TLWEAdd<P>(r[29], r[3], x[8]);
    TLWEAdd<P>(r[30], r[3], r[0]);
    TLWEAdd<P>(r[31], r[4], x[14]);
    TLWEAdd<P>(r[32], r[4], x[7]);
    TLWEAdd<P>(r[33], r[4], x[0]);
    TLWEAdd<P>(r[34], r[5], x[29]);
    TLWEAdd<P>(r[35], r[6], x[28]);
    TLWEAdd<P>(r[36], r[6], x[4]);
    TLWEAdd<P>(r[37], r[10], x[26]);
    TLWEAdd<P>(r[38], r[10], r[4]);
    TLWEAdd<P>(r[39], r[13], r[2]);
    TLWEAdd<P>(r[40], r[15], x[22]);
    TLWEAdd<P>(r[41], r[16], x[30]);
    TLWEAdd<P>(r[42], r[16], x[28]);
    TLWEAdd<P>(r[43], r[18], x[19]);
    TLWEAdd<P>(r[44], r[19], x[19]);
    TLWEAdd<P>(r[45], r[19], r[12]);
    TLWEAdd<P>(r[46], r[20], x[26]);
    TLWEAdd<P>(r[47], r[20], r[13]);
    TLWEAdd<P>(r[48], r[21], x[17]);
    TLWEAdd<P>(r[49], r[22], x[25]);
    TLWEAdd<P>(r[50], r[22], x[17]);
    TLWEAdd<P>(r[51], r[23], x[16]);
    TLWEAdd<P>(r[52], r[23], x[9]);
    TLWEAdd<P>(r[53], r[24], x[18]);
    TLWEAdd<P>(r[54], r[24], x[12]);

    // --- outputs that depend only on r0 … r54 -------------------------------
    TLWEAdd<P>(y[15], r[25], r[5]);
    TLWEAdd<P>(y[13], r[26], r[6]);
    TLWEAdd<P>(y[5], r[27], r[7]);
    TLWEAdd<P>(y[0], r[29], r[13]);
    TLWEAdd<P>(y[7], r[31], r[14]);
    TLWEAdd<P>(y[31], r[32], r[15]);
    TLWEAdd<P>(y[8], r[33], r[3]);
    TLWEAdd<P>(y[30], r[34], r[17]);
    TLWEAdd<P>(y[2], r[37], r[22]);

    // --- remaining intermediates -------------------------------------------
    TLWEAdd<P>(r[55], r[37], r[28]);
    TLWEAdd<P>(r[56], r[40], x[21]);
    TLWEAdd<P>(r[57], r[40], r[13]);
    TLWEAdd<P>(y[6], r[41], r[5]);
    TLWEAdd<P>(r[58], r[42], x[29]);
    TLWEAdd<P>(r[59], r[43], r[4]);
    TLWEAdd<P>(r[60], r[44], x[2]);
    TLWEAdd<P>(y[11], r[44], r[38]);
    TLWEAdd<P>(y[28], r[45], r[36]);
    TLWEAdd<P>(r[61], r[46], r[45]);
    TLWEAdd<P>(r[62], r[47], x[10]);
    TLWEAdd<P>(y[4], r[47], r[35]);
    TLWEAdd<P>(y[18], r[48], r[9]);
    TLWEAdd<P>(y[10], r[48], r[11]);
    TLWEAdd<P>(y[17], r[49], r[30]);
    TLWEAdd<P>(r[63], r[50], r[29]);
    TLWEAdd<P>(y[24], r[51], r[12]);
    TLWEAdd<P>(y[16], r[51], r[30]);
    TLWEAdd<P>(r[64], r[51], r[33]);
    TLWEAdd<P>(y[1], r[52], r[39]);
    TLWEAdd<P>(y[19], r[53], r[46]);
    TLWEAdd<P>(y[20], r[54], r[43]);
    TLWEAdd<P>(y[26], r[55], r[48]);
    TLWEAdd<P>(y[14], r[56], x[13]);
    TLWEAdd<P>(y[22], r[56], r[34]);
    TLWEAdd<P>(y[23], r[57], r[14]);
    TLWEAdd<P>(y[21], r[58], x[20]);
    TLWEAdd<P>(y[29], r[58], r[27]);
    TLWEAdd<P>(y[12], r[59], r[8]);
    TLWEAdd<P>(y[27], r[61], r[60]);
    TLWEAdd<P>(y[3], r[62], r[60]);

    // y25 depends on y24 already produced
    TLWEAdd<P>(y[25], y[24], r[63]);

    TLWEAdd<P>(y[9], r[64], r[28]);
}

// https://eprint.iacr.org/2019/833
template <class P>
void MixColumns(std::array<TLWE<P>, 128> &state)
{
    // The MixColumns operation is applied to each 32-bit column of the state.
    // The AES state is 128 bits, so there are 4 such columns.

    for (int col = 0; col < 4; ++col) {
        // Extract the current 32-bit column into a working array (x0 to x31)
        std::array<TLWE<P>, 32> x;  // Input bits for the current column
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 8; ++j)
                x[i * 8 + j] = state[i * 32 + col * 8 + j];
        for (int i = 0; i < 32; ++i)
            x[i][P::k * P::n] +=
                1ULL << (std::numeric_limits<typename P::T>::digits - 2);

        // Apply the MixColumn transformation
        std::array<TLWE<P>, 32> y_out;  // For final output bits y0...y31
        MixColumn<P>(y_out, x);
        // MixColumnDepth4<P>(y_out, x);

        // Place the resulting 32-bit column (y_out) back into the state array
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 8; ++j)
                state[i * 32 + col * 8 + j] = y_out[i * 8 + j];
    }
    for (int i = 0; i < 128; i++)
        state[i][P::k * P::n] -=
            (1ULL << (std::numeric_limits<typename P::T>::digits - 2));
}

template <class P>
void xxtimes(std::array<TLWE<P>, 8> &statebyte)
{
    std::array<TLWE<P>, 8> tmp;
    tmp[0] = statebyte[6];
    TLWEAdd<P>(tmp[1], statebyte[6], statebyte[7]);
    TLWEAdd<P>(tmp[2], statebyte[0], statebyte[7]);
    TLWEAdd<P>(tmp[3], statebyte[1], statebyte[6]);
    TLWEAdd<P>(tmp[4], statebyte[2], tmp[1]);
    TLWEAdd<P>(tmp[5], statebyte[3], statebyte[7]);
    tmp[6] = statebyte[4];
    tmp[7] = statebyte[5];
    for (int i = 0; i < 8; i++) statebyte[i] = tmp[i];
}

// https://doi.org/10.1007/s13389-017-0176-3
template <class P>
void InvMixColumns(std::array<TLWE<P>, 128> &state)
{
    // The Inverse MixColumns operation is applied to each 32-bit column of the
    // state. The AES state is 128 bits, so there are 4 such columns.

    for (int col = 0; col < 4; ++col) {
        // Extract the current 32-bit column into a working array (x0 to x31)
        std::array<TLWE<P>, 32> x;  // Input bits for the current column
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 8; ++j)
                x[i * 8 + j] = state[i * 32 + col * 8 + j];
        for (int i = 0; i < 32; ++i)
            x[i][P::k * P::n] +=
                1ULL << (std::numeric_limits<typename P::T>::digits -
                         2);  // Shift to suppor XOR

        // Apply the Inverse MixColumn transformation (Decomposed to fix and
        // MixColumn)
        for (int offset = 0; offset < 2; ++offset) {
            std::array<TLWE<P>, 8> statebyte;
            for (int i = 0; i < 8; i++)
                TLWEAdd<P>(statebyte[i], x[offset * 8 + i],
                           x[(2 + offset) * 8 + i]);
            xxtimes<P>(statebyte);
            for (int i = 0; i < 8; i++)
                TLWEAdd<P>(x[offset * 8 + i], x[offset * 8 + i], statebyte[i]);
            for (int i = 0; i < 8; i++)
                TLWEAdd<P>(x[(2 + offset) * 8 + i], x[(2 + offset) * 8 + i],
                           statebyte[i]);
        }

        // Apply the MixColumn transformation
        std::array<TLWE<P>, 32> y_out;  // For final output bits y0...y31
        MixColumn<P>(y_out, x);

        // Place the resulting 32-bit column (y_out) back into the state array
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 8; ++j)
                state[i * 32 + col * 8 + j] = y_out[i * 8 + j];
    }
    for (int i = 0; i < 128; i++)
        state[i][P::k * P::n] -=
            (1ULL << (std::numeric_limits<typename P::T>::digits - 2));
}

template <class P>
void AddRoundKey(std::array<TLWE<P>, 128> &state,
                 const std::array<TLWE<P>, 128> &roundkey)
{
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < Nb; j++)
            for (int k = 0; k < 8; k++) {
                TLWEAdd<P>(state[i * Nb * 8 + j * 8 + k],
                           state[i * Nb * 8 + j * 8 + k],
                           roundkey[j * 4 * 8 + i * 8 + k]);
                state[i * Nb * 8 + j * 8 + k][P::k * P::n] +=
                    1ULL << (std::numeric_limits<typename P::T>::digits - 2);
            }
}

template <class iksP, class brP, class ahP>
void AESEnc(std::array<TLWE<typename brP::targetP>, 128> &cipher,
            const std::array<TLWE<typename iksP::domainP>, 128> &plain,
            const std::array<std::array<TLWE<typename brP::targetP>, 128>,
                             Nr + 1> &expandedkey,
            EvalKey &ek)
{
    std::array<TLWE<typename iksP::domainP>, 128> state;
    // Copy plaintext to state with transposition
    // Initial AddRoundKey
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < Nb; j++)
            for (int k = 0; k < 8; k++) {
                TLWEAdd<typename iksP::domainP>(
                    state[i * Nb * 8 + j * 8 + k], plain[j * 4 * 8 + i * 8 + k],
                    expandedkey[0][j * 4 * 8 + i * 8 + k]);
                state[i * Nb * 8 + j * 8 + k]
                     [iksP::domainP::k * iksP::domainP::n] +=
                    1ULL
                    << (std::numeric_limits<typename iksP::domainP::T>::digits -
                        2);
            }

    // Rounds
    for (int round = 1; round < Nr; round++) {
        SubBytes<iksP, brP, ahP>(state, ek);
        ShiftRows<typename brP::targetP>(state);
        MixColumns<typename brP::targetP>(state);
        AddRoundKey<typename brP::targetP>(state, expandedkey[round]);
    }
    SubBytes<iksP, brP, ahP>(state, ek);
    ShiftRows<typename brP::targetP>(state);
    AddRoundKey<typename brP::targetP>(state, expandedkey[Nr]);

    // Copy state to ciphertext with transposition
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < Nb; j++)
            for (int k = 0; k < 8; k++)
                cipher[j * 4 * 8 + i * 8 + k] = state[i * Nb * 8 + j * 8 + k];
}

template <class iksP, class brP, class cbiksP, class cbbrP, class ahP>
void AESEnc(std::array<TLWE<typename brP::targetP>, 128> &cipher,
            const std::array<TLWE<typename iksP::domainP>, 128> &plain,
            const std::array<std::array<TLWE<typename brP::targetP>, 128>,
                             Nr + 1> &expandedkey,
            EvalKey &ek)
{
    std::array<TLWE<typename iksP::domainP>, 128> state;
    // Copy plaintext to state with transposition
    // Initial AddRoundKey
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < Nb; j++)
            for (int k = 0; k < 8; k++) {
                TLWEAdd<typename iksP::domainP>(
                    state[i * Nb * 8 + j * 8 + k], plain[j * 4 * 8 + i * 8 + k],
                    expandedkey[0][j * 4 * 8 + i * 8 + k]);
                state[i * Nb * 8 + j * 8 + k]
                     [iksP::domainP::k * iksP::domainP::n] +=
                    1ULL
                    << (std::numeric_limits<typename iksP::domainP::T>::digits -
                        2);
            }

    // Rounds
    for (int round = 1; round < Nr; round++) {
        SubBytes<iksP, brP, cbiksP, cbbrP, ahP>(state, ek);
        ShiftRows<typename brP::targetP>(state);
        MixColumns<typename brP::targetP>(state);
        AddRoundKey<typename brP::targetP>(state, expandedkey[round]);
    }
    SubBytes<iksP, brP, cbiksP, cbbrP, ahP>(state, ek);
    ShiftRows<typename brP::targetP>(state);
    AddRoundKey<typename brP::targetP>(state, expandedkey[Nr]);

    // Copy state to ciphertext with transposition
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < Nb; j++)
            for (int k = 0; k < 8; k++)
                cipher[j * 4 * 8 + i * 8 + k] = state[i * Nb * 8 + j * 8 + k];
}

inline uint8_t Rcon(const uint8_t n)
{
    uint8_t rcon = 1;
    for (int i = 0; i < n - 1; i++) {
        rcon = xtime(rcon);
    }
    return rcon;
}

// Currently, assuming only AES128
// AES key expansion
template <class iksP, class brP, class ahP>
void KeyExpand(std::array<TLWE<typename brP::targetP>, 128> &next,
               const std::array<TLWE<typename iksP::domainP>, 128> &prev,
               EvalKey &ek, const uint8_t n)
{
    // SubWord & RotWord
    std::array<TLWE<typename iksP::domainP>, 8> sbox;
    for (int i = 0; i < 8; i++) sbox[i] = prev[96 + i];
    AESSboxROM<iksP, brP, ahP>(sbox, sbox, ek);
    for (int i = 0; i < 8; i++)
        TLWEAdd<typename iksP::domainP>(next[8 * 3 + i], sbox[i],
                                        prev[8 * 3 + i]);
    for (int i = 0; i < 8; i++) sbox[i] = prev[96 + 8 + i];
    AESSboxROM<iksP, brP, ahP>(sbox, sbox, ek);
    for (int i = 0; i < 8; i++)
        TLWEAdd<typename iksP::domainP>(next[i], sbox[i], prev[i]);
    // Add round constant
    const uint8_t rcon = Rcon(n);
    for (int i = 0; i < 8; i++)
        if ((rcon >> i) & 0x1)
            next[i][brP::targetP::k * brP::targetP::n] +=
                1ULL << (std::numeric_limits<typename brP::targetP::T>::digits -
                         1);
    for (int pos = 2; pos < 4; pos++) {
        for (int i = 0; i < 8; i++) sbox[i] = prev[96 + pos * 8 + i];
        AESSboxROM<iksP, brP, ahP>(sbox, sbox, ek);
        for (int i = 0; i < 8; i++)
            TLWEAdd<typename iksP::domainP>(next[(pos - 1) * 8 + i], sbox[i],
                                            prev[(pos - 1) * 8 + i]);
    }
    for (int i = 0; i < 32; i++)
        next[i][brP::targetP::k * brP::targetP::n] +=
            1ULL << (std::numeric_limits<typename brP::targetP::T>::digits - 2);

    for (int i = 1; i < 4; i++) {
        for (int j = 0; j < 32; j++) {
            TLWEAdd<typename iksP::domainP>(next[i * 32 + j], prev[i * 32 + j],
                                            next[(i - 1) * 32 + j]);
            next[i * 32 + j][brP::targetP::k * brP::targetP::n] +=
                1ULL << (std::numeric_limits<typename brP::targetP::T>::digits -
                         2);
        }
    }
}

template <class iksP, class brP, class ahP>
void KeyExpansion(
    std::array<std::array<TLWE<typename brP::targetP>, 128>, 10> &expandedkey,
    const std::array<TLWE<typename iksP::domainP>, 128> &key, EvalKey &ek)
{
    KeyExpand<iksP, brP,ahP>(expandedkey[0], key, ek, 1);
    for (int i = 1; i < 10; i++)
        KeyExpand<iksP, brP, ahP>(expandedkey[i], expandedkey[i - 1], ek, i + 1);
}

}  // namespace TFHEpp