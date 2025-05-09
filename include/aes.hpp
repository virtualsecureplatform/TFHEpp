#pragma once
#include <AES.h>
// Transciphering by AES
//  Based on Hippogryph
namespace TFHEpp {

constexpr uint aesNb = 4;

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

template <class P>
inline std::array<Polynomial<P>, (1<<8)/(P::n/8)> AESInvSboxROMPoly()
{
    std::array<Polynomial<P>, (1<<8)/(P::n/8)> polys;
    for(int i = 0; i < (1<<8)/(P::n/8); i++) 
        for(int j = 0; j < P::n/8; j++) 
            for(int k = 0; k < 8; k++) {
                const uint index = i * (P::n/8) + j;
                polys[i][j*8+k] = ((inv_sbox[index>>4][index&0xf] >> k) & 0x1)?1ULL << (std::numeric_limits<typename P::T>::digits - 2):-(1ULL << (std::numeric_limits<typename P::T>::digits - 2));
            }
    return polys;
}

template <class iksP, class brP>
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
    for(int i = 0; i < 8; i++){
        TLWE<typename iksP::domainP> shifted = tlwe[i];
        // shifted[iksP::targetP::k * iksP::targetP::n] -=
            // 1ULL << (std::numeric_limits<typename iksP::domainP::T>::digits - 2);
        if(i >= width_bit) for(int j = 0; j <= iksP::domainP::k*iksP::domainP::n; j++)
            shifted[j] *= -1;
        AnnihilateCircuitBootstrappingFFT<iksP, brP>(trgsw[i], shifted, ek);
    }
    std::array<TRLWE<typename brP::targetP>, (1<<8)/(brP::targetP::n/8)> rom = {};
    for(int i = 0; i < (1<<8)/(brP::targetP::n/8); i++) rom[i][brP::targetP::k] = AESInvSboxROMPoly<typename brP::targetP>()[i];
    if constexpr (address_bit!=width_bit){
        TRLWE<typename brP::targetP> trlwe;
        UROMUX<typename brP::targetP, address_bit, width_bit>(
            trlwe, trgsw, rom);
        LROMUX<typename brP::targetP, address_bit, width_bit>(res, trgsw, trlwe);
    }else{
        LROMUX<typename brP::targetP, address_bit, width_bit>(res, trgsw, rom[0]);
    }
}

template <class P, uint index, uint n>
void ShiftRow(std::array<TLWE<P>, 128> &res)
{
    std::array<TLWE<P>, 32> tmp;
    for (int i = 0; i < aesNb; i++)
        for(int j = 0; j < 8; j++)
            tmp[i*8 + j] = res[index * aesNb * 8 + ((i+n)%aesNb)*8 + j];
    for (int i = 0; i < aesNb*8; i++)
        res[index * aesNb * 8 + i] = tmp[i];
}

template <class P>
void ShiftRows(std::array<TLWE<P>, 128> &res)
{
    ShiftRow<P,1,1>(res);
    ShiftRow<P,2,2>(res);
    ShiftRow<P,3,3>(res);
}

template <class P>
void InvShiftRows(std::array<TLWE<P>, 128> &res)
{
    ShiftRow<P,1,aesNb-1>(res);
    ShiftRow<P,2,aesNb-2>(res);
    ShiftRow<P,3,aesNb-3>(res);
}

// https://eprint.iacr.org/2019/833
template <class P>
void MixColumns(std::array<TLWE<P>, 128> &state) {
    // The MixColumns operation is applied to each 32-bit column of the state.
    // The AES state is 128 bits, so there are 4 such columns.

    for (int col = 0; col < 4; ++col) {
        // Extract the current 32-bit column into a working array (x0 to x31)
        std::array<TLWE<P>, 32> x; // Input bits for the current column
        for (int i = 0; i < 4; ++i) 
            for (int j = 0; j < 8; ++j)
                x[i*8+j] = state[i*32 + col*8 + j];
        for (int i = 0; i < 32; ++i) x[i][P::k * P::n] += 1ULL << (std::numeric_limits<typename P::T>::digits - 2);

        // Temporary variables for intermediate XOR results, based on the 92-gate circuit
        // Naming corresponds to t0...t59 and y0...y31 as per the circuit.
        // We will use 't' for all intermediate values and 'y_out' for the final 32 output bits.
        std::array<TLWE<P>, 60> t; // For t0...t59
        std::array<TLWE<P>, 32> y_out; // For final output bits y0...y31

        // Implement the 92 XOR gates from Listing 1 in the provided PDF [cite: 28, 29]
        // Note: The paper uses various symbols for XOR (e.g., ^, ~, -, ´). We interpret all as XOR.
        // The indices for x will be 0-31.

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

        // Place the resulting 32-bit column (y_out) back into the state array
        for (int i = 0; i < 4; ++i) 
            for (int j = 0; j < 8; ++j)
                state[i*32 + col*8 + j] = y_out[i*8+j];
    }
    for(int i = 0; i < 128; i++)
        state[i][P::k*P::n] -= (1ULL << (std::numeric_limits<typename P::T>::digits - 2));
}
}  // namespace TFHEpp