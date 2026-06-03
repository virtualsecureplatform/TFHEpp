#pragma once

#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <span>
#include <type_traits>

namespace TFHEpp {

// NIST SP 800-232 / ascon-c main uses little-endian byte loading. Word
// operations are represented LSB-first inside each 64-bit word.
constexpr uint ascon_words = 5;
constexpr uint ascon_word_bits = 64;
constexpr uint ascon_state_bits = ascon_words * ascon_word_bits;
constexpr uint ascon_aead_rate_bytes = 16;
constexpr uint ascon_aead_rate_bits = ascon_aead_rate_bytes * 8;
constexpr uint ascon_xof_rate_bytes = 8;
constexpr uint ascon_xof_rate_bits = ascon_xof_rate_bytes * 8;
constexpr uint ascon_tag_bytes = 16;
constexpr uint ascon_key_bytes = 16;
constexpr uint ascon_nonce_bytes = 16;

constexpr std::array<uint8_t, 12> ascon_round_constants = {
    0xf0, 0xe1, 0xd2, 0xc3, 0xb4, 0xa5, 0x96, 0x87, 0x78, 0x69, 0x5a, 0x4b};

constexpr uint64_t ascon_128_iv =
    (static_cast<uint64_t>(1) << 0) | (static_cast<uint64_t>(12) << 16) |
    (static_cast<uint64_t>(8) << 20) |
    (static_cast<uint64_t>(ascon_tag_bytes * 8) << 24) |
    (static_cast<uint64_t>(ascon_aead_rate_bytes) << 40);

constexpr uint64_t ascon_xof_iv =
    (static_cast<uint64_t>(3) << 0) | (static_cast<uint64_t>(12) << 16) |
    (static_cast<uint64_t>(12) << 20) |
    (static_cast<uint64_t>(ascon_xof_rate_bytes) << 40);

template <class P>
using ASCONState = std::array<TLWE<P>, ascon_state_bits>;

template <class P>
inline constexpr typename P::T ascon_bit_mu =
    typename P::T{1} << (std::numeric_limits<typename P::T>::digits - 2);

constexpr std::size_t ASCONBitIndex(const std::size_t word,
                                    const std::size_t bit)
{
    return word * ascon_word_bits + bit;
}

constexpr std::size_t ASCONRateByteBitIndex(const std::size_t byte,
                                            const std::size_t bit)
{
    return ASCONBitIndex(byte / 8, 8 * (byte % 8) + bit);
}

constexpr std::size_t ASCONKeyByteBitIndex(const std::size_t byte,
                                           const std::size_t bit)
{
    return (byte < 8) ? ASCONBitIndex(0, 8 * byte + bit)
                      : ASCONBitIndex(1, 8 * (byte - 8) + bit);
}

template <class P>
void ASCONTrivialEncrypt(TLWE<P> &res, const bool bit)
{
    res = {};
    res[P::k * P::n] = bit ? ascon_bit_mu<P> : -ascon_bit_mu<P>;
}

template <class P>
void ASCONSetPublicState(ASCONState<P> &state,
                         const std::array<uint64_t, ascon_words> &words)
{
    for (std::size_t word = 0; word < ascon_words; word++)
        for (std::size_t bit = 0; bit < ascon_word_bits; bit++)
            ASCONTrivialEncrypt<P>(state[ASCONBitIndex(word, bit)],
                                   (words[word] >> bit) & 1);
}

template <class P>
void ASCONSetXOFInitialState(ASCONState<P> &state)
{
    ASCONSetPublicState<P>(state, {ascon_xof_iv, 0, 0, 0, 0});
}

template <class P>
void ASCONXOR(TLWE<P> &res, const TLWE<P> &a, const TLWE<P> &b)
{
    TLWEAdd<P>(res, a, b);
    res[P::k * P::n] += ascon_bit_mu<P>;
}

template <class P>
void ASCONXORInPlace(TLWE<P> &a, const TLWE<P> &b)
{
    TLWEAdd<P>(a, a, b);
    a[P::k * P::n] += ascon_bit_mu<P>;
}

template <class P>
void ASCONXORPublicBit(TLWE<P> &a, const bool bit)
{
    if (!bit) return;
    HomNOT<P>(a, a);
}

constexpr uint8_t ASCONSboxWord(const uint8_t input)
{
    std::array<uint8_t, ascon_words> x = {
        static_cast<uint8_t>((input >> 0) & 1),
        static_cast<uint8_t>((input >> 1) & 1),
        static_cast<uint8_t>((input >> 2) & 1),
        static_cast<uint8_t>((input >> 3) & 1),
        static_cast<uint8_t>((input >> 4) & 1)};

    std::array<uint8_t, ascon_words> t = {
        static_cast<uint8_t>(x[0] ^ ((x[1] ^ 1) & x[2])),
        static_cast<uint8_t>(x[1] ^ ((x[2] ^ 1) & x[3])),
        static_cast<uint8_t>(x[2] ^ ((x[3] ^ 1) & x[4])),
        static_cast<uint8_t>(x[3] ^ ((x[4] ^ 1) & x[0])),
        static_cast<uint8_t>(x[4] ^ ((x[0] ^ 1) & x[1]))};

    uint8_t output = 0;
    for (std::size_t word = 0; word < ascon_words; word++)
        output |= static_cast<uint8_t>(t[word] << word);
    return output;
}

template <class brP>
struct ASCONDefaultAHParam {
    using type = AHlvl1param;
};

template <>
struct ASCONDefaultAHParam<lvl02param> {
    using type = AHlvl2param;
};

template <>
struct ASCONDefaultAHParam<lvlh2param> {
    using type = AHlvl2param;
};

template <>
struct ASCONDefaultAHParam<cblvl02param> {
    using type = cbAHlvl2param;
};

template <class brP>
using ASCONDefaultAHParamT = typename ASCONDefaultAHParam<brP>::type;

template <class brP, class cbbrP,
          bool direct = std::is_same_v<typename brP::targetP,
                                       typename cbbrP::targetP>>
struct ASCONSboxAHParam {
    using type = ASCONDefaultAHParamT<cbbrP>;
};

template <class brP, class cbbrP>
struct ASCONSboxAHParam<brP, cbbrP, true> {
    using type = ASCONDefaultAHParamT<brP>;
};

template <class brP, class cbbrP>
using ASCONSboxAHParamT = typename ASCONSboxAHParam<brP, cbbrP>::type;

template <class P>
inline Polynomial<P> ASCONSboxROMPoly()
{
    Polynomial<P> poly = {};
    constexpr std::size_t address_bit = ascon_words;
    constexpr std::size_t segment = P::n / (std::size_t{1} << address_bit);
    static_assert(segment >= ascon_words);

    for (std::size_t input = 0; input < (std::size_t{1} << address_bit);
         input++) {
        const uint8_t output = ASCONSboxWord(static_cast<uint8_t>(input));
        for (std::size_t word = 0; word < ascon_words; word++)
            poly[input * segment + word] =
                ((output >> word) & 1) ? ascon_bit_mu<P> : -ascon_bit_mu<P>;
    }
    return poly;
}

template <class iksP, class brP, class ahP = ASCONDefaultAHParamT<brP>>
void ASCONSboxROM(std::array<TLWE<typename brP::targetP>, ascon_words> &out,
                  const std::array<TLWE<typename iksP::domainP>, ascon_words>
                      &in,
                  const EvalKey &ek)
{
    using P = typename brP::targetP;
    static_assert(std::is_same_v<typename iksP::targetP, typename brP::domainP>,
                  "ASCONSboxROM expects iksP::targetP == brP::domainP.");
    static_assert(P::k == ahP::k,
                  "ASCONSboxROM expects ahP::k == brP::targetP::k.");

    constexpr uint32_t address_bit = ascon_words;
    constexpr uint32_t width_bit = ascon_words;
    alignas(64) std::array<TRGSWFFT<P>, address_bit> trgsw;
    for (std::size_t word = 0; word < ascon_words; word++)
        AnnihilateCircuitBootstrapping<iksP, brP, ahP>(
            trgsw[word], in[word], ek);

    TRLWE<P> rom = {};
    rom[P::k] = ASCONSboxROMPoly<P>();
    LROMUX<P, address_bit, width_bit, ascon_words>(std::span(out), trgsw, rom);
}

template <class iksP, class brP, class cbiksP = lvl20param,
          class cbbrP = lvl02param,
          class ahP = ASCONSboxAHParamT<brP, cbbrP>>
void ASCONSbox(std::array<TLWE<typename brP::targetP>, ascon_words> &out,
               const std::array<TLWE<typename iksP::domainP>, ascon_words> &in,
               const EvalKey &ek)
{
    using P = typename brP::targetP;
    using cbP = typename cbbrP::targetP;
    static_assert(std::is_same_v<typename iksP::domainP, P>,
                  "ASCONSbox expects iksP::domainP == brP::targetP.");

    if constexpr (std::is_same_v<P, cbP>) {
        ASCONSboxROM<iksP, brP, ahP>(out, in, ek);
    }
    else {
        static_assert(std::is_same_v<typename cbiksP::domainP, cbP>,
                      "ASCONSbox expects cbiksP::domainP == cbbrP::targetP.");
        static_assert(
            std::is_same_v<typename cbiksP::targetP, typename brP::domainP>,
            "ASCONSbox expects cbiksP::targetP == brP::domainP.");
        std::array<TLWE<cbP>, ascon_words> rom_out;
        ASCONSboxROM<iksP, cbbrP, ahP>(rom_out, in, ek);
        for (std::size_t word = 0; word < ascon_words; word++)
            GateBootstrapping<cbiksP, brP, ascon_bit_mu<P>>(out[word],
                                                            rom_out[word], ek);
    }
}

template <class P>
void ASCONCopyRateBytes(std::span<TLWE<P>> out, const ASCONState<P> &state,
                        const std::size_t byte_count,
                        const std::size_t out_offset_bits = 0)
{
    assert(out.size() >= out_offset_bits + byte_count * 8);
    for (std::size_t byte = 0; byte < byte_count; byte++)
        for (std::size_t bit = 0; bit < 8; bit++)
            out[out_offset_bits + byte * 8 + bit] =
                state[ASCONRateByteBitIndex(byte, bit)];
}

template <class iksP, class brP, class cbiksP = lvl20param,
          class cbbrP = lvl02param,
          class ahP = ASCONSboxAHParamT<brP, cbbrP>>
void ASCONRound(ASCONState<typename brP::targetP> &state, const uint8_t C,
                const EvalKey &ek)
{
    using P = typename brP::targetP;
    static_assert(std::is_same_v<typename iksP::domainP, P>,
                  "ASCONRound expects iksP::domainP == brP::targetP.");

    for (std::size_t bit = 0; bit < 8; bit++)
        ASCONXORPublicBit<P>(state[ASCONBitIndex(2, bit)], (C >> bit) & 1);

    for (std::size_t bit = 0; bit < ascon_word_bits; bit++) {
        ASCONXORInPlace<P>(state[ASCONBitIndex(0, bit)],
                           state[ASCONBitIndex(4, bit)]);
        ASCONXORInPlace<P>(state[ASCONBitIndex(4, bit)],
                           state[ASCONBitIndex(3, bit)]);
        ASCONXORInPlace<P>(state[ASCONBitIndex(2, bit)],
                           state[ASCONBitIndex(1, bit)]);
    }

    auto t = std::make_unique<ASCONState<P>>();
    std::array<TLWE<P>, ascon_words> column;
    std::array<TLWE<P>, ascon_words> sboxed;
    for (std::size_t bit = 0; bit < ascon_word_bits; bit++) {
        for (std::size_t word = 0; word < ascon_words; word++)
            column[word] = state[ASCONBitIndex(word, bit)];
        ASCONSbox<iksP, brP, cbiksP, cbbrP, ahP>(sboxed, column, ek);
        for (std::size_t word = 0; word < ascon_words; word++)
            (*t)[ASCONBitIndex(word, bit)] = sboxed[word];
    }

    for (std::size_t bit = 0; bit < ascon_word_bits; bit++) {
        ASCONXORInPlace<P>((*t)[ASCONBitIndex(1, bit)],
                           (*t)[ASCONBitIndex(0, bit)]);
        ASCONXORInPlace<P>((*t)[ASCONBitIndex(0, bit)],
                           (*t)[ASCONBitIndex(4, bit)]);
        ASCONXORInPlace<P>((*t)[ASCONBitIndex(3, bit)],
                           (*t)[ASCONBitIndex(2, bit)]);
        ASCONXORPublicBit<P>((*t)[ASCONBitIndex(2, bit)], true);
    }

    constexpr std::array<uint, ascon_words> rot0 = {19, 61, 1, 10, 7};
    constexpr std::array<uint, ascon_words> rot1 = {28, 39, 6, 17, 41};
    TLWE<P> tmp;
    for (std::size_t word = 0; word < ascon_words; word++) {
        for (std::size_t bit = 0; bit < ascon_word_bits; bit++) {
            ASCONXOR<P>(tmp, (*t)[ASCONBitIndex(word, bit)],
                        (*t)[ASCONBitIndex(word, (bit + rot0[word]) & 63)]);
            ASCONXOR<P>(state[ASCONBitIndex(word, bit)], tmp,
                        (*t)[ASCONBitIndex(word, (bit + rot1[word]) & 63)]);
        }
    }
}

template <class iksP, class brP>
void ASCONPermute(ASCONState<typename brP::targetP> &state,
                  const std::size_t rounds, const EvalKey &ek)
{
    assert(rounds <= ascon_round_constants.size());
    const std::size_t begin = ascon_round_constants.size() - rounds;
    for (std::size_t i = begin; i < ascon_round_constants.size(); i++)
        ASCONRound<iksP, brP>(state, ascon_round_constants[i], ek);
}

template <class iksP, class brP>
void ASCONP12(ASCONState<typename brP::targetP> &state, const EvalKey &ek)
{
    ASCONPermute<iksP, brP>(state, 12, ek);
}

template <class iksP, class brP>
void ASCONP8(ASCONState<typename brP::targetP> &state, const EvalKey &ek)
{
    ASCONPermute<iksP, brP>(state, 8, ek);
}

template <class iksP, class brP>
void ASCONP6(ASCONState<typename brP::targetP> &state, const EvalKey &ek)
{
    ASCONPermute<iksP, brP>(state, 6, ek);
}

template <class P>
void ASCONXORKeyWord(ASCONState<P> &state, const std::size_t state_word,
                     std::span<const TLWE<P>> key, const std::size_t key_word)
{
    assert(key.size() >= ascon_key_bytes * 8);
    for (std::size_t byte = 0; byte < 8; byte++)
        for (std::size_t bit = 0; bit < 8; bit++)
            ASCONXORInPlace<P>(state[ASCONBitIndex(state_word, 8 * byte + bit)],
                               key[(key_word * 8 + byte) * 8 + bit]);
}

template <class iksP, class brP>
void ASCON128InitialState(ASCONState<typename brP::targetP> &state,
                          std::span<const TLWE<typename brP::targetP>> key,
                          std::span<const TLWE<typename brP::targetP>> nonce,
                          const EvalKey &ek)
{
    using P = typename brP::targetP;
    static_assert(
        std::is_same_v<typename iksP::domainP, P>,
        "ASCON128InitialState expects iksP::domainP == brP::targetP.");
    assert(key.size() >= ascon_key_bytes * 8);
    assert(nonce.size() >= ascon_nonce_bytes * 8);

    ASCONSetPublicState<P>(state, {ascon_128_iv, 0, 0, 0, 0});
    for (std::size_t byte = 0; byte < ascon_key_bytes; byte++)
        for (std::size_t bit = 0; bit < 8; bit++) {
            const std::size_t word = byte < 8 ? 1 : 2;
            const std::size_t word_byte = byte < 8 ? byte : byte - 8;
            state[ASCONBitIndex(word, 8 * word_byte + bit)] =
                key[byte * 8 + bit];
        }
    for (std::size_t byte = 0; byte < ascon_nonce_bytes; byte++)
        for (std::size_t bit = 0; bit < 8; bit++) {
            const std::size_t word = byte < 8 ? 3 : 4;
            const std::size_t word_byte = byte < 8 ? byte : byte - 8;
            state[ASCONBitIndex(word, 8 * word_byte + bit)] =
                nonce[byte * 8 + bit];
        }

    ASCONP12<iksP, brP>(state, ek);
    ASCONXORKeyWord<P>(state, 3, key, 0);
    ASCONXORKeyWord<P>(state, 4, key, 1);
}

template <class iksP, class brP>
void ASCONAbsorbAssociatedData(ASCONState<typename brP::targetP> &state,
                               std::span<const TLWE<typename brP::targetP>> ad,
                               const EvalKey &ek)
{
    using P = typename brP::targetP;
    static_assert(
        std::is_same_v<typename iksP::domainP, P>,
        "ASCONAbsorbAssociatedData expects iksP::domainP == brP::targetP.");
    assert(ad.size() % 8 == 0);

    std::size_t byte_len = ad.size() / 8;
    std::size_t offset_bits = 0;
    if (byte_len != 0) {
        while (byte_len >= ascon_aead_rate_bytes) {
            for (std::size_t byte = 0; byte < ascon_aead_rate_bytes; byte++)
                for (std::size_t bit = 0; bit < 8; bit++)
                    ASCONXORInPlace<P>(state[ASCONRateByteBitIndex(byte, bit)],
                                       ad[offset_bits + byte * 8 + bit]);
            ASCONP8<iksP, brP>(state, ek);
            offset_bits += ascon_aead_rate_bits;
            byte_len -= ascon_aead_rate_bytes;
        }
        for (std::size_t byte = 0; byte < byte_len; byte++)
            for (std::size_t bit = 0; bit < 8; bit++)
                ASCONXORInPlace<P>(state[ASCONRateByteBitIndex(byte, bit)],
                                   ad[offset_bits + byte * 8 + bit]);
        ASCONXORPublicBit<P>(state[ASCONRateByteBitIndex(byte_len, 0)], true);
        ASCONP8<iksP, brP>(state, ek);
    }
    ASCONXORPublicBit<P>(state[ASCONBitIndex(4, 63)], true);
}

template <class iksP, class brP>
void ASCONEncrypt(ASCONState<typename brP::targetP> &state,
                  std::span<TLWE<typename brP::targetP>> cipher,
                  std::span<const TLWE<typename brP::targetP>> plain,
                  const EvalKey &ek)
{
    using P = typename brP::targetP;
    static_assert(std::is_same_v<typename iksP::domainP, P>,
                  "ASCONEncrypt expects iksP::domainP == brP::targetP.");
    assert(plain.size() % 8 == 0);
    assert(cipher.size() >= plain.size());

    std::size_t byte_len = plain.size() / 8;
    std::size_t offset_bits = 0;
    while (byte_len >= ascon_aead_rate_bytes) {
        for (std::size_t byte = 0; byte < ascon_aead_rate_bytes; byte++)
            for (std::size_t bit = 0; bit < 8; bit++) {
                const std::size_t out = offset_bits + byte * 8 + bit;
                const std::size_t st = ASCONRateByteBitIndex(byte, bit);
                ASCONXORInPlace<P>(state[st], plain[out]);
                cipher[out] = state[st];
            }
        ASCONP8<iksP, brP>(state, ek);
        offset_bits += ascon_aead_rate_bits;
        byte_len -= ascon_aead_rate_bytes;
    }

    for (std::size_t byte = 0; byte < byte_len; byte++)
        for (std::size_t bit = 0; bit < 8; bit++) {
            const std::size_t out = offset_bits + byte * 8 + bit;
            const std::size_t st = ASCONRateByteBitIndex(byte, bit);
            ASCONXORInPlace<P>(state[st], plain[out]);
            cipher[out] = state[st];
        }
    ASCONXORPublicBit<P>(state[ASCONRateByteBitIndex(byte_len, 0)], true);
}

template <class iksP, class brP>
void ASCONDecrypt(ASCONState<typename brP::targetP> &state,
                  std::span<TLWE<typename brP::targetP>> plain,
                  std::span<const TLWE<typename brP::targetP>> cipher,
                  const EvalKey &ek)
{
    using P = typename brP::targetP;
    static_assert(std::is_same_v<typename iksP::domainP, P>,
                  "ASCONDecrypt expects iksP::domainP == brP::targetP.");
    assert(cipher.size() % 8 == 0);
    assert(plain.size() >= cipher.size());

    std::size_t byte_len = cipher.size() / 8;
    std::size_t offset_bits = 0;
    while (byte_len >= ascon_aead_rate_bytes) {
        for (std::size_t byte = 0; byte < ascon_aead_rate_bytes; byte++)
            for (std::size_t bit = 0; bit < 8; bit++) {
                const std::size_t out = offset_bits + byte * 8 + bit;
                const std::size_t st = ASCONRateByteBitIndex(byte, bit);
                ASCONXOR<P>(plain[out], state[st], cipher[out]);
                state[st] = cipher[out];
            }
        ASCONP8<iksP, brP>(state, ek);
        offset_bits += ascon_aead_rate_bits;
        byte_len -= ascon_aead_rate_bytes;
    }

    for (std::size_t byte = 0; byte < byte_len; byte++)
        for (std::size_t bit = 0; bit < 8; bit++) {
            const std::size_t out = offset_bits + byte * 8 + bit;
            const std::size_t st = ASCONRateByteBitIndex(byte, bit);
            ASCONXOR<P>(plain[out], state[st], cipher[out]);
            state[st] = cipher[out];
        }
    ASCONXORPublicBit<P>(state[ASCONRateByteBitIndex(byte_len, 0)], true);
}

template <class iksP, class brP>
void ASCONFinalize128(ASCONState<typename brP::targetP> &state,
                      std::span<TLWE<typename brP::targetP>> tag,
                      std::span<const TLWE<typename brP::targetP>> key,
                      const EvalKey &ek)
{
    using P = typename brP::targetP;
    static_assert(std::is_same_v<typename iksP::domainP, P>,
                  "ASCONFinalize128 expects iksP::domainP == brP::targetP.");
    assert(key.size() >= ascon_key_bytes * 8);
    assert(tag.size() >= ascon_tag_bytes * 8);

    ASCONXORKeyWord<P>(state, 2, key, 0);
    ASCONXORKeyWord<P>(state, 3, key, 1);
    ASCONP12<iksP, brP>(state, ek);
    ASCONXORKeyWord<P>(state, 3, key, 0);
    ASCONXORKeyWord<P>(state, 4, key, 1);

    for (std::size_t byte = 0; byte < 8; byte++)
        for (std::size_t bit = 0; bit < 8; bit++) {
            tag[byte * 8 + bit] = state[ASCONBitIndex(3, 8 * byte + bit)];
            tag[(8 + byte) * 8 + bit] = state[ASCONBitIndex(4, 8 * byte + bit)];
        }
}

template <class iksP, class brP>
void ASCONFinalize128WithKeyState(
    ASCONState<typename brP::targetP> &state,
    std::span<TLWE<typename brP::targetP>> tag,
    const ASCONState<typename brP::targetP> &key_state, const EvalKey &ek)
{
    // key_state words x1/x2 hold K0/K1. They are applied to x2/x3 before P12
    // and x3/x4 after P12, matching ascon-c's Ascon-AEAD128 finalization.
    using P = typename brP::targetP;
    static_assert(
        std::is_same_v<typename iksP::domainP, P>,
        "ASCONFinalize128WithKeyState expects iksP::domainP == brP::targetP.");
    assert(tag.size() >= ascon_tag_bytes * 8);

    for (std::size_t bit = 0; bit < ascon_word_bits; bit++) {
        ASCONXORInPlace<P>(state[ASCONBitIndex(2, bit)],
                           key_state[ASCONBitIndex(1, bit)]);
        ASCONXORInPlace<P>(state[ASCONBitIndex(3, bit)],
                           key_state[ASCONBitIndex(2, bit)]);
    }
    ASCONP12<iksP, brP>(state, ek);
    for (std::size_t bit = 0; bit < ascon_word_bits; bit++) {
        ASCONXORInPlace<P>(state[ASCONBitIndex(3, bit)],
                           key_state[ASCONBitIndex(1, bit)]);
        ASCONXORInPlace<P>(state[ASCONBitIndex(4, bit)],
                           key_state[ASCONBitIndex(2, bit)]);
    }

    for (std::size_t byte = 0; byte < 8; byte++)
        for (std::size_t bit = 0; bit < 8; bit++) {
            tag[byte * 8 + bit] = state[ASCONBitIndex(3, 8 * byte + bit)];
            tag[(8 + byte) * 8 + bit] = state[ASCONBitIndex(4, 8 * byte + bit)];
        }
}

template <class iksP, class brP>
void ASCONXOFInitialize(ASCONState<typename brP::targetP> &state,
                        const EvalKey &ek)
{
    ASCONSetXOFInitialState<typename brP::targetP>(state);
    ASCONP12<iksP, brP>(state, ek);
}

template <class iksP, class brP>
void ASCONXOFAbsorb(ASCONState<typename brP::targetP> &state,
                    std::span<const TLWE<typename brP::targetP>> input,
                    const EvalKey &ek)
{
    using P = typename brP::targetP;
    static_assert(std::is_same_v<typename iksP::domainP, P>,
                  "ASCONXOFAbsorb expects iksP::domainP == brP::targetP.");
    assert(input.size() % 8 == 0);

    std::size_t byte_len = input.size() / 8;
    std::size_t offset_bits = 0;
    while (byte_len >= ascon_xof_rate_bytes) {
        for (std::size_t byte = 0; byte < ascon_xof_rate_bytes; byte++)
            for (std::size_t bit = 0; bit < 8; bit++)
                ASCONXORInPlace<P>(state[ASCONRateByteBitIndex(byte, bit)],
                                   input[offset_bits + byte * 8 + bit]);
        ASCONP12<iksP, brP>(state, ek);
        offset_bits += ascon_xof_rate_bits;
        byte_len -= ascon_xof_rate_bytes;
    }
    for (std::size_t byte = 0; byte < byte_len; byte++)
        for (std::size_t bit = 0; bit < 8; bit++)
            ASCONXORInPlace<P>(state[ASCONRateByteBitIndex(byte, bit)],
                               input[offset_bits + byte * 8 + bit]);
    ASCONXORPublicBit<P>(state[ASCONRateByteBitIndex(byte_len, 0)], true);
    ASCONP12<iksP, brP>(state, ek);
}

template <class iksP, class brP>
void ASCONXOFSqueeze(ASCONState<typename brP::targetP> &state,
                     std::span<TLWE<typename brP::targetP>> output,
                     const EvalKey &ek)
{
    using P = typename brP::targetP;
    static_assert(std::is_same_v<typename iksP::domainP, P>,
                  "ASCONXOFSqueeze expects iksP::domainP == brP::targetP.");
    assert(output.size() % 8 == 0);

    std::size_t byte_len = output.size() / 8;
    std::size_t offset_bits = 0;
    while (byte_len > ascon_xof_rate_bytes) {
        ASCONCopyRateBytes<P>(output, state, ascon_xof_rate_bytes, offset_bits);
        ASCONP12<iksP, brP>(state, ek);
        offset_bits += ascon_xof_rate_bits;
        byte_len -= ascon_xof_rate_bytes;
    }
    ASCONCopyRateBytes<P>(output, state, byte_len, offset_bits);
}

template <class iksP, class brP>
void ASCONXOF(ASCONState<typename brP::targetP> &state,
              std::span<TLWE<typename brP::targetP>> output,
              std::span<const TLWE<typename brP::targetP>> input,
              const EvalKey &ek)
{
    ASCONXOFAbsorb<iksP, brP>(state, input, ek);
    ASCONXOFSqueeze<iksP, brP>(state, output, ek);
}

}  // namespace TFHEpp
