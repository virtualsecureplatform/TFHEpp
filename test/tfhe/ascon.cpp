#include <algorithm>
#include <array>
#include <cassert>
#include <cctype>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <span>
#include <string>
#include <tfhe++.hpp>
#include <vector>

extern "C" {
int crypto_aead_encrypt(unsigned char *c, unsigned long long *clen,
                        const unsigned char *m, unsigned long long mlen,
                        const unsigned char *ad, unsigned long long adlen,
                        const unsigned char *nsec, const unsigned char *npub,
                        const unsigned char *k);

int crypto_aead_decrypt(unsigned char *m, unsigned long long *mlen,
                        unsigned char *nsec, const unsigned char *c,
                        unsigned long long clen, const unsigned char *ad,
                        unsigned long long adlen, const unsigned char *npub,
                        const unsigned char *k);

int crypto_hash(unsigned char *out, const unsigned char *in,
                unsigned long long inlen);
}

namespace {

const unsigned char *bytes_or_null(const std::vector<uint8_t> &bytes)
{
    return bytes.empty() ? nullptr : bytes.data();
}

std::vector<uint8_t> from_hex(const std::string &hex)
{
    assert(hex.size() % 2 == 0);
    std::vector<uint8_t> out(hex.size() / 2);
    for (std::size_t i = 0; i < out.size(); i++) {
        const auto hi = static_cast<uint8_t>(
            std::isdigit(hex[2 * i]) ? hex[2 * i] - '0'
                                     : std::toupper(hex[2 * i]) - 'A' + 10);
        const auto lo =
            static_cast<uint8_t>(std::isdigit(hex[2 * i + 1])
                                     ? hex[2 * i + 1] - '0'
                                     : std::toupper(hex[2 * i + 1]) - 'A' + 10);
        out[i] = static_cast<uint8_t>((hi << 4) | lo);
    }
    return out;
}

std::string to_hex(const std::vector<uint8_t> &bytes)
{
    static constexpr char table[] = "0123456789ABCDEF";
    std::string out;
    out.reserve(bytes.size() * 2);
    for (const uint8_t byte : bytes) {
        out.push_back(table[byte >> 4]);
        out.push_back(table[byte & 0xf]);
    }
    return out;
}

std::vector<uint8_t> ascon_c_aead_encrypt(const std::vector<uint8_t> &key,
                                          const std::vector<uint8_t> &nonce,
                                          const std::vector<uint8_t> &ad,
                                          const std::vector<uint8_t> &plain)
{
    std::vector<uint8_t> cipher(plain.size() + TFHEpp::ascon_tag_bytes);
    unsigned long long clen = 0;
    const int rc = crypto_aead_encrypt(
        cipher.data(), &clen, bytes_or_null(plain), plain.size(),
        bytes_or_null(ad), ad.size(), nullptr, nonce.data(), key.data());
    assert(rc == 0);
    cipher.resize(clen);
    return cipher;
}

std::vector<uint8_t> ascon_c_aead_decrypt(const std::vector<uint8_t> &key,
                                          const std::vector<uint8_t> &nonce,
                                          const std::vector<uint8_t> &ad,
                                          const std::vector<uint8_t> &cipher)
{
    assert(cipher.size() >= TFHEpp::ascon_tag_bytes);
    std::vector<uint8_t> plain(cipher.size() - TFHEpp::ascon_tag_bytes);
    unsigned long long mlen = 0;
    const int rc = crypto_aead_decrypt(
        plain.data(), &mlen, nullptr, cipher.data(), cipher.size(),
        bytes_or_null(ad), ad.size(), nonce.data(), key.data());
    assert(rc == 0);
    plain.resize(mlen);
    return plain;
}

std::vector<uint8_t> ascon_c_xof(const std::vector<uint8_t> &message,
                                 const std::size_t outlen)
{
    std::vector<uint8_t> out(64);
    const int rc =
        crypto_hash(out.data(), bytes_or_null(message), message.size());
    assert(rc == 0);
    out.resize(outlen);
    return out;
}

template <class P>
void encrypt_bytes_as_bits(std::vector<TFHEpp::TLWE<P>> &out,
                           const std::vector<uint8_t> &bytes,
                           const TFHEpp::SecretKey &sk)
{
    out.resize(bytes.size() * 8);
    for (std::size_t byte = 0; byte < bytes.size(); byte++)
        for (std::size_t bit = 0; bit < 8; bit++)
            TFHEpp::tlweSymEncrypt<P>(out[byte * 8 + bit],
                                      ((bytes[byte] >> bit) & 1)
                                          ? TFHEpp::ascon_bit_mu<P>
                                          : -TFHEpp::ascon_bit_mu<P>,
                                      0.0, sk.key.get<P>());
}

template <class P>
std::vector<uint8_t> decrypt_bits_as_bytes(
    std::span<const TFHEpp::TLWE<P>> bits, const TFHEpp::SecretKey &sk)
{
    assert(bits.size() % 8 == 0);
    std::vector<uint8_t> bytes(bits.size() / 8);
    for (std::size_t byte = 0; byte < bytes.size(); byte++)
        for (std::size_t bit = 0; bit < 8; bit++)
            if (TFHEpp::tlweSymDecrypt<P>(bits[byte * 8 + bit], sk))
                bytes[byte] |= static_cast<uint8_t>(1U << bit);
    return bytes;
}

void compare_homomorphic_ascon_xof_with_ascon_c()
{
    using brP = TFHEpp::lvlh2param;
    using iksP = TFHEpp::lvl2hparam;
    using P = typename brP::targetP;

    const std::vector<uint8_t> empty;
    const auto ref_xof = ascon_c_xof(empty, TFHEpp::ascon_xof_rate_bytes);

    TFHEpp::SecretKey sk;
    TFHEpp::EvalKey ek;
    ek.emplacebkfft<brP>(sk);
    ek.emplaceiksk<iksP>(sk);
    ek.emplaceahk<TFHEpp::ASCONDefaultAHParamT<brP>>(sk);
    ek.emplacecbsk<TFHEpp::ASCONDefaultAHParamT<brP>>(sk);

    auto state = std::make_unique<TFHEpp::ASCONState<P>>();
    std::vector<TFHEpp::TLWE<P>> empty_bits;
    std::vector<TFHEpp::TLWE<P>> cxof(TFHEpp::ascon_xof_rate_bits);

    TFHEpp::ASCONXOFInitialize<iksP, brP>(*state, ek);
    TFHEpp::ASCONXOF<iksP, brP>(*state, cxof, empty_bits, ek);
    const auto actual_xof = decrypt_bits_as_bytes<P>(cxof, sk);
    if (actual_xof != ref_xof)
        std::cerr << "ASCON XOF mismatch expected=" << to_hex(ref_xof)
                  << " actual=" << to_hex(actual_xof) << std::endl;
    assert(actual_xof == ref_xof);
}

void compare_homomorphic_ascon_with_ascon_c()
{
    using brP = TFHEpp::lvlh2param;
    using iksP = TFHEpp::lvl2hparam;
    using P = typename brP::targetP;

    if (std::getenv("TFHEPP_RUN_ASCON_HOM_TEST") == nullptr) return;

    const auto key = from_hex("000102030405060708090A0B0C0D0E0F");
    const auto nonce = from_hex("101112131415161718191A1B1C1D1E1F");
    const std::vector<uint8_t> empty;
    const auto plain = from_hex("20");
    const auto ref_one = ascon_c_aead_encrypt(key, nonce, empty, plain);

    TFHEpp::SecretKey sk;
    TFHEpp::EvalKey ek;
    ek.emplacebkfft<brP>(sk);
    ek.emplaceiksk<iksP>(sk);
    ek.emplaceahk<TFHEpp::ASCONDefaultAHParamT<brP>>(sk);
    ek.emplacecbsk<TFHEpp::ASCONDefaultAHParamT<brP>>(sk);

    bool sbox_ok = true;
    for (std::size_t input = 0; input < (std::size_t{1} << TFHEpp::ascon_words);
         input++) {
        std::array<TFHEpp::TLWE<P>, TFHEpp::ascon_words> cin;
        std::array<TFHEpp::TLWE<P>, TFHEpp::ascon_words> cout;
        for (std::size_t word = 0; word < TFHEpp::ascon_words; word++)
            TFHEpp::tlweSymEncrypt<P>(
                cin[word], ((input >> word) & 1) ? TFHEpp::ascon_bit_mu<P>
                                                 : -TFHEpp::ascon_bit_mu<P>,
                0.0, sk.key.get<P>());
        TFHEpp::ASCONSbox<iksP, brP>(cout, cin, ek);
        uint8_t actual = 0;
        for (std::size_t word = 0; word < TFHEpp::ascon_words; word++)
            if (TFHEpp::tlweSymDecrypt<P>(cout[word], sk))
                actual |= static_cast<uint8_t>(1U << word);
        const uint8_t expected =
            TFHEpp::ASCONSboxWord(static_cast<uint8_t>(input));
        if (actual != expected)
            std::cerr << "ASCON S-box mismatch input=" << input
                      << " expected=" << static_cast<int>(expected)
                      << " actual=" << static_cast<int>(actual) << std::endl;
        sbox_ok &= actual == expected;
    }
    assert(sbox_ok);

    std::vector<TFHEpp::TLWE<P>> ckey;
    std::vector<TFHEpp::TLWE<P>> cnonce;
    std::vector<TFHEpp::TLWE<P>> cplain;
    encrypt_bytes_as_bits<P>(ckey, key, sk);
    encrypt_bytes_as_bits<P>(cnonce, nonce, sk);
    encrypt_bytes_as_bits<P>(cplain, plain, sk);

    auto state = std::make_unique<TFHEpp::ASCONState<P>>();
    std::vector<TFHEpp::TLWE<P>> ccipher(cplain.size());
    std::vector<TFHEpp::TLWE<P>> ctag(TFHEpp::ascon_tag_bytes * 8);
    std::vector<TFHEpp::TLWE<P>> empty_bits;

    TFHEpp::ASCON128InitialState<iksP, brP>(*state, ckey, cnonce, ek);
    TFHEpp::ASCONAbsorbAssociatedData<iksP, brP>(*state, empty_bits, ek);
    TFHEpp::ASCONEncrypt<iksP, brP>(*state, ccipher, cplain, ek);
    TFHEpp::ASCONFinalize128<iksP, brP>(*state, ctag, ckey, ek);
    auto actual_one = decrypt_bits_as_bytes<P>(ccipher, sk);
    const auto actual_tag = decrypt_bits_as_bytes<P>(ctag, sk);
    actual_one.insert(actual_one.end(), actual_tag.begin(), actual_tag.end());
    if (actual_one != ref_one)
        std::cerr << "ASCON AEAD encrypt mismatch expected=" << to_hex(ref_one)
                  << " actual=" << to_hex(actual_one) << std::endl;
    assert(actual_one == ref_one);

    const std::vector<uint8_t> ref_cipher(ref_one.begin(),
                                          ref_one.end() -
                                              TFHEpp::ascon_tag_bytes);
    std::vector<TFHEpp::TLWE<P>> ccipher_ref;
    std::vector<TFHEpp::TLWE<P>> cdecrypted(cplain.size());
    encrypt_bytes_as_bits<P>(ccipher_ref, ref_cipher, sk);

    TFHEpp::ASCON128InitialState<iksP, brP>(*state, ckey, cnonce, ek);
    TFHEpp::ASCONAbsorbAssociatedData<iksP, brP>(*state, empty_bits, ek);
    TFHEpp::ASCONDecrypt<iksP, brP>(*state, cdecrypted, ccipher_ref, ek);
    TFHEpp::ASCONFinalize128<iksP, brP>(*state, ctag, ckey, ek);
    const auto actual_plain = decrypt_bits_as_bytes<P>(cdecrypted, sk);
    const auto actual_decrypt_tag = decrypt_bits_as_bytes<P>(ctag, sk);
    const std::vector<uint8_t> ref_tag(ref_one.end() - TFHEpp::ascon_tag_bytes,
                                       ref_one.end());
    if (actual_plain != plain || actual_decrypt_tag != ref_tag)
        std::cerr << "ASCON AEAD decrypt mismatch expected_plain="
                  << to_hex(plain) << " actual_plain=" << to_hex(actual_plain)
                  << " expected_tag=" << to_hex(ref_tag)
                  << " actual_tag=" << to_hex(actual_decrypt_tag)
                  << std::endl;
    assert(actual_plain == plain);
    assert(actual_decrypt_tag == ref_tag);
}

}  // namespace

int main()
{
    const auto key = from_hex("000102030405060708090A0B0C0D0E0F");
    const auto nonce = from_hex("101112131415161718191A1B1C1D1E1F");
    const std::vector<uint8_t> empty;

    const auto ct_empty = ascon_c_aead_encrypt(key, nonce, empty, empty);
    assert(ct_empty == from_hex("4F9C278211BEC9316BF68F46EE8B2EC6"));
    assert(ascon_c_aead_decrypt(key, nonce, empty, ct_empty).empty());

    const auto ct_one = ascon_c_aead_encrypt(key, nonce, empty, from_hex("20"));
    assert(ct_one == from_hex("E8DD576ABA1CD3E6FC704DE02AEDB79588"));
    assert(ascon_c_aead_decrypt(key, nonce, empty, ct_one) == from_hex("20"));

    assert(ascon_c_xof(empty, 64) ==
           from_hex("473D5E6164F58B39DFD84AACDB8AE42E"
                    "C2D91FED33388EE0D960D9B3993295C6"
                    "AD77855A5D3B13FE6AD9E6098988373A"
                    "F7D0956D05A8F1665D2C67D1A3AD10FF"));

    compare_homomorphic_ascon_xof_with_ascon_c();
    compare_homomorphic_ascon_with_ascon_c();

    std::cout << "Passed" << std::endl;
}
