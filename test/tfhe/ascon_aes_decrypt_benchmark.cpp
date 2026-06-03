#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cctype>
#include <cstdlib>
#include <iomanip>
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
}

namespace {

using Clock = std::chrono::steady_clock;

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

std::array<uint8_t, 16> block_from_hex(const std::string &hex)
{
    const auto bytes = from_hex(hex);
    assert(bytes.size() == 16);
    std::array<uint8_t, 16> block;
    std::copy(bytes.begin(), bytes.end(), block.begin());
    return block;
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
void encrypt_block_as_bits(std::array<TFHEpp::TLWE<P>, 128> &out,
                           const std::array<uint8_t, 16> &bytes,
                           const TFHEpp::SecretKey &sk)
{
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

template <class P>
std::array<uint8_t, 16> decrypt_block_bits(
    const std::array<TFHEpp::TLWE<P>, 128> &bits, const TFHEpp::SecretKey &sk)
{
    std::array<uint8_t, 16> bytes = {};
    for (std::size_t byte = 0; byte < bytes.size(); byte++)
        for (std::size_t bit = 0; bit < 8; bit++)
            if (TFHEpp::tlweSymDecrypt<P>(bits[byte * 8 + bit], sk))
                bytes[byte] |= static_cast<uint8_t>(1U << bit);
    return bytes;
}

template <class P>
void encrypt_expanded_aes_key(
    std::array<std::array<TFHEpp::TLWE<P>, 128>, TFHEpp::Nr + 1> &out,
    const std::array<uint8_t, 16> &key, const TFHEpp::SecretKey &sk)
{
    std::array<uint8_t, 4 * TFHEpp::Nb *(TFHEpp::Nr + 1)> expanded;
    TFHEpp::KeyExpansion(expanded, key);
    for (std::size_t round = 0; round < TFHEpp::Nr + 1; round++)
        for (std::size_t byte = 0; byte < 16; byte++)
            for (std::size_t bit = 0; bit < 8; bit++)
                TFHEpp::tlweSymEncrypt<P>(
                    out[round][byte * 8 + bit],
                    ((expanded[round * 16 + byte] >> bit) & 1)
                        ? TFHEpp::ascon_bit_mu<P>
                        : -TFHEpp::ascon_bit_mu<P>,
                    0.0, sk.key.get<P>());
}

template <class F>
double time_ms(F &&f)
{
    const auto start = Clock::now();
    f();
    const auto end = Clock::now();
    return std::chrono::duration<double, std::milli>(end - start).count();
}

}  // namespace

int main(int argc, char **argv)
{
    using brP = TFHEpp::lvlh2param;
    using iksP = TFHEpp::lvl2hparam;
    using ahP = TFHEpp::AHlvl2param;
    using P = typename brP::targetP;

    int repetitions = 1;
    if (argc > 1) repetitions = std::max(1, std::atoi(argv[1]));

    const auto aes_key = block_from_hex("000102030405060708090A0B0C0D0E0F");
    const auto plain_block =
        block_from_hex("00112233445566778899AABBCCDDEEFF");
    const auto ascon_key = from_hex("000102030405060708090A0B0C0D0E0F");
    const auto ascon_nonce = from_hex("101112131415161718191A1B1C1D1E1F");
    const std::vector<uint8_t> plain(plain_block.begin(), plain_block.end());
    const std::vector<uint8_t> empty_ad;

    AES aes(AESKeyLength::AES_128);
    std::vector<unsigned char> aes_key_vec(aes_key.begin(), aes_key.end());
    std::vector<unsigned char> plain_vec(plain_block.begin(),
                                         plain_block.end());
    const auto aes_cipher_vec = aes.EncryptECB(plain_vec, aes_key_vec);
    assert(aes.DecryptECB(aes_cipher_vec, aes_key_vec) == plain_vec);

    std::array<uint8_t, 16> aes_cipher;
    std::copy(aes_cipher_vec.begin(), aes_cipher_vec.end(), aes_cipher.begin());

    const auto ascon_ref =
        ascon_c_aead_encrypt(ascon_key, ascon_nonce, empty_ad, plain);
    assert(ascon_ref.size() == plain.size() + TFHEpp::ascon_tag_bytes);
    const std::vector<uint8_t> ascon_cipher(ascon_ref.begin(),
                                            ascon_ref.begin() + plain.size());

    TFHEpp::SecretKey sk;
    TFHEpp::EvalKey ek;
    ek.emplacebkfft<brP>(sk);
    ek.emplaceiksk<iksP>(sk);
    ek.emplaceahk<ahP>(sk);
    ek.emplacecbsk<ahP>(sk);

    using AESBlock = std::array<TFHEpp::TLWE<P>, 128>;
    using AESRoundKeys = std::array<AESBlock, TFHEpp::Nr + 1>;
    auto caes_expanded_key = std::make_unique<AESRoundKeys>();
    auto caes_cipher = std::make_unique<AESBlock>();
    encrypt_expanded_aes_key<P>(*caes_expanded_key, aes_key, sk);
    encrypt_block_as_bits<P>(*caes_cipher, aes_cipher, sk);

    std::vector<TFHEpp::TLWE<P>> cascon_key;
    std::vector<TFHEpp::TLWE<P>> cascon_nonce;
    std::vector<TFHEpp::TLWE<P>> cascon_cipher;
    std::vector<TFHEpp::TLWE<P>> empty_bits;
    encrypt_bytes_as_bits<P>(cascon_key, ascon_key, sk);
    encrypt_bytes_as_bits<P>(cascon_nonce, ascon_nonce, sk);
    encrypt_bytes_as_bits<P>(cascon_cipher, ascon_cipher, sk);

    double aes_total_ms = 0;
    double ascon_total_ms = 0;
    auto caes_plain = std::make_unique<AESBlock>();
    std::vector<TFHEpp::TLWE<P>> cascon_plain(plain.size() * 8);
    auto ascon_state = std::make_unique<TFHEpp::ASCONState<P>>();

    for (int i = 0; i < repetitions; i++) {
        aes_total_ms += time_ms([&] {
            TFHEpp::AESDec<iksP, brP, ahP>(*caes_plain, *caes_cipher,
                                           *caes_expanded_key, ek);
        });

        ascon_total_ms += time_ms([&] {
            TFHEpp::ASCON128InitialState<iksP, brP>(*ascon_state, cascon_key,
                                                    cascon_nonce, ek);
            TFHEpp::ASCONAbsorbAssociatedData<iksP, brP>(*ascon_state,
                                                         empty_bits, ek);
            TFHEpp::ASCONDecrypt<iksP, brP>(*ascon_state, cascon_plain,
                                            cascon_cipher, ek);
        });
    }

    assert(decrypt_block_bits<P>(*caes_plain, sk) == plain_block);
    assert(decrypt_bits_as_bytes<P>(cascon_plain, sk) == plain);

    const double aes_avg_ms = aes_total_ms / repetitions;
    const double ascon_avg_ms = ascon_total_ms / repetitions;

    std::cout << std::fixed << std::setprecision(3);
    std::cout << "Homomorphic 128-bit decrypt benchmark\n";
    std::cout << "parameters: iksP=lvl2hparam brP=lvlh2param\n";
    std::cout << "repetitions: " << repetitions << "\n";
    std::cout << "AES-128 block decrypt ms: " << aes_avg_ms << "\n";
    std::cout << "ASCON-128 decrypt no-tag ms: " << ascon_avg_ms << "\n";
    std::cout << "ASCON/AES ratio: " << (ascon_avg_ms / aes_avg_ms) << "\n";
    std::cout << "Passed" << std::endl;
}
