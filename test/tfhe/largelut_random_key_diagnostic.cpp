#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <tfhe++.hpp>

namespace {

using Clock = std::chrono::steady_clock;

template <class F>
double time_ms(F &&f)
{
    const auto start = Clock::now();
    f();
    const auto end = Clock::now();
    return std::chrono::duration<double, std::milli>(end - start).count();
}

template <class P>
long double plaintext_signed_error_units(const TFHEpp::TLWE<P> &ct,
                                         const uint32_t message,
                                         const TFHEpp::SecretKey &sk)
{
    constexpr uint32_t plain_modulus = 1 << 5;
    const auto phase = TFHEpp::tlweSymPhase<P>(ct, sk.key.get<P>());
    const long double modulus =
        std::ldexp(1.0L, std::numeric_limits<typename P::T>::digits);
    const long double delta = modulus / plain_modulus;
    long double diff =
        static_cast<long double>(phase) -
        static_cast<long double>(message % plain_modulus) * delta;
    diff = std::fmod(diff + modulus / 2, modulus);
    if (diff < 0) diff += modulus;
    diff -= modulus / 2;
    return diff / delta;
}

template <class T>
int64_t signed_plain(const T value)
{
    if (value > std::numeric_limits<T>::max() / 2)
        return -static_cast<int64_t>(T{} - value);
    return static_cast<int64_t>(value);
}

template <uint32_t K>
std::array<uint32_t, 1 << K> AESInvSboxNibbleTable(const bool high_nibble)
{
    std::array<uint32_t, 1 << K> table{};
    for (uint32_t x = 0; x < table.size(); x++) {
        const uint32_t byte = (x >> (K - 8)) & 0xFF;
        const uint8_t value = inv_sbox[byte >> 4][byte & 0xF];
        table[x] = high_nibble ? value >> 4 : value & 0xF;
    }
    return table;
}

template <uint32_t K>
std::array<uint32_t, 1 << K> IdentityByteNibbleTable(const bool high_nibble)
{
    std::array<uint32_t, 1 << K> table{};
    for (uint32_t x = 0; x < table.size(); x++) {
        const uint32_t byte = (x >> (K - 8)) & 0xFF;
        table[x] = high_nibble ? byte >> 4 : byte & 0xF;
    }
    return table;
}

template <class targetP, uint32_t W, uint32_t K, class Table>
TFHEpp::TRLWE<targetP> TableTRLWE(const Table &table)
{
    TFHEpp::TRLWE<targetP> trlwe;
    const auto poly = TFHEpp::LargeLUTPolynomial<targetP, W, K>(table);
    TFHEpp::TrivialTRLWEFromPolynomial<targetP>(trlwe, poly);
    return trlwe;
}

template <class P>
uint8_t decrypt_nibble_pair(const std::array<TFHEpp::TLWE<P>, 2> &ct,
                            const TFHEpp::SecretKey &sk)
{
    constexpr uint32_t plain_modulus = 1 << 5;
    const uint8_t low =
        TFHEpp::tlweSymIntDecrypt<P, plain_modulus>(ct[0], sk) & 0xF;
    const uint8_t high =
        TFHEpp::tlweSymIntDecrypt<P, plain_modulus>(ct[1], sk) & 0xF;
    return low | (high << 4);
}

uint32_t inv_sbox_preimage(const uint8_t value)
{
    for (uint32_t x = 0; x < 256; x++)
        if (inv_sbox[x >> 4][x & 0xF] == value) return x;
    return 256;
}

}  // namespace

int main()
{
    using brP = TFHEpp::lvlh2param;
    using iksP = TFHEpp::lvl2hparam;
    using ahP = TFHEpp::AHlvl2param;
    using domainP = typename brP::domainP;
    using targetP = typename brP::targetP;

    static_assert(std::is_same_v<typename iksP::domainP, targetP>);
    static_assert(std::is_same_v<typename iksP::targetP, domainP>);

    constexpr uint32_t W = 4;
    constexpr uint32_t K = 11;
    constexpr uint32_t L = TFHEpp::LargeLUTDigitCount<W, K>();
    constexpr auto offsets = TFHEpp::LargeLUTDigitOffsets<W, K>();
    constexpr uint32_t plain_modulus = 1 << (W + 1);
    constexpr uint32_t trials = 8;
    constexpr uint8_t input = 0x53;
    const uint8_t expected = inv_sbox[input >> 4][input & 0xF];
    constexpr uint32_t largelut_input =
        (static_cast<uint32_t>(input) << (K - 8)) | (uint32_t{1} << (K - 9));
    constexpr auto decomposed =
        TFHEpp::LargeLUTRadixDecompose<W, K>(largelut_input);

    const auto aes_low =
        TableTRLWE<targetP, W, K>(AESInvSboxNibbleTable<K>(false));
    const auto aes_high =
        TableTRLWE<targetP, W, K>(AESInvSboxNibbleTable<K>(true));
    const auto id_low =
        TableTRLWE<targetP, W, K>(IdentityByteNibbleTable<K>(false));
    const auto id_high =
        TableTRLWE<targetP, W, K>(IdentityByteNibbleTable<K>(true));

    uint32_t failures = 0;
    double total_keygen_ms = 0;
    double total_eval_ms = 0;

    std::cout << std::fixed << std::setprecision(3);
    std::cout << "LargeLUT random-key diagnostic, brP=lvlh2param, "
              << "iksP=lvl2hparam, ahP=AHlvl2param, W=" << W << ", K=" << K
              << ", trials=" << trials << "\n";
    std::cout << "input=0x" << std::hex << static_cast<uint32_t>(input)
              << ", largelut_input=0x" << largelut_input
              << ", expected=0x" << static_cast<uint32_t>(expected)
              << std::dec << "\n";
    std::cout << "decomposed digits:";
    for (uint32_t i = 0; i < L; i++)
        std::cout << " d" << i << "=" << decomposed[i] << "@"
                  << offsets[i];
    std::cout << "\n";

    for (uint32_t trial = 0; trial < trials; trial++) {
        auto sk = std::make_unique<TFHEpp::SecretKey>();
        TFHEpp::EvalKey ek;

        total_keygen_ms += time_ms([&] {
            ek.emplacebkfft<brP>(*sk);
            ek.emplaceiksk<iksP>(*sk);
            ek.emplaceahk<ahP>(*sk);
        });

        std::array<TFHEpp::TLWE<domainP>, L> encrypted_digits;
        for (uint32_t i = 0; i < L; i++)
            TFHEpp::tlweSymIntEncrypt<domainP, plain_modulus>(
                encrypted_digits[i], decomposed[i], 0.0,
                sk->key.get<domainP>());

        std::array<TFHEpp::TLWE<domainP>, L> calibrated;
        std::array<TFHEpp::TLWE<targetP>, 2> aes_res;
        std::array<TFHEpp::TLWE<targetP>, 2> id_res;

        total_eval_ms += time_ms([&] {
            TFHEpp::LargeLUTCalibrate<brP, iksP, W, K>(
                calibrated, encrypted_digits, ek.getbkfft<brP>(),
                ek.getiksk<iksP>());
            TFHEpp::LargeLUTWithCalibrated<brP, ahP, W, K>(
                aes_res[0], aes_low, calibrated, ek.getbkfft<brP>(),
                ek.getahk<ahP>());
            TFHEpp::LargeLUTWithCalibrated<brP, ahP, W, K>(
                aes_res[1], aes_high, calibrated, ek.getbkfft<brP>(),
                ek.getahk<ahP>());
            TFHEpp::LargeLUTWithCalibrated<brP, ahP, W, K>(
                id_res[0], id_low, calibrated, ek.getbkfft<brP>(),
                ek.getahk<ahP>());
            TFHEpp::LargeLUTWithCalibrated<brP, ahP, W, K>(
                id_res[1], id_high, calibrated, ek.getbkfft<brP>(),
                ek.getahk<ahP>());
        });

        const uint8_t aes_out = decrypt_nibble_pair<targetP>(aes_res, *sk);
        const uint8_t selected = decrypt_nibble_pair<targetP>(id_res, *sk);
        const bool ok = aes_out == expected && selected == input;
        failures += ok ? 0 : 1;

        long double max_digit_error = 0;
        long double estimated_index = 0;
        std::array<int64_t, L> cal_dec{};
        std::array<long double, L> cal_err{};
        for (uint32_t i = 0; i < L; i++) {
            const auto raw =
                TFHEpp::tlweSymIntDecrypt<domainP, plain_modulus>(
                    calibrated[i], *sk);
            cal_dec[i] = signed_plain(raw);
            cal_err[i] = plaintext_signed_error_units<domainP>(
                calibrated[i], decomposed[i], *sk);
            max_digit_error = std::max(max_digit_error, std::fabs(cal_err[i]));
            estimated_index +=
                (static_cast<long double>(decomposed[i]) + cal_err[i]) *
                static_cast<long double>(uint32_t{1} << offsets[i]);
        }

        std::cout << "trial " << trial << ": "
                  << (ok ? "ok" : "FAIL") << ", aes=0x" << std::hex
                  << static_cast<uint32_t>(aes_out) << ", selected=0x"
                  << static_cast<uint32_t>(selected) << ", inv_preimage=0x"
                  << inv_sbox_preimage(aes_out) << std::dec
                  << ", max calibrated error/Delta="
                  << static_cast<double>(max_digit_error)
                  << ", est_index=" << static_cast<double>(estimated_index)
                  << ", est_byte=0x" << std::hex
                  << static_cast<uint32_t>(
                         std::llround(estimated_index) >> (K - 8))
                  << std::dec << ", cal=[";
        for (uint32_t i = 0; i < L; i++) {
            if (i != 0) std::cout << ",";
            std::cout << cal_dec[i];
        }
        std::cout << "], cal_signed_err=[";
        for (uint32_t i = 0; i < L; i++) {
            if (i != 0) std::cout << ",";
            std::cout << static_cast<double>(cal_err[i]);
        }
        std::cout << "]\n";
    }

    std::cout << "failures: " << failures << "/" << trials << "\n";
    std::cout << "avg keygen ms: " << total_keygen_ms / trials << "\n";
    std::cout << "avg eval ms: " << total_eval_ms / trials << "\n";
}
