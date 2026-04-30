#include <array>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <tfhe++.hpp>
#include <vector>

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

double average(const std::vector<double> &values)
{
    return std::accumulate(values.begin(), values.end(), 0.0) / values.size();
}

template <class P, std::size_t N>
uint64_t checksum(const std::array<TFHEpp::TLWE<P>, N> &tlwes)
{
    uint64_t acc = 0;
    for (const auto &tlwe : tlwes)
        for (const auto coeff : tlwe) acc += static_cast<uint64_t>(coeff);
    return acc;
}

template <class P>
long double plaintext_error_units(const TFHEpp::TLWE<P> &ct,
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
    return std::fabs(diff) / delta;
}

template <class P>
long double max_nibble_pair_error_units(
    const std::array<TFHEpp::TLWE<P>, 2> &ct, const uint8_t expected,
    const TFHEpp::SecretKey &sk)
{
    return std::max(plaintext_error_units<P>(ct[0], expected & 0xF, sk),
                    plaintext_error_units<P>(ct[1], expected >> 4, sk));
}

template <class P, uint32_t t, uint32_t basebit>
constexpr std::size_t R2RKeyBytes()
{
    return P::n * t * (std::size_t{1} << (basebit - 1)) *
           sizeof(TFHEpp::TRLWE<P>);
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

template <class targetP, uint32_t W, uint32_t K>
TFHEpp::TRLWE<targetP> LargeLUTNibbleTRLWE(const bool high_nibble)
{
    TFHEpp::TRLWE<targetP> trlwe;
    const auto table = AESInvSboxNibbleTable<K>(high_nibble);
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

template <class P>
uint8_t decrypt_bits(const std::array<TFHEpp::TLWE<P>, 8> &ct,
                     const TFHEpp::SecretKey &sk)
{
    uint8_t res = 0;
    for (uint32_t i = 0; i < ct.size(); i++)
        if (TFHEpp::tlweSymDecrypt<P>(ct[i], sk)) res |= uint8_t{1} << i;
    return res;
}

}  // namespace

int main()
{
    using brP = TFHEpp::lvlh2param;
    using iksP = TFHEpp::lvl2hparam;
    using ahP = TFHEpp::AHlvl2param;
    using domainP = typename brP::domainP;
    using targetP = typename brP::targetP;
    constexpr uint32_t r2r_t = iksP::t;
    constexpr uint32_t r2r_basebit = iksP::basebit;

    static_assert(std::is_same_v<typename iksP::domainP, targetP>);
    static_assert(std::is_same_v<typename iksP::targetP, domainP>);

    constexpr uint32_t W = 4;
    constexpr uint32_t K = 11;
    constexpr uint32_t L = TFHEpp::LargeLUTDigitCount<W, K>();
    constexpr uint32_t plain_modulus = 1 << (W + 1);
    constexpr uint32_t repeats = 3;
    constexpr uint8_t input = 0x53;
    constexpr uint32_t largelut_input =
        (static_cast<uint32_t>(input) << (K - 8)) | (uint32_t{1} << (K - 9));

    auto sk = std::make_unique<TFHEpp::SecretKey>();
    TFHEpp::EvalKey ek;

    const double bk_keygen_ms =
        time_ms([&] { ek.emplacebkfft<brP>(*sk); });
    const double iksk_keygen_ms =
        time_ms([&] { ek.emplaceiksk<iksP>(*sk); });
    const double ah_keygen_ms =
        time_ms([&] { ek.emplaceahk<ahP>(*sk); });
    const double cbsk_keygen_ms =
        time_ms([&] { ek.emplacecbsk<ahP>(*sk); });

    using R2RStepKey = TFHEpp::R2RKey<targetP, r2r_t, r2r_basebit>;
    auto step_r2rks = std::make_unique<std::array<R2RStepKey, L - 1>>();
    std::array<double, L - 1> r2r_step_keygen_ms{};
    for (uint32_t i = 1; i < L; i++) {
        r2r_step_keygen_ms[i - 1] = time_ms([&] {
            const auto positions =
                TFHEpp::LargeLUTStepR2RPositions<targetP, W, K>(i);
            const uint32_t Bi =
                TFHEpp::LargeLUTStepBlockSize<targetP, W, K>(i);
            TFHEpp::R2RKeyGen<targetP, r2r_t, r2r_basebit>(
                (*step_r2rks)[i - 1], positions, Bi, sk->key.get<targetP>());
        });
    }
    const double r2r_keygen_ms =
        std::accumulate(r2r_step_keygen_ms.begin(),
                        r2r_step_keygen_ms.end(), 0.0);

    std::array<TFHEpp::TLWE<domainP>, L> largelut_digits;
    const auto decomposed = TFHEpp::LargeLUTRadixDecompose<W, K>(largelut_input);
    for (uint32_t i = 0; i < L; i++)
        TFHEpp::tlweSymIntEncrypt<domainP, plain_modulus>(
            largelut_digits[i], decomposed[i], 0.0, sk->key.get<domainP>());

    std::array<TFHEpp::TLWE<targetP>, 2> aes_nibbles;
    TFHEpp::tlweSymIntEncrypt<targetP, plain_modulus>(
        aes_nibbles[0], input & 0xF, 0.0, sk->key.get<targetP>());
    TFHEpp::tlweSymIntEncrypt<targetP, plain_modulus>(
        aes_nibbles[1], input >> 4, 0.0, sk->key.get<targetP>());

    std::array<TFHEpp::TLWE<targetP>, 8> aes_bits;
    for (uint32_t i = 0; i < aes_bits.size(); i++) {
        const auto message =
            ((input >> i) & 0x1)
                ? (typename targetP::T{1}
                   << (std::numeric_limits<typename targetP::T>::digits - 2))
                : (typename targetP::T{} -
                   (typename targetP::T{1}
                    << (std::numeric_limits<typename targetP::T>::digits - 2)));
        TFHEpp::tlweSymEncrypt<targetP>(
            aes_bits[i], message, 0.0, sk->key.get<targetP>());
    }

    const auto low_table = LargeLUTNibbleTRLWE<targetP, W, K>(false);
    const auto high_table = LargeLUTNibbleTRLWE<targetP, W, K>(true);

    std::array<TFHEpp::TLWE<domainP>, L> calibrated;
    TFHEpp::LargeLUTCalibrate<brP, iksP, W, K>(
        calibrated, largelut_digits, ek.getbkfft<brP>(), ek.getiksk<iksP>());

    auto largelut_from_calibrated = [&] {
        std::array<TFHEpp::TLWE<targetP>, 2> res;
        TFHEpp::LargeLUTWithCalibrated<brP, ahP, W, K>(
            res[0], low_table, calibrated, ek.getbkfft<brP>(),
            ek.getahk<ahP>());
        TFHEpp::LargeLUTWithCalibrated<brP, ahP, W, K>(
            res[1], high_table, calibrated, ek.getbkfft<brP>(),
            ek.getahk<ahP>());
        return res;
    };

    auto largelut_r2r_from_calibrated = [&] {
        std::array<TFHEpp::TLWE<targetP>, 2> res;
        TFHEpp::LargeLUTOptimizedWithCalibrated<brP, r2r_t, r2r_basebit, W, K>(
            res[0], low_table, calibrated, ek.getbkfft<brP>(),
            *step_r2rks);
        TFHEpp::LargeLUTOptimizedWithCalibrated<brP, r2r_t, r2r_basebit, W, K>(
            res[1], high_table, calibrated, ek.getbkfft<brP>(),
            *step_r2rks);
        return res;
    };

    auto largelut_full = [&] {
        std::array<TFHEpp::TLWE<domainP>, L> local_calibrated;
        TFHEpp::LargeLUTCalibrate<brP, iksP, W, K>(
            local_calibrated, largelut_digits, ek.getbkfft<brP>(),
            ek.getiksk<iksP>());

        std::array<TFHEpp::TLWE<targetP>, 2> res;
        TFHEpp::LargeLUTWithCalibrated<brP, ahP, W, K>(
            res[0], low_table, local_calibrated, ek.getbkfft<brP>(),
            ek.getahk<ahP>());
        TFHEpp::LargeLUTWithCalibrated<brP, ahP, W, K>(
            res[1], high_table, local_calibrated, ek.getbkfft<brP>(),
            ek.getahk<ahP>());
        return res;
    };

    auto largelut_r2r_full = [&] {
        std::array<TFHEpp::TLWE<domainP>, L> local_calibrated;
        TFHEpp::LargeLUTCalibrate<brP, iksP, W, K>(
            local_calibrated, largelut_digits, ek.getbkfft<brP>(),
            ek.getiksk<iksP>());

        std::array<TFHEpp::TLWE<targetP>, 2> res;
        TFHEpp::LargeLUTOptimizedWithCalibrated<brP, r2r_t, r2r_basebit, W, K>(
            res[0], low_table, local_calibrated, ek.getbkfft<brP>(),
            *step_r2rks);
        TFHEpp::LargeLUTOptimizedWithCalibrated<brP, r2r_t, r2r_basebit, W, K>(
            res[1], high_table, local_calibrated, ek.getbkfft<brP>(),
            *step_r2rks);
        return res;
    };

    auto aes_inv_sbox = [&] {
        std::array<TFHEpp::TLWE<targetP>, 2> res;
        TFHEpp::AESInvSbox<iksP, brP, ahP>(res, aes_nibbles, ek);
        return res;
    };

    auto aes_inv_sbox_rom = [&] {
        std::array<TFHEpp::TLWE<targetP>, 8> res;
        TFHEpp::AESInvSboxROM<iksP, brP, ahP>(res, aes_bits, ek);
        return res;
    };

    TFHEpp::TLWE<domainP> shifted_lower;
    TFHEpp::IdentityKeySwitch<iksP>(shifted_lower, aes_nibbles[0],
                                   ek.getiksk<iksP>());
    shifted_lower[domainP::k * domainP::n] +=
        1ULL << (std::numeric_limits<typename domainP::T>::digits - 6);

    TFHEpp::TLWE<domainP> shifted_upper;
    TFHEpp::IdentityKeySwitch<iksP>(shifted_upper, aes_nibbles[1],
                                   ek.getiksk<iksP>());
    shifted_upper[domainP::k * domainP::n] +=
        1ULL << (std::numeric_limits<typename domainP::T>::digits - 6);

    std::array<std::array<TFHEpp::TLWE<targetP>, 2>, 16> midtlwes;
    for (uint32_t i = 0; i < 16; i++)
        TFHEpp::GateBootstrappingManyLUT<brP, 2>(
            midtlwes[i], shifted_lower, ek.getbkfft<brP>(),
            TFHEpp::AESInvSboxPoly<targetP>(i));

    auto cmux_one_br_stage = [&] {
        std::array<std::array<TFHEpp::TLWE<targetP>, 16>, 2> tabletlwe;
        for (uint32_t i = 0; i < 2; i++)
            for (uint32_t j = 0; j < 16; j++) tabletlwe[i][j] = midtlwes[j][i];

        TFHEpp::TRLWE<targetP> trlwe;
        TFHEpp::TLWE2TablePackingManyLUT<ahP, 16, 2>(
            trlwe, tabletlwe, ek.getahk<ahP>());

        std::array<TFHEpp::TLWE<targetP>, 2> res;
        TFHEpp::GateBootstrappingManyLUT<brP, 2>(
            res, shifted_upper, ek.getbkfft<brP>(), trlwe);
        return res;
    };

    const auto expected = inv_sbox[input >> 4][input & 0xF];
    const auto large_check = largelut_full();
    const auto r2r_check = largelut_r2r_full();
    const auto cmux_check = aes_inv_sbox();
    const auto cmux_rom_check = aes_inv_sbox_rom();
    const auto cmux_stage_check = cmux_one_br_stage();
    const uint8_t large_decrypted =
        decrypt_nibble_pair<targetP>(large_check, *sk);
    const uint8_t r2r_decrypted = decrypt_nibble_pair<targetP>(r2r_check, *sk);
    const uint8_t cmux_decrypted =
        decrypt_nibble_pair<targetP>(cmux_check, *sk);
    const uint8_t cmux_rom_decrypted =
        decrypt_bits<targetP>(cmux_rom_check, *sk);
    const uint8_t cmux_stage_decrypted =
        decrypt_nibble_pair<targetP>(cmux_stage_check, *sk);
    assert(cmux_decrypted == expected);
    assert(cmux_rom_decrypted == expected);
    assert(cmux_stage_decrypted == expected);

    std::vector<double> largelut_ah_full_ms;
    std::vector<double> largelut_ah_calibrated_ms;
    std::vector<double> largelut_r2r_full_ms;
    std::vector<double> largelut_r2r_calibrated_ms;
    std::vector<double> aes_inv_sbox_ms;
    std::vector<double> aes_inv_sbox_rom_ms;
    std::vector<double> cmux_one_br_ms;

    uint64_t sink = checksum<targetP>(large_check) + checksum<targetP>(r2r_check) +
                    checksum<targetP>(cmux_check) +
                    checksum<targetP>(cmux_rom_check) +
                    checksum<targetP>(cmux_stage_check);

    for (uint32_t repeat = 0; repeat < repeats; repeat++) {
        std::array<TFHEpp::TLWE<targetP>, 2> local;

        largelut_ah_full_ms.push_back(time_ms([&] {
            local = largelut_full();
            sink += checksum<targetP>(local);
        }));

        largelut_ah_calibrated_ms.push_back(time_ms([&] {
            local = largelut_from_calibrated();
            sink += checksum<targetP>(local);
        }));

        largelut_r2r_full_ms.push_back(time_ms([&] {
            local = largelut_r2r_full();
            sink += checksum<targetP>(local);
        }));

        largelut_r2r_calibrated_ms.push_back(time_ms([&] {
            local = largelut_r2r_from_calibrated();
            sink += checksum<targetP>(local);
        }));

        aes_inv_sbox_ms.push_back(time_ms([&] {
            local = aes_inv_sbox();
            sink += checksum<targetP>(local);
        }));

        aes_inv_sbox_rom_ms.push_back(time_ms([&] {
            const auto bits = aes_inv_sbox_rom();
            sink += checksum<targetP>(bits);
        }));

        cmux_one_br_ms.push_back(time_ms([&] {
            local = cmux_one_br_stage();
            sink += checksum<targetP>(local);
        }));
    }

    const double large_ah_full = average(largelut_ah_full_ms);
    const double large_ah_calibrated = average(largelut_ah_calibrated_ms);
    const double large_r2r_full = average(largelut_r2r_full_ms);
    const double large_r2r_calibrated =
        average(largelut_r2r_calibrated_ms);
    const double aes_full = average(aes_inv_sbox_ms);
    const double aes_rom_full = average(aes_inv_sbox_rom_ms);
    const double cmux_stage = average(cmux_one_br_ms);
    const double r2r_keys_mib =
        static_cast<double>((L - 1) *
                            R2RKeyBytes<targetP, r2r_t, r2r_basebit>()) /
        (1024.0 * 1024.0);

    std::cout << std::fixed << std::setprecision(3);
    std::cout << "LargeLUT vs AES-style CMUX benchmark, brP=lvlh2param, "
              << "iksP=lvl2hparam, ahP=AHlvl2param, R2R t=" << r2r_t
              << ", basebit=" << r2r_basebit << "\n";
    std::cout << "AES inverse S-box as two 4-bit LUT outputs embedded with "
              << (uint32_t{1} << (K - 8)) << " adjacent copies per byte, W="
              << W << ", K=" << K << ", digits=" << L << "\n";
    std::cout << "input: 0x" << std::hex << static_cast<uint32_t>(input)
              << ", expected: 0x" << static_cast<uint32_t>(expected)
              << std::dec << "\n";
    std::cout << "AH decrypted: 0x" << std::hex
              << static_cast<uint32_t>(large_decrypted) << std::dec
              << (large_decrypted == expected ? " (ok)" : " (FAIL)")
              << ", max error/Delta: "
              << static_cast<double>(max_nibble_pair_error_units<targetP>(
                     large_check, expected, *sk))
              << "\n";
    std::cout << "R2R decrypted: 0x" << std::hex
              << static_cast<uint32_t>(r2r_decrypted) << std::dec
              << (r2r_decrypted == expected ? " (ok)" : " (FAIL)")
              << ", max error/Delta: "
              << static_cast<double>(max_nibble_pair_error_units<targetP>(
                     r2r_check, expected, *sk))
              << "\n";
    std::cout << "AESInvSbox decrypted: 0x" << std::hex
              << static_cast<uint32_t>(cmux_decrypted) << std::dec
              << " (ok), max error/Delta: "
              << static_cast<double>(max_nibble_pair_error_units<targetP>(
                     cmux_check, expected, *sk))
              << "\n";
    std::cout << "AESInvSboxROM+CB decrypted: 0x" << std::hex
              << static_cast<uint32_t>(cmux_rom_decrypted) << std::dec
              << " (ok)\n";
    std::cout << "CMUX+one BR decrypted: 0x" << std::hex
              << static_cast<uint32_t>(cmux_stage_decrypted) << std::dec
              << " (ok), max error/Delta: "
              << static_cast<double>(max_nibble_pair_error_units<targetP>(
                     cmux_stage_check, expected, *sk))
              << "\n";
    std::cout << "BootstrappingKeyFFT gen ms: " << bk_keygen_ms << "\n";
    std::cout << "KeySwitchingKey gen ms: " << iksk_keygen_ms << "\n";
    std::cout << "AnnihilateKey gen ms: " << ah_keygen_ms << "\n";
    std::cout << "CBswitchingKey gen ms: " << cbsk_keygen_ms << "\n";
    std::cout << "R2R step keys MiB: " << r2r_keys_mib << "\n";
    std::cout << "R2R step key gen ms: " << r2r_keygen_ms << "\n";
    std::cout << "LargeLUT AH full shared-calibration eval ms: " << large_ah_full
              << "\n";
    std::cout << "LargeLUT AH after calibration, two outputs ms: "
              << large_ah_calibrated << "\n";
    std::cout << "LargeLUT R2R full shared-calibration eval ms: "
              << large_r2r_full << "\n";
    std::cout << "LargeLUT R2R after calibration, two outputs ms: "
              << large_r2r_calibrated << "\n";
    std::cout << "AESInvSbox full CMUX/LUT eval ms: " << aes_full << "\n";
    std::cout << "AESInvSboxROM+CB full eval ms: " << aes_rom_full << "\n";
    std::cout << "AES-style TLWE2TablePackingManyLUT+one BR ms: " << cmux_stage
              << "\n";
    std::cout << "AH full LargeLUT speedup vs AESInvSbox: "
              << aes_full / large_ah_full << "x\n";
    std::cout << "R2R full LargeLUT speedup vs AESInvSbox: "
              << aes_full / large_r2r_full << "x\n";
    std::cout << "AH full LargeLUT speedup vs AESInvSboxROM+CB: "
              << aes_rom_full / large_ah_full << "x\n";
    std::cout << "R2R full LargeLUT speedup vs AESInvSboxROM+CB: "
              << aes_rom_full / large_r2r_full << "x\n";
    std::cout << "AH post-calibration speedup vs CMUX+one BR stage: "
              << cmux_stage / large_ah_calibrated << "x\n";
    std::cout << "R2R post-calibration speedup vs CMUX+one BR stage: "
              << cmux_stage / large_r2r_calibrated << "x\n";
    std::cout << "R2R/AH full runtime ratio: "
              << large_r2r_full / large_ah_full << "x\n";
    std::cout << "R2R/AH post-calibration runtime ratio: "
              << large_r2r_calibrated / large_ah_calibrated << "x\n";
    std::cout << "checksum: " << sink << "\n";
}
