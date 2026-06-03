#include <cassert>
#include <cstdint>
#include <iostream>
#include <limits>
#include <memory>
#include <vector>
#include <tfhe++.hpp>

struct TinyR2RParam {
    static constexpr int32_t key_value_max = 1;
    static constexpr int32_t key_value_min = 0;
    static constexpr uint32_t nbit = 3;
    static constexpr uint32_t n = 1U << nbit;
    static constexpr uint32_t k = 1;
    static constexpr TFHEpp::ErrorDistribution errordist =
        TFHEpp::ErrorDistribution::ModularGaussian;
    static constexpr double α = 0.0;
    using T = uint16_t;
    static constexpr uint32_t plain_modulus = 8;
    static constexpr double Δ =
        static_cast<double>(uint64_t{1} << std::numeric_limits<T>::digits) /
        plain_modulus;
};

void test_r2rpks()
{
    using P = TinyR2RParam;
    constexpr uint32_t plain_modulus = P::plain_modulus;
    constexpr uint32_t t = 4;
    constexpr uint32_t basebit = 4;

    TFHEpp::Key<P> key{};
    for (uint32_t i = 0; i < P::n; i++) key[i] = i & 1U;

    TFHEpp::Polynomial<P> plain{};
    for (uint32_t i = 0; i < P::n; i++) plain[i] = (3 * i + 1) & 7U;

    TFHEpp::TRLWE<P> input;
    TFHEpp::trlweSymIntEncrypt<P>(input, plain, key);

    const std::vector<uint32_t> positions{1, 6};
    constexpr uint32_t R = 2;
    TFHEpp::R2RKey<P, t, basebit> r2rk;
    TFHEpp::R2RKeyGen<P, t, basebit>(r2rk, positions, R, key);

    TFHEpp::TRLWE<P> output;
    TFHEpp::R2RPKS<P, t, basebit>(output, input, positions, R, r2rk);

    const auto decrypted =
        TFHEpp::trlweSymIntDecrypt<P, plain_modulus>(output, key);
    assert(decrypted[0] == plain[positions[0]]);
    assert(decrypted[1] == plain[positions[0]]);
    assert(decrypted[2] == plain[positions[1]]);
    assert(decrypted[3] == plain[positions[1]]);
    for (uint32_t i = 4; i < P::n; i++) assert(decrypted[i] == 0);
}

void instantiate_optimized_large_lut_api(
    TFHEpp::TLWE<TFHEpp::lvl1param> &result,
    const TFHEpp::TRLWE<TFHEpp::lvl1param> &table,
    const std::array<TFHEpp::TLWE<TFHEpp::lvl0param>,
                     TFHEpp::LargeLUTDigitCount<2, 3>()> &digits,
    const TFHEpp::BootstrappingKeyFFT<TFHEpp::lvl01param> &bkfft,
    const TFHEpp::KeySwitchingKey<TFHEpp::lvl10param> &iksk,
    const TFHEpp::R2RKey<TFHEpp::lvl1param, 1, 1> &batch_r2rk,
    const std::array<TFHEpp::R2RKey<TFHEpp::lvl1param, 1, 1>,
                     TFHEpp::LargeLUTDigitCount<2, 3>() - 1> &step_r2rks)
{
    TFHEpp::LargeLUTOptimized<TFHEpp::lvl01param, TFHEpp::lvl10param, 1, 1, 2,
                              3>(result, table, digits, bkfft, iksk,
                                 batch_r2rk, step_r2rks);
}

void instantiate_large_lut_r2r_keygen(
    TFHEpp::R2RKey<TFHEpp::lvl1param, 1, 1> &batch_r2rk,
    std::array<TFHEpp::R2RKey<TFHEpp::lvl1param, 1, 1>,
               TFHEpp::LargeLUTDigitCount<2, 3>() - 1> &step_r2rks,
    const TFHEpp::Key<TFHEpp::lvl1param> &key)
{
    TFHEpp::LargeLUTR2RKeyGen<TFHEpp::lvl1param, 1, 1, 2, 3>(
        batch_r2rk, step_r2rks, key);
}

int main()
{
    test_r2rpks();

    using brP = TFHEpp::lvl01param;
    using iksP = TFHEpp::lvl10param;
    using ahP = typename brP::targetP;
    using domainP = typename brP::domainP;
    using targetP = typename brP::targetP;

    constexpr uint32_t W = 2;
    constexpr uint32_t K = 3;
    constexpr uint32_t message_modulus = 1U << (W + 1);
    constexpr uint32_t table_size = 1U << K;
    constexpr uint32_t digits_count = TFHEpp::LargeLUTDigitCount<W, K>();

    TFHEpp::SecretKey sk;
    auto bkfft = std::make_unique<TFHEpp::BootstrappingKeyFFT<brP>>();
    TFHEpp::bkfftgen<brP>(*bkfft, sk);
    auto iksk = std::make_unique<TFHEpp::KeySwitchingKey<iksP>>();
    TFHEpp::ikskgen<iksP>(*iksk, sk);
    auto ahk = std::make_unique<TFHEpp::AnnihilateKey<ahP>>();
    TFHEpp::annihilatekeygen<ahP>(*ahk, sk);

    std::array<uint32_t, table_size> table{};
    for (uint32_t i = 0; i < table_size; i++) table[i] = (3 * i + 1) & 3U;

    const auto table_poly = TFHEpp::LargeLUTPolynomial<targetP, W, K>(table);
    TFHEpp::TRLWE<targetP> table_trlwe;
    TFHEpp::TrivialTRLWEFromPolynomial<targetP>(table_trlwe, table_poly);

    for (uint32_t x = 0; x < table_size; x++) {
        const auto decomposed = TFHEpp::LargeLUTRadixDecompose<W, K>(x);
        std::array<TFHEpp::TLWE<domainP>, digits_count> encrypted_digits;
        for (uint32_t i = 0; i < digits_count; i++)
            TFHEpp::tlweSymIntEncrypt<domainP, message_modulus>(
                encrypted_digits[i], decomposed[i], sk.key.get<domainP>());

        TFHEpp::TLWE<targetP> result;
        TFHEpp::LargeLUT<brP, iksP, ahP, W, K>(
            result, table_trlwe, encrypted_digits, *bkfft, *iksk, *ahk);

        const auto decrypted =
            TFHEpp::tlweSymIntDecrypt<targetP, message_modulus>(
                result, sk.key.get<targetP>());
        if (decrypted != table[x]) {
            std::cerr << "LargeLUT mismatch at x=" << x << ": expected "
                      << table[x] << ", got " << decrypted << std::endl;
        }
        assert(decrypted == table[x]);
    }

    std::cout << "Passed" << std::endl;
}
