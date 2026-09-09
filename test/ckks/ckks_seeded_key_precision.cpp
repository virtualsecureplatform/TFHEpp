#include <cmath>
#include <iostream>
#include <memory>
#include <string_view>
#include <tfhe++.hpp>

// Optional full-ring diagnostic for choosing test noise. In particular,
// half-open signed uniform noise is biased, unlike a centered Gaussian.
void rotation_noise()
{
    using P = TFHEpp::lvl6param;
    constexpr std::uint32_t q = 880, scale = 50;
    auto secret = std::make_unique<TFHEpp::Key<P>>();
    TFHEpp::keyGen<P>(*secret);
    auto slots = std::make_unique<TFHEpp::CKKSSlotVector<P>>();
    slots->fill({0.25, 0.0});
    auto input = std::make_unique<TFHEpp::CKKSCiphertext<P, q, scale>>();
    TFHEpp::ckksSlotEncrypt<P, q, scale>(*input, *slots, *secret, {0.0, 0});
    for (const auto noise : {TFHEpp::CKKSNoise{0.0, 0},
                             TFHEpp::CKKSNoise{0.0, 4},
                             TFHEpp::CKKSNoise{std::ldexp(3.2, -880), 0}}) {
        auto key = std::make_unique<TFHEpp::CKKSSeededAutoKey<P, q>>();
        TFHEpp::CKKSAutoKeyGen<P, q>(*key, 5, *secret, noise);
        auto result = std::make_unique<TFHEpp::CKKSCiphertext<P, q, scale>>();
        TFHEpp::CKKSEvalAuto<P, q>(result->ct, input->ct, 5, *key);
        TFHEpp::ckksSlotDecrypt<P, q, scale>(*slots, *result, *secret);
        double max_error = 0;
        for (const auto value : *slots)
            max_error = std::max(max_error, std::abs(value - 0.25));
        std::cout << "rotation noise uniform_bits=" << noise.uniform_bits
                  << " gaussian=" << (noise.modular_stdev != 0)
                  << " max_error=" << max_error << std::endl;
    }
}

int main(int argc, char **argv)
{
    using P = TFHEpp::lvl6param;
    constexpr std::uint32_t q = 880;
    auto key = std::make_unique<TFHEpp::Key<P>>();
    TFHEpp::keyGen<P>(*key);
    auto row = std::make_unique<TFHEpp::CKKSSeededKeySwitchRow<P, q>>();
    auto zero = std::make_unique<TFHEpp::Polynomial<P>>();
    auto mask = std::make_unique<TFHEpp::Polynomial<P>>();
    auto secret = std::make_unique<TFHEpp::Polynomial<P>>();
    auto product = std::make_unique<TFHEpp::Polynomial<P>>();
    for (std::uint32_t i = 0; i < P::n; i++) (*secret)[i] = (*key)[i];
    TFHEpp::ckks_detail::encryptPolynomialAtLevel<P, q>(
        *row, *zero, *key, {0.0, 0});
    TFHEpp::ckks_detail::fillSeededUniformPolynomialAtLevel<P, q>(
        *mask, row->mask_seeds[0]);
    TFHEpp::ckks_detail::exactActiveTorusMulBySmall<P, q>(*product, *mask, *secret);
    long double max_error = 0;
    for (std::uint32_t i = 0; i < P::n; i++) {
        const auto error = TFHEpp::ckks_detail::reduceToLevel<P, q>(row->body[i] - (*product)[i]);
        max_error = std::max(max_error, std::abs(TFHEpp::ckks_detail::levelToLongDouble<P, q>(error)));
    }
    std::cout << "Zero-noise seeded row max coefficient error=" << max_error << std::endl;
    if (max_error != 0) return 1;
    if (argc == 2 && std::string_view(argv[1]) == "--rotation-noise")
        rotation_noise();
    else if (argc != 1) return 1;
    std::cout << "Passed" << std::endl;
}
