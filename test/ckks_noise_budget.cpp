#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <tfhe++.hpp>

namespace {

using P = TFHEpp::ckkslvl3param;

void fill_key(TFHEpp::Key<P> &key)
{
    for (std::size_t i = 0; i < P::n; i++) {
        const int v = static_cast<int>(i % 3) - 1;
        key[i] = static_cast<typename P::T>(v);
    }
}

void add_product(std::array<double, P::n> &out, std::size_t ai, double av,
                 std::size_t bi, double bv)
{
    const std::size_t raw = ai + bi;
    if (raw < P::n)
        out[raw] += av * bv;
    else
        out[raw - P::n] -= av * bv;
}

void negacyclic_product(std::array<double, P::n> &out,
                        const std::array<double, P::n> &lhs,
                        const std::array<double, P::n> &rhs)
{
    out.fill(0.0);
    for (std::size_t i = 0; i < P::n; i++) {
        if (lhs[i] == 0.0) continue;
        for (std::size_t j = 0; j < P::n; j++) {
            if (rhs[j] == 0.0) continue;
            add_product(out, i, lhs[i], j, rhs[j]);
        }
    }
}

double max_abs_error(const std::array<double, P::n> &got,
                     const std::array<double, P::n> &want)
{
    double max_error = 0.0;
    for (std::size_t i = 0; i < P::n; i++)
        max_error = std::max(max_error, std::abs(got[i] - want[i]));
    return max_error;
}

void require_error_below(const char *label,
                         const std::array<double, P::n> &got,
                         const std::array<double, P::n> &want, double tol)
{
    const double error = max_abs_error(got, want);
    std::cout << label << " max_abs_error=" << error << std::endl;
    if (error > tol) {
        std::cerr << label << " exceeded tolerance " << tol << std::endl;
        std::exit(1);
    }
}

}  // namespace

int main()
{
    constexpr std::uint32_t log_delta = 50;
    constexpr std::uint32_t log_q = 128;
    const TFHEpp::CKKSNoise level_noise{P::α, 0};
    using FreshCt = TFHEpp::CKKSCiphertext<P, log_q, log_delta>;
    using ProductCt =
        TFHEpp::CKKSMultResult<P, log_q, log_delta, log_q, log_delta>;
    static_assert(ProductCt::log_q == log_q - log_delta);
    static_assert(ProductCt::log_delta == log_delta);

    auto key = std::make_unique<TFHEpp::Key<P>>();
    fill_key(*key);
    auto relinkey =
        TFHEpp::makeCKKSRelinKey<P, ProductCt::log_q>(*key, level_noise);

    auto lhs = std::make_unique<std::array<double, P::n>>();
    auto rhs = std::make_unique<std::array<double, P::n>>();
    auto expected = std::make_unique<std::array<double, P::n>>();
    auto got = std::make_unique<std::array<double, P::n>>();
    lhs->fill(0.0);
    rhs->fill(0.0);

    (*lhs)[0] = 0.25;
    (*lhs)[1] = -0.125;
    (*lhs)[7] = 0.03125;
    (*rhs)[0] = -0.5;
    (*rhs)[2] = 0.25;
    (*rhs)[5] = -0.0625;
    negacyclic_product(*expected, *lhs, *rhs);

    auto lhs_ct = std::make_unique<FreshCt>();
    auto rhs_ct = std::make_unique<FreshCt>();
    auto product = std::make_unique<ProductCt>();

    TFHEpp::ckksEncrypt<P, log_q, log_delta>(*lhs_ct, *lhs, *key, level_noise);
    TFHEpp::ckksDecrypt<P>(*got, *lhs_ct, *key);
    require_error_below("fresh ciphertext", *got, *lhs, 1e-6);

    TFHEpp::ckksEncrypt<P, log_q, log_delta>(*rhs_ct, *rhs, *key, level_noise);
    TFHEpp::CKKSMult<P>(*product, *lhs_ct, *rhs_ct, *relinkey);
    TFHEpp::ckksDecrypt<P>(*got, *product, *key);
    require_error_below("one multiply with noisy relin key", *got, *expected,
                        0.25);

    std::cout << "Passed" << std::endl;
}
