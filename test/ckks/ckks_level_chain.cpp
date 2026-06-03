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

}  // namespace

int main()
{
    constexpr std::uint32_t log_q = 128;
    constexpr std::uint32_t log_delta = 36;
    const TFHEpp::CKKSNoise level_noise{P::α, 0};

    using FreshCt = TFHEpp::CKKSCiphertext<P, log_q, log_delta>;
    using Level1Ct =
        TFHEpp::CKKSMultResult<P, FreshCt::log_q, FreshCt::log_delta,
                               FreshCt::log_q, FreshCt::log_delta>;
    using Level2Ct =
        TFHEpp::CKKSMultResult<P, Level1Ct::log_q, Level1Ct::log_delta,
                               FreshCt::log_q, FreshCt::log_delta>;

    static_assert(Level1Ct::log_q == 92);
    static_assert(Level1Ct::log_delta == log_delta);
    static_assert(Level2Ct::log_q == 56);
    static_assert(Level2Ct::log_delta == log_delta);

    auto key = std::make_unique<TFHEpp::Key<P>>();
    fill_key(*key);
    auto relin_level1 =
        TFHEpp::makeCKKSRelinKey<P, Level1Ct::log_q>(*key, level_noise);
    auto relin_level2 =
        TFHEpp::makeCKKSRelinKey<P, Level2Ct::log_q>(*key, level_noise);

    auto lhs = std::make_unique<std::array<double, P::n>>();
    auto rhs = std::make_unique<std::array<double, P::n>>();
    auto extra = std::make_unique<std::array<double, P::n>>();
    auto expected_level1 = std::make_unique<std::array<double, P::n>>();
    auto expected_level2 = std::make_unique<std::array<double, P::n>>();
    auto got = std::make_unique<std::array<double, P::n>>();
    lhs->fill(0.0);
    rhs->fill(0.0);
    extra->fill(0.0);

    (*lhs)[0] = 0.125;
    (*lhs)[1] = -0.0625;
    (*lhs)[4] = 0.03125;
    (*rhs)[0] = -0.25;
    (*rhs)[2] = 0.125;
    (*rhs)[5] = 0.0625;
    (*extra)[0] = 0.5;
    (*extra)[3] = -0.25;

    negacyclic_product(*expected_level1, *lhs, *rhs);
    negacyclic_product(*expected_level2, *expected_level1, *extra);

    auto lhs_ct = std::make_unique<FreshCt>();
    auto rhs_ct = std::make_unique<FreshCt>();
    auto extra_ct = std::make_unique<FreshCt>();
    auto level1_ct = std::make_unique<Level1Ct>();
    auto level2_ct = std::make_unique<Level2Ct>();

    TFHEpp::ckksEncrypt<P, FreshCt::log_q, FreshCt::log_delta>(
        *lhs_ct, *lhs, *key, level_noise);
    TFHEpp::ckksEncrypt<P, FreshCt::log_q, FreshCt::log_delta>(
        *rhs_ct, *rhs, *key, level_noise);
    TFHEpp::ckksEncrypt<P, FreshCt::log_q, FreshCt::log_delta>(
        *extra_ct, *extra, *key, level_noise);

    TFHEpp::CKKSMult<P>(*level1_ct, *lhs_ct, *rhs_ct, *relin_level1);
    TFHEpp::CKKSMult<P>(*level2_ct, *level1_ct, *extra_ct, *relin_level2);

    TFHEpp::ckksDecrypt<P>(*got, *level2_ct, *key);
    const double error = max_abs_error(*got, *expected_level2);

    std::cout << "CKKS chain levels " << FreshCt::log_q << " -> "
              << Level1Ct::log_q << " -> " << Level2Ct::log_q
              << " at delta=" << log_delta << '\n';
    std::cout << "max_abs_error=" << error << std::endl;

    if (error > 0.05) {
        std::cerr << "CKKS level chain exceeded tolerance" << std::endl;
        return 1;
    }
    std::cout << "Passed" << std::endl;
}
