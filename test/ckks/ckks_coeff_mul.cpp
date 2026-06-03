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

void require_close(double got, double want, double tol, const char *label)
{
    if (std::abs(got - want) > tol) {
        std::cerr << label << " got=" << got << " want=" << want
                  << " tol=" << tol << std::endl;
        std::exit(1);
    }
}

}  // namespace

int main()
{
    constexpr std::uint32_t lhs_log_delta = 42;
    constexpr std::uint32_t rhs_log_delta = 40;
    constexpr std::uint32_t log_q = 88;
    constexpr double tol = 0.25;
    using LhsCt = TFHEpp::CKKSCiphertext<P, log_q, lhs_log_delta>;
    using RhsCt = TFHEpp::CKKSCiphertext<P, log_q, rhs_log_delta>;
    using ProductCt =
        TFHEpp::CKKSMultResult<P, log_q, lhs_log_delta, log_q, rhs_log_delta>;
    static_assert(ProductCt::log_q == log_q - lhs_log_delta);
    static_assert(ProductCt::log_delta == rhs_log_delta);

    auto key = std::make_unique<TFHEpp::Key<P>>();
    fill_key(*key);
    auto relinkey =
        TFHEpp::makeCKKSRelinKey<P, ProductCt::log_q>(*key, {0.0, 0});

    auto lhs = std::make_unique<std::array<double, P::n>>();
    auto rhs = std::make_unique<std::array<double, P::n>>();
    auto expected = std::make_unique<std::array<double, P::n>>();
    auto got = std::make_unique<std::array<double, P::n>>();
    lhs->fill(0.0);
    rhs->fill(0.0);
    expected->fill(0.0);

    (*lhs)[0] = 0.25;
    (*lhs)[1] = -0.125;
    (*lhs)[3] = 0.0625;
    (*rhs)[0] = -0.5;
    (*rhs)[2] = 0.25;

    add_product(*expected, 0, (*lhs)[0], 0, (*rhs)[0]);
    add_product(*expected, 0, (*lhs)[0], 2, (*rhs)[2]);
    add_product(*expected, 1, (*lhs)[1], 0, (*rhs)[0]);
    add_product(*expected, 1, (*lhs)[1], 2, (*rhs)[2]);
    add_product(*expected, 3, (*lhs)[3], 0, (*rhs)[0]);
    add_product(*expected, 3, (*lhs)[3], 2, (*rhs)[2]);

    auto lhs_ct = std::make_unique<LhsCt>();
    auto rhs_ct = std::make_unique<RhsCt>();
    auto product = std::make_unique<ProductCt>();
    TFHEpp::ckksEncrypt<P, log_q, lhs_log_delta>(*lhs_ct, *lhs, *key);
    TFHEpp::ckksEncrypt<P, log_q, rhs_log_delta>(*rhs_ct, *rhs, *key);

    TFHEpp::CKKSMult<P>(*product, *lhs_ct, *rhs_ct, *relinkey);

    TFHEpp::ckksDecrypt<P>(*got, *product, *key);
    for (std::size_t i = 0; i < 8; i++)
        require_close((*got)[i], (*expected)[i], tol, "CKKS coefficient multiply");
    for (std::size_t i = 8; i < P::n; i++)
        require_close((*got)[i], 0.0, tol, "CKKS coefficient multiply zero tail");

    std::cout << "Passed" << std::endl;
}
