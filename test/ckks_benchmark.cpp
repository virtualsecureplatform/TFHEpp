#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <memory>
#include <tfhe++.hpp>

namespace {

using P = TFHEpp::ckkslvl3param;
using Clock = std::chrono::steady_clock;

void fill_key(TFHEpp::Key<P> &key)
{
    for (std::size_t i = 0; i < P::n; i++) {
        const int v = static_cast<int>(i % 3) - 1;
        key[i] = static_cast<typename P::T>(v);
    }
}

template <class F>
double elapsed_ms(F &&fn)
{
    const auto begin = Clock::now();
    fn();
    const auto end = Clock::now();
    return std::chrono::duration<double, std::milli>(end - begin).count();
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
    constexpr std::uint32_t lhs_log_delta = 50;
    constexpr std::uint32_t rhs_log_delta = 50;
    constexpr std::uint32_t log_q = 128;
    const TFHEpp::CKKSNoise relin_noise{P::α, 0};
    using LhsCt = TFHEpp::CKKSCiphertext<P, log_q, lhs_log_delta>;
    using RhsCt = TFHEpp::CKKSCiphertext<P, log_q, rhs_log_delta>;
    using ProductCt =
        TFHEpp::CKKSMultResult<P, log_q, lhs_log_delta, log_q, rhs_log_delta>;

    auto key = std::make_unique<TFHEpp::Key<P>>();
    const double key_ms = elapsed_ms([&] { fill_key(*key); });

    std::unique_ptr<TFHEpp::CKKSRelinKey<P, ProductCt::log_q>> relinkey;
    const double relinkey_ms = elapsed_ms([&] {
        relinkey =
            TFHEpp::makeCKKSRelinKey<P, ProductCt::log_q>(*key, relin_noise);
    });

    auto lhs = std::make_unique<std::array<double, P::n>>();
    auto rhs = std::make_unique<std::array<double, P::n>>();
    auto expected = std::make_unique<std::array<double, P::n>>();
    auto got = std::make_unique<std::array<double, P::n>>();
    lhs->fill(0.0);
    rhs->fill(0.0);

    (*lhs)[0] = 0.25;
    (*lhs)[1] = -0.125;
    (*lhs)[3] = 0.0625;
    (*rhs)[0] = -0.5;
    (*rhs)[2] = 0.25;
    negacyclic_product(*expected, *lhs, *rhs);

    auto lhs_ct = std::make_unique<LhsCt>();
    auto rhs_ct = std::make_unique<RhsCt>();
    const double encrypt_ms = elapsed_ms([&] {
        TFHEpp::ckksEncrypt<P, log_q, lhs_log_delta>(*lhs_ct, *lhs, *key);
        TFHEpp::ckksEncrypt<P, log_q, rhs_log_delta>(*rhs_ct, *rhs, *key);
    });

    auto tensor = std::make_unique<TFHEpp::TRLWE3<P>>();
    const double tensor_ms = elapsed_ms([&] {
        TFHEpp::CKKSTensorProductRescale<P, log_q, log_q, lhs_log_delta>(
            *tensor, lhs_ct->ct, rhs_ct->ct);
    });

    auto relin_product = std::make_unique<ProductCt>();
    const double relin_ms = elapsed_ms([&] {
        TFHEpp::CKKSRelinearization<P, ProductCt::log_q>(
            relin_product->ct, *tensor, *relinkey);
    });

    auto product = std::make_unique<ProductCt>();
    const double multiply_ms = elapsed_ms(
        [&] { TFHEpp::CKKSMult<P>(*product, *lhs_ct, *rhs_ct, *relinkey); });

    const double decrypt_ms =
        elapsed_ms([&] { TFHEpp::ckksDecrypt<P>(*got, *product, *key); });
    const double error = max_abs_error(*got, *expected);

    std::cout << "CKKS benchmark n=" << P::n << " lbar=" << P::l̅
              << " Bbarbit=" << P::B̅gbit << '\n';
    std::cout << "key_fill_ms=" << key_ms << '\n';
    std::cout << "relin_keygen_ms=" << relinkey_ms << '\n';
    std::cout << "encrypt_pair_ms=" << encrypt_ms << '\n';
    std::cout << "tensor_rescale_ms=" << tensor_ms << '\n';
    std::cout << "relinearization_ms=" << relin_ms << '\n';
    std::cout << "multiply_total_ms=" << multiply_ms << '\n';
    std::cout << "decrypt_ms=" << decrypt_ms << '\n';
    std::cout << "max_abs_error=" << error << std::endl;

    return error < 0.25 ? 0 : 1;
}
