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
    constexpr std::uint32_t log_delta = 42;
    constexpr std::uint32_t log_q = 84;
    constexpr double tol = 1e-3;
    using Ct = TFHEpp::CKKSCiphertext<P, log_q, log_delta>;
    static_assert(Ct::log_delta == log_delta);
    static_assert(Ct::log_q == log_q);
    static_assert(Ct::log_budget == log_q - log_delta);

    auto key = std::make_unique<TFHEpp::Key<P>>();
    fill_key(*key);

    auto input = std::make_unique<std::array<double, P::n>>();
    auto output = std::make_unique<std::array<double, P::n>>();
    input->fill(0.0);
    (*input)[0] = 0.25;
    (*input)[1] = -0.5;
    (*input)[2] = 0.75;
    (*input)[11] = -0.125;

    auto ct = std::make_unique<Ct>();
    TFHEpp::ckksEncrypt<P, log_q, log_delta>(*ct, *input, *key);
    TFHEpp::ckksDecrypt<P>(*output, *ct, *key);

    for (std::size_t i = 0; i < P::n; i++)
        require_close((*output)[i], (*input)[i], tol, "encrypt/decrypt");

    std::cout << "Passed" << std::endl;
}
