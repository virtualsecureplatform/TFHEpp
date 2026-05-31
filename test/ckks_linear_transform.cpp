#include <algorithm>
#include <cmath>
#include <complex>
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

double max_error(
    const std::array<std::complex<double>, P::n / 2> &got,
    const std::array<std::complex<double>, P::n / 2> &want)
{
    double err = 0.0;
    for (std::size_t i = 0; i < P::n / 2; i++)
        err = std::max(err, std::abs(got[i] - want[i]));
    return err;
}

void require_close(
    const std::array<std::complex<double>, P::n / 2> &got,
    const std::array<std::complex<double>, P::n / 2> &want, double tol,
    const char *label)
{
    const double err = max_error(got, want);
    std::cout << label << " max_error=" << err << std::endl;
    if (err > tol) {
        std::cerr << label << " exceeded tolerance " << tol << std::endl;
        std::exit(1);
    }
}

}  // namespace

int main()
{
    constexpr std::uint32_t log_q = 110;
    constexpr std::uint32_t log_delta = 40;
    constexpr std::uint32_t plain_log_delta = 20;
    constexpr double tol = 0.05;
    using InCt = TFHEpp::CKKSCiphertext<P, log_q, log_delta>;
    using OutCt =
        TFHEpp::CKKSPlainMulResult<P, log_q, log_delta, plain_log_delta>;
    static_assert(OutCt::log_q == 90);
    static_assert(OutCt::log_delta == log_delta);

    auto key = std::make_unique<TFHEpp::Key<P>>();
    fill_key(*key);

    auto slots = std::make_unique<TFHEpp::CKKSSlotVector<P>>();
    auto rotated_expected = std::make_unique<TFHEpp::CKKSSlotVector<P>>();
    auto decoded = std::make_unique<TFHEpp::CKKSSlotVector<P>>();
    slots->fill({0.0, 0.0});
    for (std::size_t i = 0; i < 64; i++) {
        const double re = static_cast<double>(static_cast<int>(i % 11) - 5) /
                          64.0;
        const double im = static_cast<double>(static_cast<int>(i % 7) - 3) /
                          128.0;
        (*slots)[i] = {re, im};
    }
    (*slots)[P::n / 2 - 2] = {-0.0625, 0.03125};
    (*slots)[P::n / 2 - 1] = {0.125, -0.0625};

    auto ct = std::make_unique<InCt>();
    TFHEpp::ckksSlotEncrypt<P, log_q, log_delta>(*ct, *slots, *key, {0.0, 0});

    auto input_gk = std::make_unique<TFHEpp::CKKSGaloisKey<P, log_q>>();
    TFHEpp::CKKSGaloisKeyGen<P, log_q>(*input_gk, *key, {0.0, 0});
    auto output_gk = std::make_unique<TFHEpp::CKKSGaloisKey<P, OutCt::log_q>>();
    TFHEpp::CKKSGaloisKeyGen<P, OutCt::log_q>(*output_gk, *key, {0.0, 0});

    {
        constexpr int steps = 3;
        constexpr int half = static_cast<int>(P::n) / 2;
        auto rotated_ct = std::make_unique<InCt>();
        TFHEpp::CKKSRotateSlots<P, log_q>(rotated_ct->ct, ct->ct, steps,
                                          *input_gk);
        TFHEpp::ckksSlotDecrypt<P, log_q, log_delta>(*decoded, *rotated_ct,
                                                     *key);
        for (int i = 0; i < half; i++)
            (*rotated_expected)[i] = (*slots)[(i + steps) % half];
        require_close(*decoded, *rotated_expected, tol, "CKKS rotation");
    }

    std::vector<TFHEpp::CKKSSlotVector<P>> diagonals(3);
    for (auto &d : diagonals) d.fill({0.0, 0.0});
    for (std::size_t i = 0; i < P::n / 2; i++) {
        diagonals[0][i] = {0.5, 0.0};
        diagonals[1][i] = {0.125, -0.0625};
        diagonals[2][i] = {i % 2 == 0 ? -0.25 : 0.25, 0.03125};
    }
    const std::vector<int> offsets{0, 1, -2};

    auto expected = std::make_unique<TFHEpp::CKKSSlotVector<P>>();
    constexpr int half = static_cast<int>(P::n) / 2;
    for (int i = 0; i < half; i++) {
        (*expected)[i] =
            diagonals[0][i] * (*slots)[i] +
            diagonals[1][i] * (*slots)[(i + 1) % half] +
            diagonals[2][i] * (*slots)[(i - 2 + half) % half];
    }

    auto transformed = std::make_unique<OutCt>();
    TFHEpp::CKKSLinearTransformBSGS<P, log_q, log_delta, plain_log_delta>(
        *transformed, *ct, diagonals, offsets, 2, *input_gk, *output_gk);
    TFHEpp::ckksSlotDecrypt<P, OutCt::log_q, OutCt::log_delta>(
        *decoded, *transformed, *key);

    require_close(*decoded, *expected, tol, "CKKS linear transform");
    std::cout << "Passed" << std::endl;
}
