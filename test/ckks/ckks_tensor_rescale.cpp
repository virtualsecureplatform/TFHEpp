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

TFHEpp::Polynomial<P> trlwe3Phase(const TFHEpp::TRLWE3<P> &ct,
                                  const TFHEpp::Key<P> &key)
{
    TFHEpp::Polynomial<P> partkey{};
    for (std::size_t i = 0; i < P::n; i++) partkey[i] = key[i];

    TFHEpp::Polynomial<P> keysquare{};
    TFHEpp::PolyMul<P>(keysquare, partkey, partkey);

    TFHEpp::Polynomial<P> cross{};
    TFHEpp::PolyMul<P>(cross, ct[0], partkey);

    TFHEpp::Polynomial<P> square{};
    TFHEpp::PolyMul<P>(square, ct[2], keysquare);

    TFHEpp::Polynomial<P> phase{};
    for (std::size_t i = 0; i < P::n; i++)
        phase[i] = ct[1][i] - cross[i] + square[i];
    return phase;
}

}  // namespace

int main()
{
    constexpr std::uint32_t lhs_log_delta = 42;
    constexpr std::uint32_t rhs_log_delta = 40;
    constexpr std::uint32_t out_log_delta = rhs_log_delta;
    constexpr std::uint32_t lhs_log_q = 128;
    constexpr std::uint32_t rhs_log_q = 128;
    constexpr std::uint32_t out_log_q = lhs_log_q - lhs_log_delta;
    constexpr double tol = 1e-6;

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

    TFHEpp::Polynomial<P> lhs_poly{};
    TFHEpp::Polynomial<P> rhs_poly{};
    TFHEpp::ckksEncodePolynomial<P, lhs_log_q, lhs_log_delta>(lhs_poly, *lhs);
    TFHEpp::ckksEncodePolynomial<P, rhs_log_q, rhs_log_delta>(rhs_poly, *rhs);

    TFHEpp::TRLWE<P> lhs_ct{};
    TFHEpp::TRLWE<P> rhs_ct{};
    lhs_ct[P::k] = lhs_poly;
    rhs_ct[P::k] = rhs_poly;

    TFHEpp::TRLWE3<P> product{};
    TFHEpp::CKKSTensorProductRescale<P, lhs_log_q, rhs_log_q, lhs_log_delta>(
        product, lhs_ct, rhs_ct);
    TFHEpp::ckksDecodePolynomial<P, out_log_q, out_log_delta>(*got, product[1]);

    for (std::size_t i = 0; i < 8; i++)
        require_close((*got)[i], (*expected)[i], tol, "CKKS tensor rescale");
    for (std::size_t i = 8; i < P::n; i++)
        require_close((*got)[i], 0.0, tol, "CKKS tensor rescale zero tail");
    for (std::size_t i = 0; i < P::n; i++) {
        require_close(TFHEpp::ckksDecodeCoeff<P, out_log_q, out_log_delta>(
                          product[0][i]),
                      0.0, tol, "CKKS tensor cross zero");
        require_close(TFHEpp::ckksDecodeCoeff<P, out_log_q, out_log_delta>(
                          product[2][i]),
                      0.0, tol, "CKKS tensor square zero");
    }

    auto key = std::make_unique<TFHEpp::Key<P>>();
    fill_key(*key);

    constexpr std::uint32_t enc_log_q = 88;
    constexpr std::uint32_t enc_out_log_q = enc_log_q - lhs_log_delta;
    using LhsEnc = TFHEpp::CKKSCiphertext<P, enc_log_q, lhs_log_delta>;
    using RhsEnc = TFHEpp::CKKSCiphertext<P, enc_log_q, rhs_log_delta>;
    auto lhs_enc = std::make_unique<LhsEnc>();
    auto rhs_enc = std::make_unique<RhsEnc>();
    TFHEpp::ckksEncrypt<P, enc_log_q, lhs_log_delta>(*lhs_enc, *lhs, *key);
    TFHEpp::ckksEncrypt<P, enc_log_q, rhs_log_delta>(*rhs_enc, *rhs, *key);

    TFHEpp::TRLWE3<P> encrypted_product{};
    TFHEpp::CKKSTensorProductRescale<P, enc_log_q, enc_log_q, lhs_log_delta>(
        encrypted_product, lhs_enc->ct, rhs_enc->ct);
    const TFHEpp::Polynomial<P> encrypted_phase =
        trlwe3Phase(encrypted_product, *key);
    TFHEpp::ckksDecodePolynomial<P, enc_out_log_q, out_log_delta>(
        *got, encrypted_phase);

    for (std::size_t i = 0; i < 8; i++)
        require_close((*got)[i], (*expected)[i], 0.25,
                      "CKKS encrypted tensor phase");
    for (std::size_t i = 8; i < P::n; i++)
        require_close((*got)[i], 0.0, 0.25,
                      "CKKS encrypted tensor phase zero tail");

    std::cout << "Passed" << std::endl;
}
