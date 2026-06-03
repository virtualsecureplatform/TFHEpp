#include <cmath>
#include <iostream>
#include <memory>
#include <tfhe++.hpp>

namespace {

using P = TFHEpp::ckkslvl3param;

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
    constexpr std::uint32_t log_delta = 40;
    constexpr std::uint32_t log_q = 80;
    constexpr double tol = std::ldexp(1.0, -static_cast<int>(log_delta) + 2);

    for (const double value :
         {0.0, 1.0, -1.0, 0.125, -0.375, 123.5, -17.25}) {
        const auto encoded =
            TFHEpp::ckksEncodeCoeff<P, log_q, log_delta>(value);
        const double decoded =
            TFHEpp::ckksDecodeCoeff<P, log_q, log_delta>(encoded);
        require_close(decoded, value, tol, "coefficient encode/decode");
    }

    auto input = std::make_unique<std::array<double, P::n>>();
    auto output = std::make_unique<std::array<double, P::n>>();
    input->fill(0.0);
    (*input)[0] = 0.25;
    (*input)[1] = -0.5;
    (*input)[7] = 1.75;
    (*input)[P::n - 1] = -0.125;

    TFHEpp::Polynomial<P> poly;
    TFHEpp::ckksEncodePolynomial<P, log_q, log_delta>(poly, *input);
    TFHEpp::ckksDecodePolynomial<P, log_q, log_delta>(*output, poly);

    for (std::size_t i = 0; i < P::n; i++)
        require_close((*output)[i], (*input)[i], tol, "polynomial encode/decode");

    std::cout << "Passed" << std::endl;
}
