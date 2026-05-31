#include <complex>
#include <cmath>
#include <iostream>
#include <memory>
#include <tfhe++.hpp>

namespace {

using P = TFHEpp::ckkslvl3param;

void require_close(std::complex<double> got, std::complex<double> want,
                   double tol, const char *label)
{
    if (std::abs(got - want) > tol) {
        std::cerr << label << " got=(" << got.real() << ',' << got.imag()
                  << ") want=(" << want.real() << ',' << want.imag()
                  << ") tol=" << tol << std::endl;
        std::exit(1);
    }
}

}  // namespace

int main()
{
    constexpr std::uint32_t log_delta = 40;
    constexpr std::uint32_t log_q = 96;
    constexpr double tol = 1e-5;

    auto slots = std::make_unique<std::array<std::complex<double>, P::n / 2>>();
    auto decoded = std::make_unique<std::array<std::complex<double>, P::n / 2>>();
    slots->fill({0.0, 0.0});
    (*slots)[0] = {0.25, -0.125};
    (*slots)[1] = {-0.5, 0.375};
    (*slots)[7] = {0.0625, 0.03125};
    (*slots)[P::n / 2 - 1] = {-0.125, -0.25};

    TFHEpp::Polynomial<P> poly;
    TFHEpp::ckksSlotEncode<P, log_q, log_delta>(poly, *slots);
    TFHEpp::ckksSlotDecode<P, log_q, log_delta>(*decoded, poly);

    for (std::size_t i = 0; i < P::n / 2; i++)
        require_close((*decoded)[i], (*slots)[i], tol, "slot encode/decode");

    std::cout << "Passed" << std::endl;
}
