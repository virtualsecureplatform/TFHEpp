#include <cstdint>
#include <iostream>
#include <tfhe++.hpp>

namespace {

template <class P, int validbit>
bool check_digits(const std::array<int, P::n> &digits,
                  const typename P::T expected, const char *label)
{
    const auto decoded =
        TFHEpp::decodeHatEncoderInt2T<P, validbit, P::n>(digits);

    if (decoded != expected) {
        std::cerr << label << " decoded as "
                  << static_cast<std::uint64_t>(decoded) << ", expected "
                  << static_cast<std::uint64_t>(expected) << "; first digits:";
        for (int i = 0; i < validbit; i++) std::cerr << ' ' << digits[i];
        std::cerr << std::endl;
        return false;
    }
    return true;
}

template <class P, int validbit>
bool check_clpx_roundtrip(const typename P::T value, const TFHEpp::Key<P> &key)
{
    const auto poly = TFHEpp::EncodeHatEncoderP<P>(value);
    const auto ct = TFHEpp::clpxSymIntEncrypt<P>(poly, key);
    const auto digits = TFHEpp::clpxSymIntDecrypt<P>(ct, key);
    return check_digits<P, validbit>(digits, value, "CLPX roundtrip");
}

template <class P, int validbit>
bool check_clpx_mult(const typename P::T lhs, const typename P::T rhs,
                     const TFHEpp::Key<P> &key,
                     const TFHEpp::relinKeyFFT<P> &relinkeyfft)
{
    const auto lhs_poly = TFHEpp::EncodeHatEncoderP<P>(lhs);
    const auto rhs_poly = TFHEpp::EncodeHatEncoderP<P>(rhs);
    const auto lhs_ct = TFHEpp::clpxSymIntEncrypt<P>(lhs_poly, key);
    const auto rhs_ct = TFHEpp::clpxSymIntEncrypt<P>(rhs_poly, key);

    TFHEpp::TRLWE<P> product_ct;
    TFHEpp::CLPXMult<P>(product_ct, lhs_ct, rhs_ct, relinkeyfft);

    const auto product_digits = TFHEpp::clpxSymIntDecrypt<P>(product_ct, key);
    return check_digits<P, validbit>(product_digits, lhs * rhs, "CLPXMult");
}

}  // namespace

int main()
{
    using P = TFHEpp::lvl2param;
    constexpr int validbit = 8;

    TFHEpp::SecretKey sk;
    const auto &key = sk.key.get<P>();
    const auto relinkeyfft = TFHEpp::relinKeyFFTgen<P>(key);

    if (!check_clpx_roundtrip<P, validbit>(2, key)) return 1;
    if (!check_clpx_roundtrip<P, validbit>(3, key)) return 1;

    if (!check_clpx_mult<P, validbit>(2, 3, key, relinkeyfft)) return 1;
    if (!check_clpx_mult<P, validbit>(5, 6, key, relinkeyfft)) return 1;
    if (!check_clpx_mult<P, validbit>(7, 7, key, relinkeyfft)) return 1;

    std::cout << "Passed" << std::endl;
}
