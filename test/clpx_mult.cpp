#include <cstdint>
#include <iostream>
#include <tfhe++.hpp>

namespace {

template <class P, int validbit>
typename P::T decode_digits(const std::array<int, P::n> &digits)
{
    return TFHEpp::decodeHatEncoderInt2T<P, validbit, P::n>(digits);
}

template <class P, int validbit>
bool check_digits(const std::array<int, P::n> &digits,
                  const typename P::T expected, const char *label)
{
    const auto decoded = decode_digits<P, validbit>(digits);
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

template <class iksP, class bkP, class sskP, int num_multi, int shift, int w,
          int validbit>
bool switch_tlwes_to_clpx(typename bkP::targetP::T value,
                          TFHEpp::TRLWE<typename bkP::targetP> &switched,
                          const TFHEpp::SecretKey &sk, TFHEpp::EvalKey &ek,
                          const TFHEpp::AnnihilateKey<typename bkP::targetP> &ahk)
{
    std::vector<uint8_t> digits(validbit);
    for (int i = 0; i < validbit; i++)
        digits[i] = static_cast<uint8_t>((value >> i) & 1);
    std::vector<TFHEpp::TLWE<typename iksP::domainP>> tlwes;
    TFHEpp::bootsSymEncrypt<typename iksP::domainP>(tlwes, digits, sk);
    if (tlwes.size() != digits.size()) {
        std::cerr << "TLWE input size mismatch" << std::endl;
        return false;
    }

    TFHEpp::TLWES2CLPXIKS<iksP, bkP, sskP, num_multi, shift, w>(
        switched, tlwes, ahk, ek);
    return true;
}

template <int validbit>
bool check_switched_clpx_mult()
{
    using iksP = TFHEpp::lvl1hparam;
    using bkP = TFHEpp::SS2CLPXlvlh2param;
    using sskP = TFHEpp::SS2CLPXlvl22param;
    using P = TFHEpp::SS2CLPXlvl2param;
    constexpr int lutnum = 4;
    constexpr int shiftnum = 5;
    constexpr int shift = shiftnum - 1;
    constexpr int w = 20;

    TFHEpp::SecretKey sk;
    TFHEpp::EvalKey ek;
    ek.emplaceiksk<iksP>(sk);
    ek.emplacebkfft<TFHEpp::lvlh2param>(sk);

    const auto key = sk.key.get<TFHEpp::lvl2param>();
    TFHEpp::AnnihilateKey<P> ahk;
    TFHEpp::annihilatekeygen<P>(ahk, key);

    TFHEpp::TRLWE<P> lhs_ct;
    TFHEpp::TRLWE<P> rhs_ct;
    if (!switch_tlwes_to_clpx<iksP, bkP, sskP, lutnum, shift, w, validbit>(
            2, lhs_ct, sk, ek, ahk))
        return false;
    if (!switch_tlwes_to_clpx<iksP, bkP, sskP, lutnum, shift, w, validbit>(
            3, rhs_ct, sk, ek, ahk))
        return false;

    const auto lhs_digits = TFHEpp::clpxSymIntDecrypt<P>(lhs_ct, key);
    if (!check_digits<P, validbit>(lhs_digits, 2, "TLWES2CLPXIKS lhs"))
        return false;
    const auto rhs_digits = TFHEpp::clpxSymIntDecrypt<P>(rhs_ct, key);
    if (!check_digits<P, validbit>(rhs_digits, 3, "TLWES2CLPXIKS rhs"))
        return false;

    const auto relinkeyfft = TFHEpp::relinKeyFFTgen<P>(key);
    TFHEpp::TRLWE<P> product_ct;
    TFHEpp::CLPXMult<P>(product_ct, lhs_ct, rhs_ct, relinkeyfft);

    const auto product_digits = TFHEpp::clpxSymIntDecrypt<P>(product_ct, key);
    return check_digits<P, validbit * 2>(product_digits, 6,
                                         "TLWES2CLPXIKS + CLPXMult");
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
    if (!check_switched_clpx_mult<validbit>()) return 1;

    std::cout << "Passed" << std::endl;
}
