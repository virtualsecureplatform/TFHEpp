// bfv_gbfv_scaffold.cpp -- regression coverage for the first GBFV-style
// recryption primitives ported from the reference implementation.

#include <cstdint>
#include <iostream>
#include <random>
#include <tfhe++.hpp>

template <class P>
bool checkPoly(const std::array<uint64_t, P::n> &got,
               const std::array<uint64_t, P::n> &want, const char *label)
{
    int mismatches = 0;
    int first_bad = -1;
    for (int i = 0; i < static_cast<int>(P::n); i++) {
        if (got[i] != want[i]) {
            if (first_bad < 0) first_bad = i;
            mismatches++;
        }
    }
    if (mismatches == 0) return true;

    std::cerr << "  FAIL " << label << " mismatches=" << mismatches
              << " first_bad=" << first_bad
              << " expected=" << want[first_bad]
              << " got=" << got[first_bad] << std::endl;
    return false;
}

template <class P>
void transparentPolyCiphertext(TFHEpp::TRLWE<P> &ct,
                               const std::array<uint64_t, P::n> &poly)
{
    for (int c = 0; c <= static_cast<int>(P::k); c++)
        for (uint32_t i = 0; i < P::n; i++) ct[c][i] = 0;

    for (uint32_t i = 0; i < P::n; i++)
        ct[P::k][i] = TFHEpp::bfvEncodeCoeff<P>(poly[i]);
}

int main()
{
    using P = TFHEpp::lvl3simdparam;
    using BP = TFHEpp::bfvboot::PrimePower2Param<P>;

    constexpr uint64_t p = static_cast<uint64_t>(P::plain_modulus);
    constexpr uint64_t p2 = static_cast<uint64_t>(BP::plain_modulus);
    static_assert(p2 == p * p);

    std::cout << "BFV GBFV scaffold test" << std::endl;
    std::cout << "  n=" << P::n << " p=" << p << " p^2=" << p2
              << std::endl;

    std::default_random_engine engine(0x47424656);
    std::uniform_int_distribution<int> binary(0, 1);

    TFHEpp::Key<P> key{};
    for (auto &v : key)
        v = binary(engine) ? static_cast<typename P::T>(1)
                           : static_cast<typename P::T>(-1);

    TFHEpp::Key<BP> boot_key{};
    TFHEpp::bfvboot::ConvertSecretKey<P, BP>(boot_key, key);

    std::array<uint64_t, P::n> base_plain{};
    for (uint32_t i = 0; i < P::n; i++)
        base_plain[i] = (7 * static_cast<uint64_t>(i) + 3) % p;

    TFHEpp::TRLWE<P> base_ct{};
    transparentPolyCiphertext<P>(base_ct, base_plain);

    TFHEpp::TRLWE<BP> lifted_ct{};
    TFHEpp::bfvboot::gbfv::ScaleAndRoundPlainMod<P, BP>(lifted_ct, base_ct);

    std::array<uint64_t, BP::n> lifted_out{};
    TFHEpp::bfvboot::BfvPolyDecrypt<BP>(lifted_out, lifted_ct, boot_key);
    if (!checkPoly<BP>(lifted_out, base_plain, "transparent p-to-p^2 scale"))
        return 1;
    std::cout << "." << std::flush;

    std::array<uint64_t, BP::n> boot_plain{};
    for (uint32_t i = 0; i < BP::n; i++)
        boot_plain[i] = (11 * static_cast<uint64_t>(i) + 5) % p;

    TFHEpp::TRLWE<BP> boot_ct{};
    TFHEpp::bfvboot::BfvPolyEncrypt<BP>(boot_ct, boot_plain, boot_key);

    TFHEpp::TRLWE<P> projected_ct{};
    TFHEpp::bfvboot::gbfv::ScaleAndRoundPlainMod<BP, P>(projected_ct,
                                                        boot_ct);

    std::array<uint64_t, P::n> projected_out{};
    TFHEpp::bfvboot::BfvPolyDecrypt<P>(projected_out, projected_ct, key);
    if (!checkPoly<P>(projected_out, boot_plain, "p^2-to-p scale"))
        return 1;
    std::cout << "." << std::flush;

    std::array<uint64_t, P::n> trace_plain{};
    for (uint32_t i = 0; i < P::n; i++)
        trace_plain[i] = (13 * static_cast<uint64_t>(i) + 7) % 29;

    std::array<uint64_t, P::n> trace_want{};
    TFHEpp::bfvboot::gbfv::EvaluateTracePlainMod<P>(trace_want, trace_plain,
                                                    1, 1);

    TFHEpp::TRLWE<P> trace_ct{};
    TFHEpp::bfvboot::BfvPolyEncrypt<P>(trace_ct, trace_plain, key);

    TFHEpp::bfvboot::gbfv::TraceKey<P> trace_key;
    TFHEpp::bfvboot::gbfv::TraceKeyGen<P>(trace_key, key, 1, 1);

    TFHEpp::TRLWE<P> trace_res{};
    TFHEpp::bfvboot::gbfv::EvaluateTrace<P>(trace_res, trace_ct, trace_key,
                                            1, 1);

    std::array<uint64_t, P::n> trace_out{};
    TFHEpp::bfvboot::BfvPolyDecrypt<P>(trace_out, trace_res, key);
    if (!checkPoly<P>(trace_out, trace_want, "single trace stage"))
        return 1;
    std::cout << " PASS" << std::endl;

    return 0;
}
