#include <tfhe++.hpp>
#include <algorithm>
#include <chrono>
#include <iostream>
#include <string_view>
#include <random>

namespace {
struct SmallCKKS : TFHEpp::CKKSInteger4Param {
    static constexpr std::uint32_t nbit = 4, n = 1U << nbit;
};
struct SmallTFHE : TFHEpp::lvl2param {
    static constexpr std::uint32_t nbit = 3, n = 1U << nbit;
};

void invalid_input_check()
{
    const auto require_rejection = [](auto &&call) {
        bool rejected = false;
        try { call(); }
        catch (const std::invalid_argument &) { rejected = true; }
        if (!rejected) throw std::runtime_error("invalid secret was accepted");
    };
    TFHEpp::Key<SmallCKKS> ckks{};
    TFHEpp::Key<SmallTFHE> tfhe{};
    auto direct = std::make_unique<TFHEpp::CKKSToTFHEKeySwitchKey<SmallCKKS, SmallTFHE, 2, 24>>();
    auto combined = std::make_unique<TFHEpp::TFHECKKSInteger4Key<SmallCKKS, SmallTFHE>>();
    ckks[0] = 2;
    require_rejection([&] { TFHEpp::CKKSToTFHEKeyGen<SmallCKKS, SmallTFHE, 2, 24>(*direct, ckks, tfhe); });
    require_rejection([&] { TFHEpp::TFHECKKSInteger4KeyGen(*combined, ckks, tfhe); });
    ckks[0] = 0;
    tfhe[0] = 2;
    require_rejection([&] { TFHEpp::TFHECKKSInteger4KeyGen(*combined, ckks, tfhe); });
    TFHEpp::TLWE<SmallTFHE> ct{};
    bool rejected = false;
    try { TFHEpp::TFHEEncryptInteger4<SmallTFHE>(ct, 16, tfhe); }
    catch (const std::out_of_range &) { rejected = true; }
    if (!rejected) throw std::runtime_error("out-of-range message was accepted");
}

void exact_product_check()
{
    using P = SmallCKKS;
    TFHEpp::Polynomial<P> a{}, b{}, expected{}, got{};
    for (std::uint32_t i = 0; i < P::n; i++) {
        a[i] = TFHEpp::ckks_detail::uniformAtLevel<P, 192>();
        b[i] = TFHEpp::ckks_detail::uniformAtLevel<P, 192>();
    }
    for (std::uint32_t i = 0; i < P::n; i++)
        for (std::uint32_t j = 0; j < P::n; j++) {
            const auto term = a[i] * b[j];
            if (i+j < P::n) expected[i+j] += term;
            else expected[i+j-P::n] -= term;
        }
    TFHEpp::ckks_integer_detail::exactProduct<P, 192>(got, a, b);
    for (std::uint32_t i = 0; i < P::n; i++)
        if (got[i] != TFHEpp::ckks_detail::reduceToLevel<P, 192>(expected[i]))
            throw std::runtime_error("NTT product disagrees with schoolbook oracle");
}

template <class P, class L, unsigned Bits = 4>
int run(bool exhaustive, bool sampled = false)
{
    auto ckks_key = std::make_unique<TFHEpp::Key<P>>();
    TFHEpp::Key<L> tfhe_key{};
    TFHEpp::keyGen<P>(*ckks_key);
    TFHEpp::keyGen<L>(tfhe_key);
    auto key = std::make_unique<TFHEpp::TFHECKKSInteger4Key<P, L>>();
    std::cout << "N=" << P::n << " TFHE n=" << L::n
              << " direct key bytes=" << sizeof(*key->to_tfhe) << std::endl;
    const auto start = std::chrono::steady_clock::now();
    TFHEpp::TFHECKKSInteger4KeyGen<P, L>(*key, *ckks_key, tfhe_key);
    const auto ready = std::chrono::steady_clock::now();
    std::cout << "Key generation seconds=" << std::chrono::duration<double>(ready-start).count() << std::endl;
    std::vector<std::pair<unsigned, unsigned>> pairs;
    if (exhaustive) {
        for (unsigned x = 0; x < (1U << Bits); x++)
            for (unsigned y = 0; y < (1U << Bits); y++) pairs.emplace_back(x,y);
    }
    else {
        constexpr unsigned m = (1U << Bits)-1;
        pairs = {{0,m},{1,m},{m,m},{(m+1)/2,(m+1)/2},{m/2,m/2+2},{m-1,m},{0,0},{m,m-1}};
        if (sampled) {
            std::mt19937 rng(20260905);
            while (pairs.size() < 256) pairs.emplace_back(rng() & m, rng() & m);
        }
    }
    long double max_error = 0;
    std::size_t done = 0;
    for (const auto &[x,y] : pairs) {
        TFHEpp::TLWE<L> lhs{}, rhs{}, product{};
        TFHEpp::TFHEEncryptInteger<Bits, L>(lhs, x, tfhe_key);
        TFHEpp::TFHEEncryptInteger<Bits, L>(rhs, y, tfhe_key);
        TFHEpp::TFHECKKSIntegerMultiply<Bits, P, L>(product, lhs, rhs, *key);
        const std::uint64_t phase = TFHEpp::tlweSymPhase<L>(product, tfhe_key);
        const unsigned decoded = TFHEpp::TFHEDecryptIntegerProduct<Bits, L>(product, tfhe_key);
        const std::uint64_t difference = phase - (std::uint64_t{x*y} << (64-2*Bits));
        const long double error = std::ldexp(std::abs(static_cast<long double>(
            static_cast<std::int64_t>(difference))), -(64-2*static_cast<int>(Bits)));
        max_error = std::max(max_error, error);
        if (decoded != x*y) throw std::runtime_error("incorrect integer product");
        if (!exhaustive && !sampled) std::cout << x << "*" << y << "=" << decoded << " error=" << error << std::endl;
        if ((exhaustive || sampled) && P::n >= 128 && ++done % 16 == 0)
            std::cout << "Checked " << done << "/" << pairs.size() << "; max error=" << max_error << std::endl;
    }
    std::cout << "Products=" << pairs.size() << " maximum integer-unit error=" << max_error
              << " evaluation seconds=" << std::chrono::duration<double>(std::chrono::steady_clock::now()-ready).count()
              << "\nPassed" << std::endl;
    return 0;
}
}

int main(int argc, char **argv)
try {
    invalid_input_check();
    exact_product_check();
    if (argc == 1) return run<SmallCKKS, SmallTFHE>(true);
    if (argc == 2 && std::string_view(argv[1]) == "--8bit-small")
        return run<SmallCKKS, SmallTFHE, 8>(true);
    if (argc == 2 && std::string_view(argv[1]) == "--8bit")
        return run<TFHEpp::CKKSInteger4Param, TFHEpp::lvl2param, 8>(false, true);
    if (argc == 2 && std::string_view(argv[1]) == "--full")
        return run<TFHEpp::CKKSInteger4Param, TFHEpp::lvl2param>(false);
    if (argc == 2 && std::string_view(argv[1]) == "--full-exhaustive")
        return run<TFHEpp::CKKSInteger4Param, TFHEpp::lvl2param>(true);
    throw std::invalid_argument("usage: ckks_integer_multiplication [--full | --full-exhaustive | --8bit | --8bit-small]");
}
catch (const std::exception &e) {
    std::cerr << e.what() << std::endl;
    return 1;
}
