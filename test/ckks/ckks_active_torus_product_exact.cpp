#include <ckks/ckks.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <chrono>
#include <fstream>
#include <iostream>
#include <random>

namespace {
using namespace TFHEpp;
using boost::multiprecision::cpp_int;
template <std::uint32_t Nbit, std::size_t Limbs>
struct Param : lvl5param {
    static constexpr std::uint32_t nbit = Nbit, n = 1U << Nbit, B̅gbit = 20;
    using T = MultiLimbUInt<Limbs>;
};

template <class T> cpp_int integer(const T &value)
{
    cpp_int result = 0;
    for (std::size_t bit = 0; bit < T::digits; bit += 64)
        result += cpp_int(static_cast<std::uint64_t>(value >> bit)) << bit;
    return result;
}

// Frozen pre-change 16-bit algorithm. The separate multiprecision oracle
// below checks both algorithms rather than accepting their agreement alone.
template <class P, std::uint32_t Q>
void oldProduct(Polynomial<P> &out, const Polynomial<P> &torus, const Polynomial<P> &small)
{
    static const modular_ntt::NegacyclicNTTPlan plan(P::n, {4611686018427322369ULL, 7});
    const auto modulus = plan.modulus();
    out.fill({});
    std::vector<std::uint64_t> lhs(P::n), rhs(P::n);
    for (std::size_t shift = 0; shift < Q; shift += 16) {
        for (std::size_t i = 0; i < P::n; ++i) {
            lhs[i] = static_cast<std::uint64_t>(ckks_detail::reduceToLevel<P,Q>(torus[i]) >> shift) & 65535;
            const auto value = multilimb_to_signed_i64(small[i]);
            rhs[i] = value < 0 ? modulus - (-value) : value;
        }
        plan.forward(lhs); plan.forward(rhs);
        for (std::size_t i = 0; i < P::n; ++i) lhs[i] = modular_ntt::multiply(lhs[i], rhs[i], modulus);
        plan.inverse(lhs);
        for (std::size_t i = 0; i < P::n; ++i) {
            const auto value = lhs[i] <= modulus / 2 ? static_cast<std::int64_t>(lhs[i])
                                                   : -static_cast<std::int64_t>(modulus - lhs[i]);
            out[i] += P::T::from_signed_i64(value) << shift;
        }
    }
    for (auto &value : out) value = ckks_detail::reduceToLevel<P,Q>(value);
}

template <class P, std::uint32_t Q>
void check()
{
    std::mt19937_64 random(20262119 + P::n);
    auto torus = std::make_unique<Polynomial<P>>(), small = std::make_unique<Polynomial<P>>();
    auto actual = std::make_unique<Polynomial<P>>(), baseline = std::make_unique<Polynomial<P>>();
    const cpp_int mask = (cpp_int(1) << Q) - 1;
    std::size_t checked = 0;
    for (const int kind : {0, 1, 2}) {
        for (std::size_t i = 0; i < P::n; ++i) {
            (*torus)[i] = kind == 0 ? typename P::T{} : typename P::T{} - typename P::T{1};
            if (kind == 0) for (std::size_t bit = 0; bit < P::T::digits; bit += 64)
                (*torus)[i] |= typename P::T{random()} << bit;
            const auto value = kind == 2 ? (i < P::n / 4 ? 524287 : -524288)
                                        : static_cast<std::int64_t>(random() % 3) - 1;
            (*small)[i] = P::T::from_signed_i64(value);
        }
        if (kind == 0 && P::n > 64) {
            small->fill({});
            for (const std::size_t i : std::array<std::size_t,4>{0,1,P::n / 2,P::n - 1})
                (*small)[i] = P::T::from_signed_i64(i % 2 ? -524288 : 524287);
        }
        ckks_detail::exactActiveTorusMulBySmall<P,Q>(*actual, *torus, *small);
        oldProduct<P,Q>(*baseline, *torus, *small);
        if (*actual != *baseline) throw std::runtime_error("active product changed old ciphertext arithmetic");
        std::int64_t total = 0, prefix = 0;
        std::vector<std::size_t> nonzero;
        for (std::size_t i = 0; i < P::n; ++i) {
            total += multilimb_to_signed_i64((*small)[i]);
            if ((*small)[i] != typename P::T{}) nonzero.push_back(i);
        }
        for (std::size_t i = 0; i < P::n; ++i) {
            cpp_int expected = 0;
            if (kind != 0) {
                prefix += multilimb_to_signed_i64((*small)[i]);
                expected = mask * (2 * prefix - total);
            }
            else for (const auto j : nonzero) {
                const auto value = multilimb_to_signed_i64((*small)[j]);
                if (value == 0) continue;
                const auto index = (i + P::n - j) % P::n;
                expected += (integer((*torus)[index]) & mask) * value * (i < j ? -1 : 1);
            }
            if (integer((*actual)[i]) != (expected & mask))
                throw std::runtime_error("active product disagrees with independent integer convolution");
            ++checked;
        }
    }
    small->fill({});
    ckks_detail::exactActiveTorusMulBySmall<P,Q>(*actual, *torus, *small);
    for (const auto &value : *actual) if (value != typename P::T{}) throw std::runtime_error("zero multiplier mismatch");
    for (const auto invalid : {typename P::T{1} << 100,
                                P::T::from_signed_i64(std::numeric_limits<std::int64_t>::min())}) {
        small->fill({}); (*small)[0] = invalid;
        bool rejected = false;
        try { ckks_detail::exactActiveTorusMulBySmall<P,Q>(*actual, *torus, *small); }
        catch (const std::invalid_argument &) { rejected = true; }
        if (!rejected) throw std::runtime_error("out-of-bound active multiplier accepted");
    }
    std::cout << "{\"n\":" << P::n << ",\"q\":" << Q << ",\"independent_coefficients\":" << checked
              << ",\"baseline16_exact\":true,\"zero_multiplier_exact\":true}" << std::endl;
}

void retained(const char *state, const char *group)
{
    using P = Param<14,7>;
    std::ifstream file(state, std::ios::binary);
    cereal::PortableBinaryInputArchive archive(file);
    std::array<std::uint64_t,10> metadata{};
    std::string function; bool finished = false;
    archive(metadata,function,finished);
    const std::array<std::uint64_t,10> expected{0x4155583850523031ULL,16384,85,44,385,48,37,108,8,448};
    if (metadata != expected || function != "sigmoid" || !finished) throw std::runtime_error("bad retained active-product profile");
    auto secret = std::make_unique<Polynomial<P>>(); archive(*secret); file.close();
    std::ifstream source(group,std::ios::binary);
    cereal::PortableBinaryInputArchive group_archive(source);
    std::array<std::uint64_t,9> header{};group_archive(header);
    const std::array<std::uint64_t,9> expected_header{0x4e4152524f573031ULL,16384,385,48,37,529,448,0,19};
    if (header != expected_header) throw std::runtime_error("bad retained active-product selector header");
    auto selector = std::make_unique<CKKSCiphertext<P,385,48>>();
    source.seekg(73+466*static_cast<std::streamoff>(sizeof(*selector)));
    group_archive(*selector);source.close();
    auto actual = std::make_unique<Polynomial<P>>(), baseline = std::make_unique<Polynomial<P>>();
    ckks_detail::exactActiveTorusMulBySmall<P,385>(*actual,selector->ct[0],*secret);
    oldProduct<P,385>(*baseline,selector->ct[0],*secret);
    if (*actual != *baseline) throw std::runtime_error("retained mask product changed");
    std::vector<Polynomial<P>> concurrent(16);
#pragma omp parallel for schedule(static)
    for (std::size_t i = 0; i < concurrent.size(); ++i)
        ckks_detail::exactActiveTorusMulBySmall<P,385>(concurrent[i],selector->ct[0],*secret);
    for (const auto &value : concurrent)
        if (value != *baseline) throw std::runtime_error("concurrent active product differs");
    std::cout << "{\"concurrent_products\":16,\"exact_coefficients\":262144,\"passed\":true}" << std::endl;
    for (const bool old : {false,true,true,false}) {
        const auto start = std::chrono::steady_clock::now();
        for (int repeat = 0; repeat < 8; ++repeat) {
            if (old) oldProduct<P,385>(*actual,selector->ct[0],*secret);
            else ckks_detail::exactActiveTorusMulBySmall<P,385>(*actual,selector->ct[0],*secret);
        }
        if (*actual != *baseline) throw std::runtime_error("retained timing product mismatch");
        std::cout << "{\"old16\":" << (old ? "true" : "false") << ",\"repeats\":8,\"seconds\":"
                  << std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()
                  << ",\"retained_mask_coefficients_exact\":16384}" << std::endl;
    }
}
}

int main(int argc,char **argv)
{
    try {
        if (argc == 4 && std::string(argv[1]) == "--retained-mask") retained(argv[2],argv[3]);
        else if (argc == 1) { check<Param<6,9>,518>(); check<Param<14,7>,128>(); check<Param<14,7>,385>(); check<Param<14,7>,448>(); check<Param<15,9>,518>(); }
        else return 2;
    }
    catch (const std::exception &error) { std::cerr << error.what() << '\n'; return 2; }
}
