// Deterministic noncryptographic rows; no encryption or key generation.
// The legacy function preserves the previous full-torus digit path verbatim
// for complete-operation timing and byte comparison outside each timer.
#include "../ckks_narrow_basis_benchmark.hpp"
#include <ckks/prepared_auto.hpp>
#include <iomanip>
#include <iostream>
using namespace TFHEpp;
template <class P, std::uint32_t LogQ, std::uint32_t SecondaryBits, std::uint32_t PrimaryBits>
inline void legacyKeySwitchRows(TRLWE<P> &result, const Polynomial<P> &poly,
    const CKKSPreparedAutoKeyFFT<P, LogQ, SecondaryBits, PrimaryBits> &key, std::size_t source)
{
    static_assert(use_multilimb_digit_fft_v<P>);
    using Cache = CKKSPreparedAutoKeyFFT<P, LogQ, SecondaryBits, PrimaryBits>;
    // Sum all primary rows before rounding the inverse FFT. Integer products
    // have ample mantissa headroom even after accumulation across those rows.
    static_assert(PrimaryBits + SecondaryBits + P::nbit + std::bit_width(Cache::rows) + 3 <
                  std::numeric_limits<double>::digits);
    if (source >= P::k || key.spectra.size() != Cache::spectrum_count)
        throw std::invalid_argument("invalid prepared auto key");
    auto centered = std::make_unique<Polynomial<P>>();
    ckks_detail::centeredPolynomialAtLevel<P, LogQ>(*centered, poly);
    auto digits = std::make_unique<std::array<Polynomial<P>, Cache::rows>>();
    ckks_detail::baseBbarDecomposePolynomialRows<ckks_detail::AutoPrimaryParam<P, PrimaryBits>, Cache::first>(*digits, *centered);
    std::vector<PolynomialInFD<P>> transformed(Cache::rows);
#pragma omp parallel for schedule(static) if(P::n >= 8192)
    for (std::size_t row = 0; row < Cache::rows; row++)
        TwistIFFTDigit<P>(transformed[row], (*digits)[row]);
    // FFT products fit signed64; expand only when recombining torus coefficients.
    std::vector<std::array<std::array<std::uint64_t, P::n>, P::k + 1>> products(Cache::secondary_rows);
#pragma omp parallel
    {
        auto accumulated = std::make_unique<PolynomialInFD<P>>();
#pragma omp for collapse(2) schedule(static)
        for (std::size_t secondary = 0; secondary < Cache::secondary_rows; secondary++)
            for (std::size_t c = 0; c <= P::k; c++) {
                accumulated->fill(0.0);
                for (std::size_t row = 0; row < Cache::rows; row++)
                    FMAInFD<P::n>(*accumulated, transformed[row], key.at(source, row, secondary, c));
                TwistFFTDigitProduct64<P>(products[secondary][c], *accumulated);
            }
    }
    if constexpr (P::n >= 8192) {
        // Each worker owns complete coefficients. Keep the secondary order,
        // but avoid streaming the full torus output once per secondary row.
#pragma omp parallel for collapse(2) schedule(static)
        for (std::size_t c = 0; c <= P::k; c++)
            for (std::size_t i = 0; i < P::n; i++) {
                typename P::T coefficient{};
                for (std::size_t secondary = 0; secondary < Cache::secondary_rows; secondary++) {
                    const auto shift = std::numeric_limits<typename P::T>::digits -
                        (Cache::secondary_first + secondary + 1) * SecondaryBits;
                    coefficient += static_cast<typename P::T>(static_cast<std::int64_t>(products[secondary][c][i])) << shift;
                }
                result[c][i] = ckks_detail::reduceToLevel<P, LogQ>(coefficient);
            }
    }
    else {
        for (auto &component : result) component.fill(typename P::T{});
        for (std::size_t secondary = 0; secondary < Cache::secondary_rows; secondary++) {
            const auto shift = std::numeric_limits<typename P::T>::digits -
                (Cache::secondary_first + secondary + 1) * SecondaryBits;
            for (std::size_t c = 0; c <= P::k; c++)
                for (std::size_t i = 0; i < P::n; i++)
                    result[c][i] += static_cast<typename P::T>(static_cast<std::int64_t>(products[secondary][c][i])) << shift;
        }
        ckks_detail::reduceTRLWEToLevel<P, LogQ>(result);
    }
}

template<std::uint32_t Q,std::uint32_t Bits>void benchmark()
{
    using P=narrow_benchmark::BenchmarkParam<14,7>;
    auto raw=std::make_unique<CKKSAutoKey<P,Q>>();
    std::mt19937_64 random(20261024+Q+Bits);
    for(auto&source:raw->data)for(auto&row:source)for(auto&component:row)
        for(auto&v:component){for(auto&w:v.limb)w=random();v=ckks_detail::reduceToLevel<P,Q>(v);}
    auto prepared=makeCKKSPreparedAutoKeyFFT<P,Q,16,Bits>(*raw);
    auto input=std::make_unique<Polynomial<P>>();
    for(auto&v:*input)for(auto&w:v.limb)w=random();
    auto reference=std::make_unique<TRLWE<P>>(),actual=std::make_unique<TRLWE<P>>();
    legacyKeySwitchRows<P,Q>(*reference,*input,*prepared,0);
    const auto run=[&](bool compact){
        const auto start=std::chrono::steady_clock::now();
        if(compact)CKKSPreparedKeySwitchRows<P,Q>(*actual,*input,*prepared,0);
        else legacyKeySwitchRows<P,Q>(*actual,*input,*prepared,0);
        const auto seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        if(*actual!=*reference)throw std::runtime_error("compact digit key switch differs from legacy");
        return seconds;
    };
    run(false);run(true);
    for(bool compact:{false,true,true,false})for(int trial=0;trial<5;++trial)
        std::cout<<"q="<<Q<<" bits="<<Bits<<" compact="<<compact<<" trial="<<trial<<" seconds="<<run(compact)<<" identical\n";
}
int main(){std::cout<<std::setprecision(17);benchmark<93,4>();benchmark<56,4>();benchmark<93,8>();benchmark<56,8>();}
