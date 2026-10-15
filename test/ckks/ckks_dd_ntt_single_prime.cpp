#include "ckks_narrow_basis_benchmark.hpp"
#include <ckks/prepared_dd_ntt.hpp>
#include <iostream>
#include <random>

using namespace TFHEpp;
struct P : narrow_benchmark::BenchmarkParam<9,7> {
    static constexpr std::uint32_t nbit=4,n=16;
};

template <std::uint32_t Bg,std::uint32_t Bbar>
void check(std::mt19937_64 &random)
{
    using K=CKKSCompactDDRelinKey<P,433,Bg,Bbar>;
    auto key=std::make_unique<K>();
    for(auto &row:key->data)for(auto &secondary:row)for(auto &component:secondary)
        for(auto &coefficient:component)
            coefficient=K::pack(P::T::from_signed_i64(
                std::int64_t(random()%(std::uint64_t{1}<<Bbar))-(std::int64_t{1}<<(Bbar-1))));
    CKKSPreparedDDRelinKeyNTT<P,433,Bg,Bbar> prepared(*key);
    static_assert((Bg==20 && decltype(prepared)::prime_count==1) ||
                  (Bg==64 && decltype(prepared)::prime_count==2));
    auto input=std::make_unique<Polynomial<P>>();
    auto digits=std::make_unique<std::array<Polynomial<P>,K::primary_rows>>();
    auto result=std::make_unique<TRLWE<P>>();
    for(std::size_t trial=0;trial<32;++trial) {
        for(auto &value:*input) {
            for(auto &word:value.limb)word=random();
            value=ckks_detail::reduceToLevel<P,433>(value);
        }
        ckks_detail::activeBaseDecomposePolynomialRows<P,433,Bg,K::primary_rows>(*digits,*input);
        prepared.apply(*result,*input);
        for(std::size_t component=0;component<2;++component)
            for(std::size_t i=0;i<P::n;++i) {
                P::T expected{};
                for(std::size_t secondary=0;secondary<K::bbar_rows;++secondary) {
                    __int128 convolution=0;
                    for(std::size_t row=0;row<K::primary_rows;++row)
                        for(std::size_t j=0;j<P::n;++j) {
                            const auto digit=std::bit_cast<std::int64_t>((*digits)[row][j].limb[0]);
                            const auto index=(i+P::n-j)%P::n;
                            convolution+=(j<=i?1:-1)*static_cast<__int128>(digit)*key->data[row][secondary][component][index];
                        }
                    expected+=decltype(prepared)::signedTorus(convolution)<<((K::bbar_rows-secondary-1)*Bbar);
                }
                if((*result)[component][i]!=ckks_detail::reduceToLevel<P,433>(expected))
                    throw std::runtime_error("NTT DD differs from independent coefficient convolution");
            }
    }
    std::cout<<"DD primary="<<Bg<<" secondary="<<Bbar<<" primes="<<decltype(prepared)::prime_count<<" exact coefficients="<<32*2*P::n<<'\n';
}
int main()
{
    std::mt19937_64 random(20262119);
    check<20,20>(random);check<64,32>(random);
}
