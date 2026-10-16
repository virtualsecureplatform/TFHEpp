#include "ckks_narrow_basis_benchmark.hpp"
#include <ckks/tlwe_coefficient_packing.hpp>

// Exact tiny-ring oracle for the packing arithmetic, including negacyclic
// wrap, both ciphertext components, all source digits, and inactive slots.
namespace {struct TinyL:TFHEpp::lvl1param {static constexpr std::uint32_t n=4;};}
int main()try {
    using namespace TFHEpp;
    using P=narrow_benchmark::BenchmarkParam<3,7>;
    using L=TinyL;
    using K=TLWECoefficientPackingKey48<P,L>;
    K key;key.data.resize(K::rows*2*P::n);
    std::uint64_t state=12345;
    const auto next=[&](){state^=state<<13;state^=state>>7;state^=state<<17;return state;};
    for(auto &v:key.data)v=next()&K::mask;
    PreparedTLWECoefficientPacking48<P,L> prepared(key);
    for(std::size_t count=1;count<=4;++count){
        std::vector<TLWE<L>> in(count);
        for(auto &ct:in)for(auto &v:ct)v=next();
        CKKSCiphertext<P,48,44> got;prepared.apply(got,in);
        std::array<std::array<std::uint64_t,P::n>,2> expected{};
        for(std::size_t row=0;row<K::rows;++row)for(std::size_t j=0;j<count;++j){
            const auto digit=(in[j][row/4]>>(8*(row%4)))&255;
            for(std::size_t c=0;c<2;++c)for(std::size_t i=0;i<P::n;++i){
                const auto value=key.data[K::index(row,c,i)]*digit;
                if(i+j<P::n)expected[c][i+j]-=value;else expected[c][i+j-P::n]+=value;
            }
        }
        for(std::size_t j=0;j<count;++j)expected[1][j]+=std::uint64_t(in[j][L::n])<<16;
        for(std::size_t c=0;c<2;++c)for(std::size_t i=0;i<P::n;++i)
            if(got.ct[c][i]!=P::T{expected[c][i]&K::mask})throw std::runtime_error("exact packing oracle mismatch");
    }
    std::cout<<"Exact coefficient packing oracle passed"<<std::endl;
}catch(const std::exception &e){std::cerr<<e.what()<<'\n';return 1;}
