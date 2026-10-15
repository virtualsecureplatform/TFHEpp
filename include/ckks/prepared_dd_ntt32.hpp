#pragma once
#include "compact_dd.hpp"
#include <modular_ntt32.hpp>

#ifdef TFHEPP_HAS_EXTENDED_MULTILIMB_PARAMS
namespace TFHEpp {
template<class P,std::uint32_t Q,std::uint32_t Bg,std::uint32_t Bbar>
constexpr bool CKKSFitsDDNTT32()
{
    constexpr auto rows=CKKSCompactDDRelinKey<P,Q,Bg,Bbar>::primary_rows;
    if constexpr((modular_ntt32::primes[0].value-1)%(2*P::n) ||
                 (modular_ntt32::primes[1].value-1)%(2*P::n))return false;
    constexpr auto limit=std::uint64_t(modular_ntt32::primes[0].value)*modular_ntt32::primes[1].value/2;
    return (static_cast<unsigned __int128>(1)<<(Bg-1)) <=
        (limit-1)/(static_cast<unsigned __int128>(1)<<(Bbar-1))/P::n/rows;
}

// Exact small-prime representation of the existing encrypted DD matrix.
// Primary/secondary gadgets, decomposition and recomposition are unchanged.
template<class P,std::uint32_t Q,std::uint32_t Bg,std::uint32_t Bbar>
class CKKSPreparedDDRelinKeyNTT32 {
    static_assert(CKKSFitsDDNTT32<P,Q,Bg,Bbar>());
    using Key=CKKSCompactDDRelinKey<P,Q,Bg,Bbar>;
    static constexpr auto primary_rows=Key::primary_rows,secondary_rows=Key::bbar_rows;
    std::array<std::unique_ptr<modular_ntt32::NegacyclicNTTPlan>,2> plans_;
    std::vector<std::uint32_t> spectra_;
    static std::size_t index(std::size_t prime,std::size_t row,std::size_t secondary,
                             std::size_t component,std::size_t i)
    { return ((((prime*primary_rows+row)*secondary_rows+secondary)*2+component)*P::n+i); }
public:
    explicit CKKSPreparedDDRelinKeyNTT32(const Key&key)
        :spectra_(2*primary_rows*secondary_rows*2*P::n)
    {
        for(std::size_t p=0;p<2;++p)
            plans_[p]=std::make_unique<modular_ntt32::NegacyclicNTTPlan>(P::n,modular_ntt32::primes[p]);
#pragma omp parallel
        {
            std::vector<std::uint32_t> values(P::n);
#pragma omp for collapse(4) schedule(static)
            for(std::size_t p=0;p<2;++p)
                for(std::size_t row=0;row<primary_rows;++row)
                    for(std::size_t secondary=0;secondary<secondary_rows;++secondary)
                        for(std::size_t c=0;c<2;++c){
                            const auto &arithmetic=plans_[p]->arithmetic();const auto modulus=arithmetic.modulus();
                            for(std::size_t i=0;i<P::n;++i){
                                const std::int64_t v=key.data[row][secondary][c][i];
                                const auto magnitude=v<0?std::uint64_t{0}-std::uint64_t(v):std::uint64_t(v);
                                const auto reduced=magnitude%modulus;
                                values[i]=arithmetic.encode(v<0&&reduced?modulus-reduced:reduced);
                            }
                            plans_[p]->forward(values);
                            std::copy(values.begin(),values.end(),spectra_.begin()+index(p,row,secondary,c,0));
                        }
        }
    }
    std::size_t byteSize()const{return spectra_.size()*sizeof(spectra_[0]);}
    void apply(TRLWE<P>&result,const Polynomial<P>&input)const
    {
        if(spectra_.size()!=2*primary_rows*secondary_rows*2*P::n || !plans_[0] || !plans_[1])
            throw std::invalid_argument("invalid small-prime DD cache");
        auto digits=std::make_unique<std::array<Polynomial<P>,primary_rows>>();
        ckks_detail::activeBaseDecomposePolynomialRows<P,Q,Bg,primary_rows>(*digits,input);
        std::vector<std::uint32_t> transformed(2*primary_rows*P::n);
#pragma omp parallel for collapse(2) schedule(static)
        for(std::size_t p=0;p<2;++p)
            for(std::size_t row=0;row<primary_rows;++row){
                const auto modulus=plans_[p]->arithmetic().modulus();
                auto values=std::span(transformed).subspan((p*primary_rows+row)*P::n,P::n);
                for(std::size_t i=0;i<P::n;++i){
                    const auto v=std::bit_cast<std::int64_t>((*digits)[row][i].limb[0]);
                    const auto magnitude=v<0?std::uint64_t{0}-std::uint64_t(v):std::uint64_t(v);
                    const auto reduced=magnitude%modulus;
                    values[i]=v<0&&reduced?modulus-reduced:reduced;
                }
                plans_[p]->forward(values);
            }
        const auto position=[](std::size_t p,std::size_t s,std::size_t c,std::size_t i){return ((p*secondary_rows+s)*2+c)*P::n+i;};
        std::vector<std::uint32_t> products(2*secondary_rows*2*P::n);
#pragma omp parallel for collapse(3) schedule(static)
        for(std::size_t p=0;p<2;++p)
            for(std::size_t s=0;s<secondary_rows;++s)
                for(std::size_t c=0;c<2;++c){
                    const auto &arithmetic=plans_[p]->arithmetic();
                    auto values=std::span(products).subspan(position(p,s,c,0),P::n);
                    for(std::size_t begin=0;begin<P::n;begin+=256)
                        for(std::size_t row=0;row<primary_rows;++row){
                            const auto end=std::min(begin+256,std::size_t(P::n));std::size_t i=begin;
                            const auto *key=spectra_.data()+index(p,row,s,c,0);
                            const auto *data=transformed.data()+(p*primary_rows+row)*P::n;
#if defined(__AVX2__)
                            for(;i+8<=end;i+=8){
                                const auto sum=arithmetic.accumulateVector(
                                    _mm256_loadu_si256(reinterpret_cast<const __m256i*>(values.data()+i)),
                                    _mm256_loadu_si256(reinterpret_cast<const __m256i*>(data+i)),
                                    _mm256_loadu_si256(reinterpret_cast<const __m256i*>(key+i)));
                                _mm256_storeu_si256(reinterpret_cast<__m256i*>(values.data()+i),sum);
                            }
#endif
                            for(;i<end;++i){
                                const auto sum=values[i]+arithmetic.multiplyLazy(data[i],key[i]);
                                values[i]=sum>=2*arithmetic.modulus()?sum-2*arithmetic.modulus():sum;
                            }
                        }
                    for(auto &v:values)v=arithmetic.canonical(v);
                    plans_[p]->inverse(values);
                }
        constexpr auto p0=modular_ntt32::primes[0].value,p1=modular_ntt32::primes[1].value;
        constexpr auto modulus=p0*p1;
        const auto inverse=modular_ntt::invert(p0%p1,p1);
#pragma omp parallel for collapse(2) schedule(static)
        for(std::size_t c=0;c<2;++c)
            for(std::size_t i=0;i<P::n;++i){
                WideSignedLimbAccumulator<P::T::limbs> sum;
                for(std::size_t s=0;s<secondary_rows;++s){
                    const auto a=products[position(0,s,c,i)],b=products[position(1,s,c,i)];
                    const auto reduced=a>=p1?a-p1:a;
                    const auto digit=(b+p1-reduced)*inverse%p1;
                    const auto value=a+p0*digit;
                    const std::int64_t exact=value>modulus/2?-std::int64_t(modulus-value):std::int64_t(value);
                    sum.add_shifted_i64(exact,(secondary_rows-s-1)*Bbar);
                }
                typename P::T value;value.limb=sum.limb;
                result[c][i]=ckks_detail::reduceToLevel<P,Q>(value);
            }
    }
};
} // namespace TFHEpp
#endif
