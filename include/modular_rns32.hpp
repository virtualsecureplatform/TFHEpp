#pragma once
#include "modular_ntt32.hpp"
#include "bfv/bfv-multilimb.hpp"

namespace TFHEpp::modular_rns32 {
inline constexpr std::size_t prime_count=17;
inline constexpr std::array<modular_ntt::PrimeModulus,prime_count> primes{{
    {1073479681ULL,11},
    {1071513601ULL,23},
    {1070727169ULL,11},
    {1068236801ULL,3},
    {1065484289ULL,3},
    {1064697857ULL,3},
    {1062862849ULL,7},
    {1062469633ULL,5},
    {1060765697ULL,3},
    {1056440321ULL,3},
    {1056178177ULL,5},
    {1055260673ULL,3},
    {1054212097ULL,5},
    {1053818881ULL,7},
    {1052508161ULL,3},
    {1051721729ULL,6},
    {1049100289ULL,7},
}};
static_assert(primes.front().value < 2*primes.back().value);

// A 32-bit Shoup constant accepts the complete uint32 input range. Its lazy
// result is below 2p, including when input is larger than the modulus.
struct Multiplier {
    std::uint32_t value=0,quotient=0;
    Multiplier()=default;
    Multiplier(std::uint32_t x,std::uint32_t p):value(x) {
        if(p<3 || p>=(std::uint32_t{1}<<30) || x>=p)
            throw std::invalid_argument("invalid direct RNS multiplier");
        quotient=(std::uint64_t(x)<<32)/p;
    }
    std::uint32_t lazy(std::uint32_t x,std::uint32_t p) const {
        const auto q=(std::uint64_t(x)*quotient)>>32;
        return std::uint32_t(std::uint64_t(x)*value-q*p);
    }
    std::uint32_t apply(std::uint32_t x,std::uint32_t p) const {
        const auto r=lazy(x,p);return r>=p?r-p:r;
    }
};

template<std::size_t Limbs>
constexpr MultiLimbUInt<Limbs> multiplySmall(MultiLimbUInt<Limbs> x,std::uint32_t scalar) {
    unsigned __int128 carry=0;
    for(auto &word:x.limb){const auto v=static_cast<unsigned __int128>(word)*scalar+carry;word=std::uint64_t(v);carry=v>>64;}
    return x;
}
inline constexpr auto product=[] {
    MultiLimbUInt<8> x{1};for(auto p:primes)x=multiplySmall(x,p.value);return x;
}();

template<std::size_t Limbs,std::size_t Q>
class CanonicalEncoder {
    static_assert(Q>0 && Q<=64*Limbs);
    static constexpr std::size_t words=(Q+31)/32;
    std::array<std::array<Multiplier,words>,prime_count> radix_;
public:
    CanonicalEncoder(){
        for(std::size_t i=0;i<prime_count;++i){
            const auto p=primes[i].value;std::uint64_t radix=1;
            for(auto &v:radix_[i]){v=Multiplier(radix,p);radix=(radix<<32)%p;}
        }
    }
    // Interpret only the active Q bits, matching the torus modulus exactly.
    std::uint32_t residue(const MultiLimbUInt<Limbs>& x,std::size_t prime)const {
        const auto p=std::uint32_t(primes[prime].value);std::uint64_t sum=0;
        for(std::size_t w=0;w<words;++w){
            auto digit=std::uint32_t(x.limb[w/2]>>(32*(w%2)));
            if constexpr(Q%32)if(w==words-1)digit&=(std::uint32_t{1}<<(Q%32))-1;
            sum+=radix_[prime][w].lazy(digit,p);
        }
        return sum%p;
    }
};

// Garner mixed radix with exact centering. The digits of floor(M/2) are
// floor(p_i/2); compare from the highest radix to determine the sign without
// floating point or a full-width reconstruction. Accumulate only the target
// torus limbs, then subtract M modulo that same torus when negative.
template<std::size_t Limbs>
class CenteredCRT {
    using T=MultiLimbUInt<Limbs>;
    std::array<std::array<Multiplier,prime_count>,prime_count> inverse_{};
    std::array<std::array<std::uint32_t,prime_count>,prime_count> inverse_mont_{};
    std::vector<modular_ntt32::Montgomery> arithmetic_;
    std::array<T,prime_count> prefixes_;
    T modulus_;
public:
    CenteredCRT(){
        T product{1};arithmetic_.reserve(prime_count);
        for(std::size_t i=0;i<prime_count;++i){
            prefixes_[i]=product;product=multiplySmall(product,primes[i].value);
            arithmetic_.emplace_back(primes[i].value);
            for(std::size_t j=0;j<i;++j){
                const auto v=modular_ntt::invert(primes[j].value%primes[i].value,primes[i].value);
                inverse_[i][j]=Multiplier(v,primes[i].value);
                inverse_mont_[i][j]=arithmetic_[i].encode(v);
            }
        }
        modulus_=product;
    }
    // Canonical residue rows are overwritten with their mixed-radix digits.
    void reconstruct(std::array<std::span<std::uint32_t>,prime_count> residues,std::span<T> output)const {
        for(auto row:residues)if(row.size()!=output.size())throw std::invalid_argument("direct CRT lengths differ");
        for(std::size_t i=1;i<prime_count;++i){
            const auto p=std::uint32_t(primes[i].value);
            for(std::size_t j=0;j<i;++j){
                std::size_t k=0;
#if defined(__AVX2__)
                const auto prime=_mm256_set1_epi32(p);
                const auto fixed=_mm256_set1_epi32(inverse_mont_[i][j]);
                for(;k+8<=output.size();k+=8){
                    const auto previous=_mm256_loadu_si256(reinterpret_cast<const __m256i*>(residues[j].data()+k));
                    // Every previous modulus is less than twice this modulus.
                    const auto reduced=arithmetic_[i].canonicalVector(previous);
                    const auto current=_mm256_loadu_si256(reinterpret_cast<const __m256i*>(residues[i].data()+k));
                    const auto input=_mm256_sub_epi32(_mm256_add_epi32(current,prime),reduced);
                    const auto value=arithmetic_[i].canonicalVector(arithmetic_[i].multiplyLazyVector(input,fixed));
                    _mm256_storeu_si256(reinterpret_cast<__m256i*>(residues[i].data()+k),value);
                }
#endif
                for(;k<output.size();++k){
                    auto previous=residues[j][k];if(previous>=p)previous-=p;
                    residues[i][k]=inverse_[i][j].apply(residues[i][k]+p-previous,p);
                }
            }
        }
        for(std::size_t k=0;k<output.size();++k){
            T sum{};
            for(std::size_t i=0;i<prime_count;++i){
                unsigned __int128 carry=0;
                for(std::size_t limb=0;limb<Limbs;++limb){
                    const auto value=static_cast<unsigned __int128>(prefixes_[i].limb[limb])*residues[i][k]+sum.limb[limb]+carry;
                    sum.limb[limb]=std::uint64_t(value);carry=value>>64;
                }
            }
            for(std::size_t i=prime_count;i-->0;){
                if(residues[i][k]>primes[i].value/2){sum-=modulus_;break;}
                if(residues[i][k]<primes[i].value/2)break;
            }
            output[k]=sum;
        }
    }
};
} // namespace TFHEpp::modular_rns32
