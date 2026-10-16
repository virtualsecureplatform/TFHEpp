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
template<std::size_t Count=prime_count>
constexpr auto basisProduct() {
    static_assert(Count>0 && Count<=prime_count);
    MultiLimbUInt<8> value{1};
    for(std::size_t i=0;i<Count;++i)value=multiplySmall(value,primes[i].value);
    return value;
}
inline constexpr auto product=[] {
    MultiLimbUInt<8> x{1};for(auto p:primes)x=multiplySmall(x,p.value);return x;
}();

template<std::size_t Limbs,std::size_t Q,std::size_t Count=prime_count>
class CanonicalEncoder {
    static_assert(Count>0 && Count<=prime_count);
    static_assert(Q>0 && Q<=64*Limbs);
    static constexpr std::size_t words=(Q+31)/32;
    std::array<std::array<Multiplier,words>,Count> radix_;
    std::array<std::uint32_t,Count> squared_,cubed_;
public:
    CanonicalEncoder(){
        for(std::size_t i=0;i<Count;++i){
            const auto p=primes[i].value;std::uint64_t radix=1;
            const auto r=(std::uint64_t{1}<<32)%p;
            squared_[i]=r*r%p;cubed_[i]=std::uint64_t(squared_[i])*r%p;
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
    // Encode a coefficient batch after a public shift and source-level mask.
    // The result may use ordinary or Montgomery residues. SourceQ/Shift make
    // the odd-sized high tensor chunk exact for uncanonical input storage.
    template<std::size_t SourceQ=64*Limbs,std::size_t Shift=0,bool Montgomery=false>
    void residues(std::span<std::uint32_t> output,
                  std::span<const MultiLimbUInt<Limbs>> input,std::size_t prime)const {
        static_assert(Shift<=SourceQ && SourceQ<=64*Limbs);
        static_assert(sizeof(MultiLimbUInt<Limbs>)==Limbs*sizeof(std::uint64_t));
        constexpr auto active_bits=std::min(Q,SourceQ-Shift);
        constexpr auto active_words=(active_bits+31)/32;
        // Each Shoup residue is below 2p, so the wide sums stay below p*2^32.
        static_assert(active_words<(std::uint64_t{1}<<31));
        if(output.size()!=input.size())throw std::invalid_argument("RNS batch lengths differ");
        if(prime>=Count)throw std::out_of_range("RNS batch prime index");
        const auto p=std::uint32_t(primes[prime].value);
        const modular_ntt32::Montgomery arithmetic(p);
        std::size_t i=0;
#if defined(__AVX2__)
        const auto offsets=_mm256_setr_epi32(0,2*Limbs,4*Limbs,6*Limbs,8*Limbs,10*Limbs,12*Limbs,14*Limbs);
        const auto modulus=_mm256_set1_epi32(p);
        for(;i+8<=input.size();i+=8){
            auto even=_mm256_setzero_si256(),odd=even;
            const auto*data=reinterpret_cast<const int*>(input.data()+i);
            for(std::size_t w=0;w<active_words;++w){
                const auto word=Shift/32+w;
                auto digit=_mm256_i32gather_epi32(data+word,offsets,4);
                if constexpr(Shift%32){
                    digit=_mm256_srli_epi32(digit,Shift%32);
                    if(word+1<2*Limbs)digit=_mm256_or_si256(digit,
                        _mm256_slli_epi32(_mm256_i32gather_epi32(data+word+1,offsets,4),32-Shift%32));
                }
                if constexpr(active_bits%32)if(w==active_words-1)
                    digit=_mm256_and_si256(digit,_mm256_set1_epi32((std::uint32_t{1}<<(active_bits%32))-1));
                const auto value=_mm256_set1_epi32(radix_[prime][w].value);
                const auto quotient=_mm256_set1_epi32(radix_[prime][w].quotient);
                const auto high=_mm256_srli_epi64(digit,32);
                const auto eq=_mm256_srli_epi64(_mm256_mul_epu32(digit,quotient),32);
                const auto oq=_mm256_srli_epi64(_mm256_mul_epu32(high,quotient),32);
                even=_mm256_add_epi64(even,_mm256_sub_epi64(_mm256_mul_epu32(digit,value),_mm256_mul_epu32(eq,modulus)));
                odd=_mm256_add_epi64(odd,_mm256_sub_epi64(_mm256_mul_epu32(high,value),_mm256_mul_epu32(oq,modulus)));
            }
            const auto reduced=arithmetic.canonicalVector(arithmetic.reduceWideVector(even,odd));
            const auto radix=_mm256_set1_epi32(Montgomery?cubed_[prime]:squared_[prime]);
            const auto result=arithmetic.canonicalVector(arithmetic.multiplyLazyVector(reduced,radix));
            _mm256_storeu_si256(reinterpret_cast<__m256i*>(output.data()+i),result);
        }
#endif
        for(;i<input.size();++i){
            std::uint64_t sum=0;
            for(std::size_t w=0;w<active_words;++w){
                const auto bit=Shift+32*w,word=bit/64,offset=bit%64;
                auto value=input[i].limb[word]>>offset;
                if(offset>32 && word+1<Limbs)value|=input[i].limb[word+1]<<(64-offset);
                auto digit=std::uint32_t(value);
                if constexpr(active_bits%32)if(w==active_words-1)digit&=(std::uint32_t{1}<<(active_bits%32))-1;
                sum+=radix_[prime][w].lazy(digit,p);
            }
            const auto value=std::uint32_t(sum%p);
            output[i]=Montgomery?arithmetic.encode(value):value;
        }
    }

};

// Canonical residues of signed 64-bit coefficients, including INT64_MIN.
// The fixed RNS primes exceed 2^31/3. For magnitude t<=2^63 and m<2^32,
// t+m*p<2^64 and REDC(t)<2^31+p<4p. Fold that wider result below p,
// then a multiplication by R^2 removes REDC's inverse-radix factor.
template<std::size_t Count=prime_count>
class SignedEncoder {
    static_assert(Count>0 && Count<=prime_count);
    static_assert([]{
        for(std::size_t p=0;p<Count;++p)
            if(3*primes[p].value<=(std::uint64_t{1}<<31))return false;
        return true;
    }());
    std::vector<modular_ntt32::Montgomery> arithmetic_;
    std::array<std::uint32_t,Count> squared_;
public:
    SignedEncoder(){
        arithmetic_.reserve(Count);
        for(std::size_t p=0;p<Count;++p){
            arithmetic_.emplace_back(primes[p].value);
            const auto radix=(std::uint64_t{1}<<32)%primes[p].value;
            squared_[p]=radix*radix%primes[p].value;
        }
    }
    void residues(std::span<std::uint32_t> output,std::span<const std::int64_t> input,
                  std::size_t prime,bool batch=true)const {
        if(output.size()!=input.size())throw std::invalid_argument("signed RNS batch lengths differ");
        if(prime>=Count)throw std::out_of_range("signed RNS batch prime index");
        const auto p=std::uint32_t(primes[prime].value);
        [[maybe_unused]] const auto&arithmetic=arithmetic_[prime];
        std::size_t i=0;
#if defined(__AVX2__)
        if(batch){
            const auto zero=_mm256_setzero_si256(),modulus=_mm256_set1_epi32(p);
            const auto twice=_mm256_set1_epi32(2*p),square=_mm256_set1_epi32(squared_[prime]);
            for(;i+8<=input.size();i+=8){
                const auto a=_mm256_loadu_si256(reinterpret_cast<const __m256i*>(input.data()+i));
                const auto b=_mm256_loadu_si256(reinterpret_cast<const __m256i*>(input.data()+i+4));
                auto even=_mm256_permute4x64_epi64(_mm256_unpacklo_epi64(a,b),_MM_SHUFFLE(3,1,2,0));
                auto odd=_mm256_permute4x64_epi64(_mm256_unpackhi_epi64(a,b),_MM_SHUFFLE(3,1,2,0));
                const auto eneg=_mm256_cmpgt_epi64(zero,even),oneg=_mm256_cmpgt_epi64(zero,odd);
                even=_mm256_sub_epi64(_mm256_xor_si256(even,eneg),eneg);
                odd=_mm256_sub_epi64(_mm256_xor_si256(odd,oneg),oneg);
                auto value=arithmetic.reduceWideVector(even,odd);
                value=_mm256_min_epu32(value,_mm256_sub_epi32(value,twice));
                value=arithmetic.canonicalVector(value);
                value=arithmetic.canonicalVector(arithmetic.multiplyLazyVector(value,square));
                const auto negative=_mm256_or_si256(_mm256_srli_epi64(eneg,32),_mm256_slli_epi64(oneg,32));
                value=_mm256_blendv_epi8(value,_mm256_sub_epi32(modulus,value),negative);
                _mm256_storeu_si256(reinterpret_cast<__m256i*>(output.data()+i),arithmetic.canonicalVector(value));
            }
        }
#else
        (void)batch;
#endif
        for(;i<input.size();++i){
            const auto v=input[i];
            const auto magnitude=v<0?std::uint64_t{0}-std::uint64_t(v):std::uint64_t(v);
            const auto reduced=magnitude%p;
            output[i]=v<0&&reduced?p-reduced:reduced;
        }
    }
};

// Garner mixed radix with exact centering. The digits of floor(M/2) are
// floor(p_i/2); compare from the highest radix to determine the sign without
// floating point or a full-width reconstruction. Accumulate only the target
// torus limbs, then subtract M modulo that same torus when negative.
template<std::size_t Limbs,std::size_t Count=prime_count>
class CenteredCRT {
    static_assert(Count>0 && Count<=prime_count);
    using T=MultiLimbUInt<Limbs>;
    std::array<std::array<Multiplier,Count>,Count> inverse_{};
    std::array<std::array<std::uint32_t,Count>,Count> inverse_mont_{};
    std::vector<modular_ntt32::Montgomery> arithmetic_;
    std::array<T,Count> prefixes_;
    T modulus_;
public:
    CenteredCRT(){
        T product{1};arithmetic_.reserve(Count);
        for(std::size_t i=0;i<Count;++i){
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
    void reconstruct(std::array<std::span<std::uint32_t>,Count> residues,std::span<T> output)const {
        for(auto row:residues)if(row.size()!=output.size())throw std::invalid_argument("direct CRT lengths differ");
        for(std::size_t i=1;i<Count;++i){
            const auto p=std::uint32_t(primes[i].value);
            for(std::size_t j=0;j<i;++j){
                std::size_t k=0;
#if defined(__AVX2__)
                const auto prime=_mm256_set1_epi32(p);
                const auto fixed=_mm256_set1_epi32(inverse_mont_[i][j]);
                const auto vector_end=output.size()&~std::size_t{7};
                for(;k<vector_end;k+=8){
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
            for(std::size_t i=0;i<Count;++i){
                unsigned __int128 carry=0;
                for(std::size_t limb=0;limb<Limbs;++limb){
                    const auto value=static_cast<unsigned __int128>(prefixes_[i].limb[limb])*residues[i][k]+sum.limb[limb]+carry;
                    sum.limb[limb]=std::uint64_t(value);carry=value>>64;
                }
            }
            for(std::size_t i=Count;i-->0;){
                if(residues[i][k]>primes[i].value/2){sum-=modulus_;break;}
                if(residues[i][k]<primes[i].value/2)break;
            }
            output[k]=sum;
        }
    }
};
// Direct CRT with an integer-certified centering quotient. For gamma_i<p_i,
// gamma_i*floor(2^64/p_i) underestimates 2^64*gamma_i/p_i by less than
// gamma_i. Thus the sum and sum+sum(gamma_i) enclose the exact quotient.
// Round both ends; if they disagree, use the original exact Garner method.
// This avoids quadratic mixed-radix conversion in the usual case without
// floating point, a bound on the represented integer, or a precision loss.
template<std::size_t Limbs,std::size_t Count=prime_count>
class DirectCenteredCRT {
    static_assert(Count>0 && Count<=prime_count);
    using T=MultiLimbUInt<Limbs>;
    static constexpr std::size_t block_size=256;
    CenteredCRT<Limbs,Count> fallback_;
    std::array<T,Count> basis_;
    std::array<T,Count+1> multiples_{};
    std::array<Multiplier,Count> inverse_;
    std::array<std::uint32_t,Count> inverse_mont_;
    std::array<std::uint64_t,Count> reciprocal_;
    std::vector<modular_ntt32::Montgomery> arithmetic_;
public:
    DirectCenteredCRT(){
        const auto full_modulus=basisProduct<Count>();
        T modulus{};
        std::copy_n(full_modulus.limb.begin(),std::min(Limbs,full_modulus.limbs),modulus.limb.begin());
        for(std::size_t q=1;q<=Count;++q)multiples_[q]=multiples_[q-1]+modulus;
        arithmetic_.reserve(Count);
        for(std::size_t i=0;i<Count;++i){
            MultiLimbUInt<8> basis{1};
            for(std::size_t j=0;j<Count;++j)if(j!=i)basis=multiplySmall(basis,primes[j].value);
            basis_[i]={};
            std::copy_n(basis.limb.begin(),std::min(Limbs,basis.limbs),basis_[i].limb.begin());
            const auto p=std::uint32_t(primes[i].value);
            std::uint64_t residue=0;
            for(std::size_t limb=basis.limbs;limb-->0;)
                residue=((static_cast<unsigned __int128>(residue)<<64)|basis.limb[limb])%p;
            const auto inverse=modular_ntt::invert(residue,p);
            inverse_[i]=Multiplier(inverse,p);
            arithmetic_.emplace_back(p);
            inverse_mont_[i]=arithmetic_[i].encode(inverse);
            reciprocal_[i]=(static_cast<unsigned __int128>(1)<<64)/p;
        }
    }
    // Canonical input rows are consumed (the fallback may overwrite entries).
    // Return the number of coefficients requiring exact Garner conversion.
    std::size_t reconstruct(std::array<std::span<std::uint32_t>,Count> residues,
                            std::span<T> output)const {
        for(auto row:residues)if(row.size()!=output.size())throw std::invalid_argument("direct CRT lengths differ");
        std::size_t fallbacks=0;
        std::array<std::array<std::uint32_t,block_size>,Count> gamma;
        for(std::size_t begin=0;begin<output.size();begin+=block_size){
            const auto length=std::min(block_size,output.size()-begin);
            for(std::size_t p=0;p<Count;++p){
                std::size_t k=0;
#if defined(__AVX2__)
                const auto fixed=_mm256_set1_epi32(inverse_mont_[p]);
                for(;k+8<=length;k+=8){
                    const auto input=_mm256_loadu_si256(reinterpret_cast<const __m256i*>(residues[p].data()+begin+k));
                    const auto value=arithmetic_[p].canonicalVector(arithmetic_[p].multiplyLazyVector(input,fixed));
                    _mm256_storeu_si256(reinterpret_cast<__m256i*>(gamma[p].data()+k),value);
                }
#endif
                for(;k<length;++k)gamma[p][k]=inverse_[p].apply(residues[p][begin+k],primes[p].value);
            }
            for(std::size_t k=0;k<length;++k){
                T sum{};
                unsigned __int128 fraction=0;
                std::uint64_t error=0;
                for(std::size_t p=0;p<Count;++p){
                    const auto digit=gamma[p][k];
                    // digit<p, so this product is strictly below 2^64.
                    fraction+=std::uint64_t(digit)*reciprocal_[p];
                    error+=digit;
                    unsigned __int128 carry=0;
                    for(std::size_t limb=0;limb<Limbs;++limb){
                        const auto value=static_cast<unsigned __int128>(basis_[p].limb[limb])*digit+sum.limb[limb]+carry;
                        sum.limb[limb]=std::uint64_t(value);carry=value>>64;
                    }
                }
                constexpr auto half=static_cast<unsigned __int128>(1)<<63;
                const auto lower=std::size_t((fraction+half)>>64);
                const auto upper=std::size_t((fraction+half+error)>>64);
                if(lower==upper)output[begin+k]=sum-multiples_[lower];
                else {
                    std::array<std::span<std::uint32_t>,Count> single;
                    for(std::size_t p=0;p<Count;++p)single[p]=residues[p].subspan(begin+k,1);
                    fallback_.reconstruct(single,output.subspan(begin+k,1));
                    ++fallbacks;
                }
            }
        }
        return fallbacks;
    }
};
} // namespace TFHEpp::modular_rns32
