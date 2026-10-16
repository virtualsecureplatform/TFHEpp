#pragma once
#include <modular_rns32.hpp>

namespace TFHEpp::modular_ntt32 {
// Diagnostic access shares every root/twist table with the canonical plan.
// No production call selects these kernels.
struct LazyNTTAccess {
    static void normalize(const NegacyclicNTTPlan&p,std::span<std::uint32_t>x){
        std::size_t i=0;const auto modulus=p.arithmetic_.modulus();
#if defined(__AVX2__)
        for(;i+8<=x.size();i+=8){
            auto v=_mm256_loadu_si256(reinterpret_cast<const __m256i*>(x.data()+i));
            v=p.arithmetic_.canonicalVector(v);
            _mm256_storeu_si256(reinterpret_cast<__m256i*>(x.data()+i),v);
        }
#endif
        for(;i<x.size();++i)if(x[i]>=modulus)x[i]-=modulus;
    }
    static void range(std::span<const std::uint32_t>x,std::uint64_t bound,bool check){
        if(check)for(auto v:x)if(v>=bound)throw std::runtime_error("lazy NTT range invariant failed");
    }
#if defined(__AVX2__)
    template<std::size_t Length>
    static void shortDIF(const NegacyclicNTTPlan&p,std::span<std::uint32_t>x,
        const std::vector<std::uint32_t>&twiddles){
        const auto prime=_mm256_set1_epi32(p.arithmetic_.modulus());
        const auto twice=_mm256_add_epi32(prime,prime);
        const auto fixed=[&]{
            if constexpr(Length==4)return _mm256_setr_epi32(twiddles[0],twiddles[1],twiddles[0],twiddles[1],twiddles[0],twiddles[1],twiddles[0],twiddles[1]);
            else if constexpr(Length==8)return _mm256_broadcastsi128_si256(_mm_loadu_si128(reinterpret_cast<const __m128i*>(twiddles.data())));
            else return _mm256_setzero_si256();
        }();
        for(std::size_t i=0;i<x.size();i+=8){
            auto input=_mm256_loadu_si256(reinterpret_cast<const __m256i*>(x.data()+i));
            if constexpr(Length==2){
                input=p.arithmetic_.canonicalVector(input);
                const auto swapped=_mm256_shuffle_epi32(input,_MM_SHUFFLE(2,3,0,1));
                const auto sum=p.arithmetic_.canonicalVector(_mm256_add_epi32(input,swapped));
                const auto diff=p.arithmetic_.canonicalVector(_mm256_sub_epi32(_mm256_add_epi32(swapped,prime),input));
                _mm256_storeu_si256(reinterpret_cast<__m256i*>(x.data()+i),_mm256_blend_epi32(sum,diff,0xaa));
            }else{
                const auto even=[&]{if constexpr(Length==4)return _mm256_shuffle_epi32(input,_MM_SHUFFLE(1,0,1,0));else return _mm256_permute2x128_si256(input,input,0x00);}();
                const auto odd=[&]{if constexpr(Length==4)return _mm256_shuffle_epi32(input,_MM_SHUFFLE(3,2,3,2));else return _mm256_permute2x128_si256(input,input,0x11);}();
                auto sum=_mm256_add_epi32(even,odd);sum=_mm256_min_epu32(sum,_mm256_sub_epi32(sum,twice));
                const auto diff=p.arithmetic_.multiplyLazyVector(_mm256_sub_epi32(_mm256_add_epi32(even,twice),odd),fixed);
                _mm256_storeu_si256(reinterpret_cast<__m256i*>(x.data()+i),_mm256_blend_epi32(sum,diff,Length==4?0xcc:0xf0));
            }
        }
    }
#endif
    static void dif(const NegacyclicNTTPlan&p,std::span<std::uint32_t>x,bool check){
        const auto modulus=p.arithmetic_.modulus();auto stage=p.forward_.size();
        for(std::size_t length=p.size_;length>=2;length>>=1){
            --stage;
#if defined(__AVX2__)
            if(p.vector_short_stages_&&p.size_>=8&&length<=8){
                if(length==8)shortDIF<8>(p,x,p.forward_[stage]);
                else if(length==4)shortDIF<4>(p,x,p.forward_[stage]);
                else shortDIF<2>(p,x,p.forward_[stage]);
                range(x,std::uint64_t(length==2?1:2)*modulus,check);continue;
            }
#endif
            const auto half=length/2;
            for(std::size_t block=0;block<p.size_;block+=length){
                std::size_t j=0;
#if defined(__AVX2__)
                const auto twice=_mm256_set1_epi32(2*modulus);
                for(;j+8<=half;j+=8){
                    const auto even=_mm256_loadu_si256(reinterpret_cast<const __m256i*>(x.data()+block+j));
                    const auto odd=_mm256_loadu_si256(reinterpret_cast<const __m256i*>(x.data()+block+j+half));
                    auto sum=_mm256_add_epi32(even,odd);
                    sum=_mm256_min_epu32(sum,_mm256_sub_epi32(sum,twice));
                    const auto difference=_mm256_sub_epi32(_mm256_add_epi32(even,twice),odd);
                    const auto product=p.arithmetic_.multiplyLazyVector(difference,_mm256_loadu_si256(reinterpret_cast<const __m256i*>(p.forward_[stage].data()+j)));
                    _mm256_storeu_si256(reinterpret_cast<__m256i*>(x.data()+block+j),sum);
                    _mm256_storeu_si256(reinterpret_cast<__m256i*>(x.data()+block+j+half),product);
                }
#endif
                for(;j<half;++j){
                    const auto even=x[block+j],odd=x[block+j+half],sum=even+odd;
                    x[block+j]=sum>=2*modulus?sum-2*modulus:sum;
                    x[block+j+half]=p.arithmetic_.multiplyLazy(even+2*modulus-odd,p.forward_[stage][j]);
                }
            }
            range(x,std::uint64_t(2)*modulus,check);
        }
#if defined(__AVX2__)
        if(p.vector_short_stages_&&p.size_>=8)return;
#endif
        normalize(p,x);
    }
    static void dit(const NegacyclicNTTPlan&p,std::span<std::uint32_t>x,bool inverse,bool permute,bool check){
        if(permute)for(auto[a,b]:p.swaps_)std::swap(x[a],x[b]);
        const auto modulus=p.arithmetic_.modulus();const auto&stages=inverse?p.inverse_:p.forward_;
        std::size_t stage=0;
        for(std::size_t length=2;length<=p.size_;length<<=1,++stage){
#if defined(__AVX2__)
            if(p.vector_short_stages_&&p.size_>=8&&length<=8){
                p.shortStage<false>(x,length,stages[stage]);range(x,modulus,check);continue;
            }
#endif
            const auto half=length/2;
            for(std::size_t block=0;block<p.size_;block+=length){
                std::size_t j=0;
#if defined(__AVX2__)
                const auto twice=_mm256_set1_epi32(2*modulus);
                for(;j+8<=half;j+=8){
                    auto even=_mm256_loadu_si256(reinterpret_cast<const __m256i*>(x.data()+block+j));
                    even=_mm256_min_epu32(even,_mm256_sub_epi32(even,twice));
                    const auto odd=p.arithmetic_.multiplyLazyVector(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(x.data()+block+j+half)),_mm256_loadu_si256(reinterpret_cast<const __m256i*>(stages[stage].data()+j)));
                    _mm256_storeu_si256(reinterpret_cast<__m256i*>(x.data()+block+j),_mm256_add_epi32(even,odd));
                    _mm256_storeu_si256(reinterpret_cast<__m256i*>(x.data()+block+j+half),_mm256_sub_epi32(_mm256_add_epi32(even,twice),odd));
                }
#endif
                for(;j<half;++j){
                    auto even=x[block+j];if(even>=2*modulus)even-=2*modulus;
                    const auto odd=p.arithmetic_.multiplyLazy(x[block+j+half],stages[stage][j]);
                    x[block+j]=even+odd;x[block+j+half]=even+2*modulus-odd;
                }
            }
            range(x,std::uint64_t(4)*modulus,check);
        }
    }
    static void forward(const NegacyclicNTTPlan&p,std::span<std::uint32_t>x,bool bit_reversed,bool check=false){
        if(x.size()!=p.size_)throw std::invalid_argument("NTT input length mismatch");
        p.multiplyRows(x,p.twist_);
        if(bit_reversed)dif(p,x,check);
        else{
            dit(p,x,false,true,check);
            // DIT may return <4p. Reduce to <2p before canonicalization.
            for(auto&v:x)if(v>=2*p.arithmetic_.modulus())v-=2*p.arithmetic_.modulus();
            normalize(p,x);
        }
        range(x,p.arithmetic_.modulus(),check);
    }
    static void inverse(const NegacyclicNTTPlan&p,std::span<std::uint32_t>x,bool bit_reversed,bool check=false){
        if(x.size()!=p.size_)throw std::invalid_argument("NTT input length mismatch");
        dit(p,x,true,!bit_reversed,check);
        // x<4p and each twiddle<p imply t<4p^2<p*2^32, so REDC<2p.
        p.multiplyRows(x,p.untwist_);range(x,p.arithmetic_.modulus(),check);
    }
};
}
