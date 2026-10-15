#pragma once
#include "compact_dd.hpp"
#include <modular_rns32.hpp>

#ifdef TFHEPP_HAS_EXTENDED_MULTILIMB_PARAMS
namespace TFHEpp {
// Exact encrypted-gadget multiplication using wide canonical key chunks.
// The primary gadget and its balanced input decomposition stay unchanged.
template<class P,std::uint32_t Q,std::uint32_t Bits,std::uint32_t ChunkBits=80>
class CKKSPreparedGadgetRNS32 {
    static_assert(P::k==1 && is_multilimb_uint_v<typename P::T>);
    static_assert(Bits>0 && Bits<=64 && Q<=P::T::digits);
    static constexpr auto rows=(Q+Bits-1)/Bits;
    static constexpr auto chunk_bits=std::min(Q,ChunkBits);
    static_assert(chunk_bits>0);
    static constexpr auto chunks=(Q+chunk_bits-1)/chunk_bits;
    static_assert(rows*Bits<=P::T::digits);
    static_assert(chunk_bits+Bits+P::nbit+std::bit_width(rows)<512);
    static constexpr auto count=[] {
        auto bound=MultiLimbUInt<8>{1}<<(chunk_bits+Bits);
        bound=modular_rns32::multiplySmall(bound,P::n);
        bound=modular_rns32::multiplySmall(bound,rows);
        MultiLimbUInt<8> product{1};
        for(std::size_t p=0;p<modular_rns32::prime_count;++p){
            product=modular_rns32::multiplySmall(product,modular_rns32::primes[p].value);
            if(product>bound)return p+1;
        }
        return std::size_t{0};
    }();
    using Small=MultiLimbUInt<(chunk_bits+Bits+P::nbit+std::bit_width(rows)+63)/64>;
    static_assert(count>0,"full gadget convolution exceeds centered CRT range");
    static const auto &plans(){
        static const auto value=[] {
            std::vector<modular_ntt32::NegacyclicNTTPlan> out;
            for(std::size_t p=0;p<count;++p)out.emplace_back(P::n,modular_rns32::primes[p]);
            return out;
        }();return value;
    }
    std::array<std::vector<std::uint32_t>,count> spectra_;
public:
    static constexpr auto prime_count=count;
    template<class GadgetRows> explicit CKKSPreparedGadgetRNS32(const GadgetRows&key){
        static_assert(GadgetRows::rows==rows);
        static const modular_rns32::CanonicalEncoder<P::T::limbs,chunk_bits,count> encoder;
        (void)plans();for(auto &row:spectra_)row.resize(chunks*rows*2*P::n);
#pragma omp parallel for collapse(4) schedule(static)
        for(std::size_t p=0;p<count;++p)for(std::size_t chunk=0;chunk<chunks;++chunk)for(std::size_t r=0;r<rows;++r)for(std::size_t c=0;c<2;++c){
            const auto &arithmetic=plans()[p].arithmetic();
            auto values=std::span(spectra_[p]).subspan(((chunk*rows+r)*2+c)*P::n,P::n);
            for(std::size_t i=0;i<P::n;++i)values[i]=arithmetic.encode(encoder.residue(ckks_detail::reduceToLevel<P,Q>(key.ciphertexts[r][c][i])>>(chunk*chunk_bits),p));
            plans()[p].forward(values);
        }
    }
    std::size_t byteSize()const{return count*chunks*rows*2*P::n*sizeof(std::uint32_t);}
    void apply(TRLWE<P>&result,const Polynomial<P>&input)const{
        auto digits=std::make_unique<std::array<Polynomial<P>,rows>>();
        ckks_detail::activeBaseDecomposePolynomialRows<P,Q,Bits,rows>(*digits,input);
        std::array<bool,rows> active{};
        for(std::size_t r=0;r<rows;++r)active[r]=std::any_of((*digits)[r].begin(),(*digits)[r].end(),[](const auto&v){return v!=typename P::T{};});
        std::array<std::vector<std::uint32_t>,count> transformed;
        for(auto &row:transformed)row.resize(rows*P::n);
#pragma omp parallel for collapse(2) schedule(static)
        for(std::size_t p=0;p<count;++p)for(std::size_t r=0;r<rows;++r)if(active[r]){
            auto values=std::span(transformed[p]).subspan(r*P::n,P::n);
            const auto modulus=modular_rns32::primes[p].value;
            for(std::size_t i=0;i<P::n;++i){
                const auto v=std::bit_cast<std::int64_t>((*digits)[r][i].limb[0]);
                const auto magnitude=v<0?std::uint64_t{0}-std::uint64_t(v):std::uint64_t(v);
                const auto reduced=magnitude%modulus;values[i]=v<0&&reduced?modulus-reduced:reduced;
            }
            plans()[p].forward(values);
        }
        digits.reset();
        std::array<std::array<std::vector<std::uint32_t>,count>,chunks*2> sums;
        for(auto &c:sums)for(auto &row:c)row.assign(P::n,0);
#pragma omp parallel for collapse(3) schedule(static)
        for(std::size_t c=0;c<chunks*2;++c)for(std::size_t p=0;p<count;++p)for(std::size_t begin=0;begin<P::n;begin+=1024){
            const auto &arithmetic=plans()[p].arithmetic();auto*out=sums[c][p].data();
            const auto end=std::min(begin+1024,std::size_t(P::n));
            for(std::size_t r=0;r<rows;++r)if(active[r]){
                const auto*x=transformed[p].data()+r*P::n,*y=spectra_[p].data()+((c/2*rows+r)*2+c%2)*P::n;std::size_t i=begin;
#if defined(__AVX2__)
                for(;i+8<=end;i+=8){
                    const auto v=arithmetic.accumulateVector(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(out+i)),_mm256_loadu_si256(reinterpret_cast<const __m256i*>(x+i)),_mm256_loadu_si256(reinterpret_cast<const __m256i*>(y+i)));
                    _mm256_storeu_si256(reinterpret_cast<__m256i*>(out+i),v);
                }
#endif
                for(;i<end;++i){const auto v=out[i]+arithmetic.multiplyLazy(x[i],y[i]);out[i]=v>=2*arithmetic.modulus()?v-2*arithmetic.modulus():v;}
            }
        }
#pragma omp parallel for collapse(2) schedule(static)
        for(std::size_t c=0;c<chunks*2;++c)for(std::size_t p=0;p<count;++p){
            for(auto &v:sums[c][p])v=plans()[p].arithmetic().canonical(v);
            plans()[p].inverse(sums[c][p]);
        }
        static const modular_rns32::CenteredCRT<Small::limbs,count> crt;
        for(auto&c:result)c.fill({});
#pragma omp parallel
        {
            std::vector<Small> recovered(512);
#pragma omp for collapse(2) schedule(static)
            for(std::size_t c=0;c<2;++c)for(std::size_t begin=0;begin<P::n;begin+=512){
                const auto length=std::min(std::size_t{512},std::size_t(P::n)-begin);
                for(std::size_t chunk=0;chunk<chunks;++chunk){
                    std::array<std::span<std::uint32_t>,count> residues;
                    for(std::size_t p=0;p<count;++p)residues[p]=std::span(sums[chunk*2+c][p]).subspan(begin,length);
                    crt.reconstruct(residues,std::span(recovered).first(length));
                    for(std::size_t i=0;i<length;++i){
                        typename P::T value{};
                        if constexpr(Small::limbs<P::T::limbs)
                            if(recovered[i].limb.back()>>63)value.limb.fill(~std::uint64_t{0});
                        std::copy_n(recovered[i].limb.begin(),std::min(Small::limbs,P::T::limbs),value.limb.begin());
                        result[c][begin+i]=ckks_detail::reduceToLevel<P,Q>(result[c][begin+i]+(value<<(chunk*chunk_bits)));
                    }
                }
            }
        }
    }
};
} // namespace TFHEpp
#endif
