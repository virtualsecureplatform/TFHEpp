#pragma once
// Rejected CPU experiment: component savings did not improve retained FBS.
// Kept here for reproducible arithmetic checks and diagnostic benchmarks.
#include <ckks/tensor_rns32.hpp>
#ifdef TFHEPP_HAS_EXTENDED_MULTILIMB_PARAMS
namespace TFHEpp::ckks_experiment {
template<class P,std::uint32_t Q,bool CompactCRT=true,bool DirectCRT=true,bool BatchEncode=true,bool BitReversed=true,bool PublicStructure=true>
void tensor(TRLWE3<P>&result,const TRLWE<P>&lhs,const TRLWE<P>&rhs)
{
    static_assert(P::k==1 && is_multilimb_uint_v<typename P::T>);
    static_assert(Q>0 && Q<=P::T::digits && Q<490);
    constexpr auto split=(Q+1)/2;
    static_assert(2*split+P::nbit+3<512);
    constexpr auto count=[] {
        MultiLimbUInt<8> product{1};
        // Four products can contribute to a cross-component diagonal; signed
        // recovery needs another factor two. All chunk coefficients <2^split.
        const auto bound=MultiLimbUInt<8>{1}<<(2*split+P::nbit+3);
        for(std::size_t p=0;p<modular_rns32::prime_count;++p){
            product=modular_rns32::multiplySmall(product,modular_rns32::primes[p].value);
            if(product>bound)return p+1;
        }
        return std::size_t{0};
    }();
    static_assert(count>0,"tensor exceeds exact RNS basis range");
    // Ciphertext coefficients and operand identity are public. Zero component
    // tests stop at the first active coefficient on ordinary dense products.
    const auto active=[](const Polynomial<P>&poly){
        return std::any_of(poly.begin(),poly.end(),[](const auto&v){
            return ckks_detail::reduceToLevel<P,Q>(v)!=typename P::T{};
        });
    };
    std::array<bool,2> left_active{true,true},right_active{true,true};
    const bool square=PublicStructure && std::addressof(lhs)==std::addressof(rhs);
    if constexpr(PublicStructure){
        for(std::size_t c=0;c<2;++c)left_active[c]=active(lhs[c]);
        for(std::size_t c=0;c<2;++c)right_active[c]=square?left_active[c]:active(rhs[c]);
    }
    if((!left_active[0]&&!left_active[1])||(!right_active[0]&&!right_active[1])){
        for(auto &component:result)component.fill({});
        return;
    }
    // Keep the dense non-aliasing kernel free of dynamic component flags.
    if constexpr(PublicStructure){
        if(!square&&left_active[0]&&left_active[1]&&right_active[0]&&right_active[1]){
            tensor<P,Q,CompactCRT,DirectCRT,BatchEncode,BitReversed,false>(result,lhs,rhs);
            return;
        }
    }
    static const auto plans=[] {
        std::vector<modular_ntt32::NegacyclicNTTPlan> out;
        for(std::size_t p=0;p<count;++p)out.emplace_back(P::n,modular_rns32::primes[p]);
        return out;
    }();
    static const modular_rns32::CanonicalEncoder<P::T::limbs,split,count> encoder;
    using Rows=std::array<std::vector<std::uint32_t>,count>;
    std::array<Rows,4> left,right;
#pragma omp parallel for collapse(2) schedule(static)
    for(std::size_t row=0;row<4;++row)for(std::size_t p=0;p<count;++p){
        const auto &arithmetic=plans[p].arithmetic();
        const auto encode=[&]<bool Montgomery>(auto&out,const auto&input){
            if constexpr(BatchEncode){
                if(row%2){
                    encoder.template residues<Q,split,Montgomery>(out,input,p);
                }else{
                    encoder.template residues<Q,0,Montgomery>(out,input,p);
                }
            }else for(std::size_t i=0;i<P::n;++i){
                auto value=ckks_detail::reduceToLevel<P,Q>(input[i]);
                if(row%2)value>>=split;
                const auto residue=encoder.residue(value,p);
                if constexpr(Montgomery)out[i]=arithmetic.encode(residue);
                else out[i]=residue;
            }
        };
        // Preserve the original dense allocation/encoding order before either
        // NTT touches its twiddle tables and work buffers.
        if(!PublicStructure||left_active[row/2])left[row][p].resize(P::n);
        if(!PublicStructure||right_active[row/2])right[row][p].resize(P::n);
        if(!PublicStructure||left_active[row/2])encode.template operator()<false>(left[row][p],lhs[row/2]);
        if((!PublicStructure||right_active[row/2])&&(!PublicStructure||!square))encode.template operator()<true>(right[row][p],rhs[row/2]);
        if(!PublicStructure||left_active[row/2]){
            if constexpr(BitReversed)plans[p].forwardBitReversed(left[row][p]);
            else plans[p].forward(left[row][p]);
        }
        if(!PublicStructure||right_active[row/2]){
            if(PublicStructure&&square){
                // NTT is linear over the residue field. Convert its completed
                // ordinary spectrum to Montgomery form instead of repeating it.
                const auto squared_radix=arithmetic.encode(arithmetic.encode(1));
                std::size_t i=0;
#if defined(__AVX2__)
                const auto factor=_mm256_set1_epi32(squared_radix);
                for(;i+8<=P::n;i+=8){
                    const auto value=arithmetic.canonicalVector(arithmetic.multiplyLazyVector(
                        _mm256_loadu_si256(reinterpret_cast<const __m256i*>(left[row][p].data()+i)),factor));
                    _mm256_storeu_si256(reinterpret_cast<__m256i*>(right[row][p].data()+i),value);
                }
#endif
                for(;i<P::n;++i)right[row][p][i]=arithmetic.multiply(left[row][p][i],squared_radix);
            }else{
                if constexpr(BitReversed)plans[p].forwardBitReversed(right[row][p]);
                else plans[p].forward(right[row][p]);
            }
        }
    }
    for(auto &component:result)component.fill({});
    // Row mapping matches CKKS tensor convention: cross, constant, square.
    constexpr std::array<std::array<int,2>,3> components{{{0,1},{1,1},{0,0}}};
    const std::array<bool,3> output_active{
        (left_active[0]&&right_active[1])||(left_active[1]&&right_active[0]),
        left_active[1]&&right_active[1],left_active[0]&&right_active[0]};
    const auto diagonal_product=[&]<std::size_t Diagonal>(){
        // Only Q-Diagonal*split low bits survive the shift and reduction.
        // CenteredCRT reconstructs modulo its output width, even when the
        // complete CRT product is wider. No sign extension is needed: every
        // omitted bit would lie at or above Q after shifting.
        constexpr auto output_bits=CompactCRT?Q-Diagonal*split:P::T::digits;
        using Recovered=MultiLimbUInt<std::max<std::size_t>(1,(output_bits+63)/64)>;
        using CRT=std::conditional_t<DirectCRT && count>=8,
            modular_rns32::DirectCenteredCRT<Recovered::limbs,count>,
            modular_rns32::CenteredCRT<Recovered::limbs,count>>;
        static const CRT crt;
        std::array<Rows,3> sums;
#pragma omp parallel for collapse(2) schedule(static)
        for(std::size_t c=0;c<3;++c)for(std::size_t p=0;p<count;++p)if(!PublicStructure||output_active[c]){
            auto &out=sums[c][p];out.assign(P::n,0);
            const auto &arithmetic=plans[p].arithmetic();
            const auto accumulate=[&](int a,int b){
                if constexpr(PublicStructure)
                    if(!left_active[a/2]||!right_active[b/2])return;
                const auto *x=left[a][p].data(),*y=right[b][p].data();std::size_t i=0;
#if defined(__AVX2__)
                for(;i+8<=P::n;i+=8){
                    const auto value=arithmetic.accumulateVector(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(out.data()+i)),_mm256_loadu_si256(reinterpret_cast<const __m256i*>(x+i)),_mm256_loadu_si256(reinterpret_cast<const __m256i*>(y+i)));
                    _mm256_storeu_si256(reinterpret_cast<__m256i*>(out.data()+i),value);
                }
#endif
                for(;i<P::n;++i){const auto value=out[i]+arithmetic.multiplyLazy(x[i],y[i]);out[i]=value>=2*arithmetic.modulus()?value-2*arithmetic.modulus():value;}
            };
            for(std::size_t high=0;high<=Diagonal;++high){
                const auto a=2*components[c][0]+high,b=2*components[c][1]+Diagonal-high;
                accumulate(a,b);if(c==0)accumulate(2*components[c][1]+high,2*components[c][0]+Diagonal-high);
            }
            for(auto &value:out)value=arithmetic.canonical(value);
            if constexpr(BitReversed)plans[p].inverseBitReversed(out);
            else plans[p].inverse(out);
        }
#pragma omp parallel
        {
            std::vector<Recovered> recovered(512);
#pragma omp for collapse(2) schedule(static)
            for(std::size_t c=0;c<3;++c)for(std::size_t begin=0;begin<P::n;begin+=512)if(!PublicStructure||output_active[c]){
                const auto length=std::min(std::size_t{512},std::size_t(P::n)-begin);
                std::array<std::span<std::uint32_t>,count> residues;
                for(std::size_t p=0;p<count;++p)residues[p]=std::span(sums[c][p]).subspan(begin,length);
                crt.reconstruct(residues,std::span(recovered).first(length));
                for(std::size_t i=0;i<length;++i){
                    typename P::T value{};
                    std::copy_n(recovered[i].limb.begin(),Recovered::limbs,value.limb.begin());
                    result[c][begin+i]=ckks_detail::reduceToLevel<P,Q>(result[c][begin+i]+(value<<(Diagonal*split)));
                }
            }
        }
    };
    diagonal_product.template operator()<0>();
    diagonal_product.template operator()<1>();
}
} // namespace TFHEpp::ckks_experiment
#endif
