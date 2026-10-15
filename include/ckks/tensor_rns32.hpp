#pragma once
#include "ckks.hpp"
#include <modular_rns32.hpp>

#ifdef TFHEPP_HAS_EXTENDED_MULTILIMB_PARAMS
namespace TFHEpp {
// Exact unscaled tensor modulo 2^Q. Split canonical coefficients at ceil(Q/2):
// the high*high contribution vanishes, leaving two convolution diagonals.
template<class P,std::uint32_t Q>
void CKKSTensorProductRNS32(TRLWE3<P>&result,const TRLWE<P>&lhs,const TRLWE<P>&rhs)
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
    static const auto plans=[] {
        std::vector<modular_ntt32::NegacyclicNTTPlan> out;
        for(std::size_t p=0;p<count;++p)out.emplace_back(P::n,modular_rns32::primes[p]);
        return out;
    }();
    static const modular_rns32::CanonicalEncoder<P::T::limbs,split,count> encoder;
    static const modular_rns32::CenteredCRT<P::T::limbs,count> crt;
    using Rows=std::array<std::vector<std::uint32_t>,count>;
    std::array<Rows,4> left,right;
#pragma omp parallel for collapse(2) schedule(static)
    for(std::size_t row=0;row<4;++row)for(std::size_t p=0;p<count;++p){
        left[row][p].resize(P::n);right[row][p].resize(P::n);
        const auto &arithmetic=plans[p].arithmetic();
        for(std::size_t i=0;i<P::n;++i){
            auto a=ckks_detail::reduceToLevel<P,Q>(lhs[row/2][i]);
            auto b=ckks_detail::reduceToLevel<P,Q>(rhs[row/2][i]);
            if(row%2){a>>=split;b>>=split;}
            left[row][p][i]=encoder.residue(a,p);
            right[row][p][i]=arithmetic.encode(encoder.residue(b,p));
        }
        plans[p].forward(left[row][p]);plans[p].forward(right[row][p]);
    }
    for(auto &component:result)component.fill({});
    // Row mapping matches CKKS tensor convention: cross, constant, square.
    constexpr std::array<std::array<int,2>,3> components{{{0,1},{1,1},{0,0}}};
    for(std::size_t diagonal=0;diagonal<2;++diagonal){
        std::array<Rows,3> sums;
#pragma omp parallel for collapse(2) schedule(static)
        for(std::size_t c=0;c<3;++c)for(std::size_t p=0;p<count;++p){
            auto &out=sums[c][p];out.assign(P::n,0);
            const auto &arithmetic=plans[p].arithmetic();
            const auto accumulate=[&](int a,int b){
                const auto *x=left[a][p].data(),*y=right[b][p].data();std::size_t i=0;
#if defined(__AVX2__)
                for(;i+8<=P::n;i+=8){
                    const auto value=arithmetic.accumulateVector(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(out.data()+i)),_mm256_loadu_si256(reinterpret_cast<const __m256i*>(x+i)),_mm256_loadu_si256(reinterpret_cast<const __m256i*>(y+i)));
                    _mm256_storeu_si256(reinterpret_cast<__m256i*>(out.data()+i),value);
                }
#endif
                for(;i<P::n;++i){const auto value=out[i]+arithmetic.multiplyLazy(x[i],y[i]);out[i]=value>=2*arithmetic.modulus()?value-2*arithmetic.modulus():value;}
            };
            for(std::size_t high=0;high<=diagonal;++high){
                const auto a=2*components[c][0]+high,b=2*components[c][1]+diagonal-high;
                accumulate(a,b);if(c==0)accumulate(2*components[c][1]+high,2*components[c][0]+diagonal-high);
            }
            for(auto &value:out)value=arithmetic.canonical(value);
            plans[p].inverse(out);
        }
#pragma omp parallel
        {
            std::vector<typename P::T> recovered(512);
#pragma omp for collapse(2) schedule(static)
            for(std::size_t c=0;c<3;++c)for(std::size_t begin=0;begin<P::n;begin+=512){
                const auto length=std::min(std::size_t{512},std::size_t(P::n)-begin);
                std::array<std::span<std::uint32_t>,count> residues;
                for(std::size_t p=0;p<count;++p)residues[p]=std::span(sums[c][p]).subspan(begin,length);
                crt.reconstruct(residues,std::span(recovered).first(length));
                for(std::size_t i=0;i<length;++i)result[c][begin+i]=ckks_detail::reduceToLevel<P,Q>(result[c][begin+i]+(recovered[i]<<(diagonal*split)));
            }
        }
    }
}
// Borrow an immutable existing relinearization key; only tensor arithmetic
// changes. Key switching and its auxiliary rounding remain unchanged.
template<class P,class Key>struct CKKSRelinKeyTensorRNS32 {
    static constexpr auto log_q=Key::log_q;
    const Key &key;
    void tensorProduct(TRLWE3<P>&out,const TRLWE<P>&a,const TRLWE<P>&b)const
    {CKKSTensorProductRNS32<P,log_q>(out,a,b);}
};
namespace ckks_detail {
template<class P,class Key>struct is_ckks_dd_relin_key<CKKSRelinKeyTensorRNS32<P,Key>> : std::true_type {};
}
template<class P,class Key>
void CKKSRelinKeySwitch(TRLWE<P>&out,const Polynomial<P>&in,const CKKSRelinKeyTensorRNS32<P,Key>&key)
{CKKSRelinKeySwitch<P>(out,in,key.key);}
} // namespace TFHEpp
#endif
