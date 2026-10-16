#pragma once

#include "scheme_switching.hpp"
#include <modular_ntt.hpp>

namespace TFHEpp {

// Functional packing switch: each public row encrypts a source-secret
// coefficient times a radix weight as a constant polynomial. Public digit
// polynomials place independent TLWE phases in consecutive coefficients.
// The destination is the bootstrap's actual structured sparse secret, never
// a zero-padded TFHE secret. The owner generates these rows once.
template <class P, class L>
struct TLWECoefficientPackingKey48 {
    static_assert(P::k == 1 && L::k == 1);
    static_assert(std::numeric_limits<typename L::T>::digits == 32);
    static constexpr std::uint64_t mask = (std::uint64_t{1} << 48) - 1;
    static constexpr std::size_t rows = L::n * 4;
    std::vector<std::uint64_t> data;
    static constexpr std::size_t index(std::size_t row, std::size_t c, std::size_t i)
    { return (row * 2 + c) * P::n + i; }
    template <class Archive> void serialize(Archive &ar) { ar(data); }
};

template <class P, class L>
void TLWECoefficientPackingKeyGen48(TLWECoefficientPackingKey48<P,L> &key,
    const Key<L> &source, const Key<P> &destination)
{
    using K = TLWECoefficientPackingKey48<P,L>;
    key.data.resize(K::rows * 2 * P::n);
#pragma omp parallel for schedule(static)
    for (std::size_t row=0; row<K::rows; ++row) {
        auto message=std::make_unique<Polynomial<P>>();
        auto ct=std::make_unique<TRLWE<P>>();
        const auto secret=ckks_scheme_switch_detail::ternarySecretToInt<L>(source[row/4]);
        (*message)[0]=P::T::from_signed_i64(secret) << (16+8*(row%4));
        ckks_detail::encryptPolynomialAtLevel<P,48>(*ct,*message,destination,
                                                  {std::ldexp(3.2,-48),0});
        for(std::size_t c=0;c<2;++c)for(std::size_t i=0;i<P::n;++i)
            key.data[K::index(row,c,i)]=(*ct)[c][i].limb[0]&K::mask;
    }
}

template <class P, class L>
class PreparedTLWECoefficientPacking48 {
    using K=TLWECoefficientPackingKey48<P,L>;
    std::array<std::unique_ptr<modular_ntt::NegacyclicNTTPlan>,2> plans_;
    std::vector<modular_ntt::ShoupMultiplier> spectra_;
    static constexpr std::size_t index(std::size_t p,std::size_t row,std::size_t c,std::size_t i)
    {return ((p*K::rows+row)*2+c)*P::n+i;}
public:
    explicit PreparedTLWECoefficientPacking48(const K &key)
      :spectra_(2*K::rows*2*P::n)
    {
        if(key.data.size()!=K::rows*2*P::n)throw std::invalid_argument("packing key length mismatch");
        for(std::size_t p=0;p<2;++p)plans_[p]=std::make_unique<modular_ntt::NegacyclicNTTPlan>(P::n,modular_ntt::degree65536_primes[p]);
#pragma omp parallel
        {
            std::vector<std::uint64_t> values(P::n);
#pragma omp for collapse(3) schedule(static)
            for(std::size_t p=0;p<2;++p)for(std::size_t row=0;row<K::rows;++row)for(std::size_t c=0;c<2;++c){
                for(std::size_t i=0;i<P::n;++i)values[i]=key.data[K::index(row,c,i)];
                plans_[p]->forward(values);
                for(std::size_t i=0;i<P::n;++i)spectra_[index(p,row,c,i)]={values[i],plans_[p]->modulus()};
            }
        }
    }
    std::size_t byteSize()const{return spectra_.size()*sizeof(spectra_[0]);}
    void apply(CKKSCiphertext<P,48,44> &out,const std::vector<TLWE<L>> &in)const
    {
        if(in.empty()||in.size()>P::n/2)throw std::invalid_argument("invalid TLWE packing batch");
        // Even summing all rows before CRT, |convolution| is below 2^84:
        // 48-bit rows * 8-bit digits * N * (4*L::n). Two 61-bit primes
        // provide ample exact centered reconstruction range.
        static_assert(P::nbit+std::bit_width(K::rows)+56<120);
        std::array<std::vector<std::uint64_t>,4> products;
        for(auto &v:products)v.resize(P::n);
        for(std::size_t p=0;p<2;++p){
            const auto modulus=plans_[p]->modulus();
#pragma omp parallel
            {
                std::vector<std::uint64_t> digits(P::n),sum0(P::n),sum1(P::n);
#pragma omp for schedule(static)
                for(std::size_t row=0;row<K::rows;++row){
                    std::fill(digits.begin(),digits.end(),0);
                    for(std::size_t j=0;j<in.size();++j)digits[j]=(in[j][row/4]>>(8*(row%4)))&255;
                    plans_[p]->forward(digits);
                    for(std::size_t i=0;i<P::n;++i){
                        sum0[i]=modular_ntt::add(sum0[i],spectra_[index(p,row,0,i)].apply(digits[i],modulus),modulus);
                        sum1[i]=modular_ntt::add(sum1[i],spectra_[index(p,row,1,i)].apply(digits[i],modulus),modulus);
                    }
                }
#pragma omp critical
                for(std::size_t i=0;i<P::n;++i){
                    products[2*p][i]=modular_ntt::add(products[2*p][i],sum0[i],modulus);
                    products[2*p+1][i]=modular_ntt::add(products[2*p+1][i],sum1[i],modulus);
                }
            }
            plans_[p]->inverse(products[2*p]);plans_[p]->inverse(products[2*p+1]);
        }
        const modular_ntt::TwoPrimeCRT crt(modular_ntt::degree65536_primes[0],modular_ntt::degree65536_primes[1]);
        for(std::size_t c=0;c<2;++c)for(std::size_t i=0;i<P::n;++i){
            const auto product=crt.reconstructSigned(products[c][i],products[2+c][i]);
            out.ct[c][i]=typename P::T{(std::uint64_t{0}-static_cast<std::uint64_t>(product))&K::mask};
        }
        for(std::size_t j=0;j<in.size();++j)
            out.ct[1][j]=ckks_detail::reduceToLevel<P,48>(out.ct[1][j]+(typename P::T{in[j][L::n]}<<16));
    }
};
} // namespace TFHEpp
