#pragma once
#include "prepared_gadget_rns32.hpp"
#include <span>
#include <type_traits>
#ifdef TFHEPP_HAS_EXTENDED_MULTILIMB_PARAMS
namespace TFHEpp {
namespace ckks_detail {
// Access immutable prepared spectra without changing the single-call kernel.
struct GadgetBatchRNS32Access {
    template<std::size_t Columns,class P,std::uint32_t Q,std::uint32_t Bits,std::uint32_t ChunkBits>
    static void apply(const CKKSPreparedGadgetRNS32<P,Q,Bits,ChunkBits>&key,
                      std::span<TRLWE<P>> output,std::span<const Polynomial<P>> input,std::size_t window)
    {
        using Cache=CKKSPreparedGadgetRNS32<P,Q,Bits,ChunkBits>;
        constexpr std::size_t rows=Cache::rows,count=Cache::count,chunks=Cache::chunks;
        constexpr auto chunk_bits=Cache::chunk_bits;
        using Small=typename Cache::Small;
        const auto&plans=Cache::plans();
        static_assert(Columns>0);
        if(output.size()!=input.size()||!window)
            throw std::invalid_argument("invalid gadget batch shape/window");
        using Digits=std::array<std::array<std::int64_t,P::n>,rows>;
        struct Work {
            std::unique_ptr<Digits> digits;
            std::array<bool,rows> active{};
            bool nonzero=false;
            std::array<std::vector<std::uint32_t>,count> transformed;
            std::array<std::array<std::vector<std::uint32_t>,count>,chunks*2> sums;
        };
        static const modular_rns32::SignedEncoder<count> encoder;
        static const modular_rns32::CenteredCRT<Small::limbs,count> crt;
        for(std::size_t offset=0;offset<input.size();offset+=window){
            const auto size=std::min(window,input.size()-offset);
            std::vector<Work> work(size);
            std::array<std::size_t,rows> active_rows{};
            std::size_t active_count=0;
#pragma omp parallel
            {
#pragma omp for schedule(static)
                for(std::size_t j=0;j<size;++j){
                    auto&w=work[j];w.digits=std::make_unique_for_overwrite<Digits>();
                    ckks_detail::activeBaseDecomposeSignedRows<P,Q,Bits,rows>(*w.digits,input[offset+j]);
                    for(std::size_t r=0;r<rows;++r){
                        w.active[r]=std::any_of((*w.digits)[r].begin(),(*w.digits)[r].end(),[](auto v){return v!=0;});
                        w.nonzero|=w.active[r];
                    }
                    for(auto&v:w.transformed)v.resize(rows*P::n);
                    for(auto&c:w.sums)for(auto&v:c)v.assign(P::n,0);
                }
#pragma omp single
                {
                    for(std::size_t r=0;r<rows;++r){
                        bool active=false;for(const auto&w:work)active|=w.active[r];
                        if(active)active_rows[active_count++]=r;
                    }
                }
#pragma omp for collapse(3) schedule(dynamic,1)
                for(std::size_t j=0;j<size;++j)for(std::size_t p=0;p<count;++p)for(std::size_t r=0;r<rows;++r){
                    auto&w=work[j];if(!w.active[r])continue;
                    auto values=std::span(w.transformed[p]).subspan(r*P::n,P::n);
                    encoder.residues(values,(*w.digits)[r],p,true);
                    plans[p].forwardBitReversed(values);
                }
#pragma omp for schedule(static)
                for(std::size_t j=0;j<size;++j)work[j].digits.reset();
#pragma omp for collapse(4) schedule(dynamic,1)
                for(std::size_t c=0;c<chunks*2;++c)for(std::size_t p=0;p<count;++p)
                    for(std::size_t begin=0;begin<P::n;begin+=256)for(std::size_t column=0;column<size;column+=Columns){
                        const auto end=std::min(begin+256,std::size_t(P::n));
                        const auto last=std::min(column+Columns,size);
                        const auto&arithmetic=plans[p].arithmetic();std::size_t first=0;
                        for(;first+4<=active_count;first+=4){
                            std::array<const std::uint32_t*,4> y;
                            std::array<std::array<const std::uint32_t*,4>,Columns> x;
                            for(std::size_t t=0;t<4;++t){
                                const auto r=active_rows[first+t];
                                y[t]=key.spectra_[p].data()+((c/2*rows+r)*2+c%2)*P::n;
                                for(auto j=column;j<last;++j)x[j-column][t]=work[j].transformed[p].data()+r*P::n;
                            }
                            std::size_t i=begin;
#if defined(__AVX2__)
                            for(;i+8<=end;i+=8){
                                const std::array<__m256i,4> key{
                                    _mm256_loadu_si256(reinterpret_cast<const __m256i*>(y[0]+i)),
                                    _mm256_loadu_si256(reinterpret_cast<const __m256i*>(y[1]+i)),
                                    _mm256_loadu_si256(reinterpret_cast<const __m256i*>(y[2]+i)),
                                    _mm256_loadu_si256(reinterpret_cast<const __m256i*>(y[3]+i))};
                                for(auto j=column;j<last;++j){
                                    if(!work[j].nonzero)continue;
                                    auto*out=work[j].sums[c][p].data();
                                    const auto value=arithmetic.accumulateDot4Vector(
                                        _mm256_loadu_si256(reinterpret_cast<const __m256i*>(out+i)),
                                        [&](std::size_t t){return _mm256_loadu_si256(reinterpret_cast<const __m256i*>(x[j-column][t]+i));},
                                        [&](std::size_t t){return key[t];});
                                    _mm256_storeu_si256(reinterpret_cast<__m256i*>(out+i),value);
                                }
                            }
#endif
                            for(;i<end;++i)for(auto j=column;j<last;++j){
                                std::uint64_t dot=0;for(std::size_t t=0;t<4;++t)dot+=std::uint64_t(x[j-column][t][i])*y[t][i];
                                auto&out=work[j].sums[c][p][i];const auto value=out+arithmetic.reduceLazy(dot);
                                out=value>=2*arithmetic.modulus()?value-2*arithmetic.modulus():value;
                            }
                        }
                        for(;first<active_count;++first){
                            const auto r=active_rows[first];const auto*y=key.spectra_[p].data()+((c/2*rows+r)*2+c%2)*P::n;
                            std::size_t i=begin;
#if defined(__AVX2__)
                            for(;i+8<=end;i+=8){
                                const auto key=_mm256_loadu_si256(reinterpret_cast<const __m256i*>(y+i));
                                for(auto j=column;j<last;++j){
                                    if(!work[j].active[r])continue;
                                    const auto*x=work[j].transformed[p].data()+r*P::n;auto*out=work[j].sums[c][p].data();
                                    const auto value=arithmetic.accumulateVector(_mm256_loadu_si256(reinterpret_cast<const __m256i*>(out+i)),_mm256_loadu_si256(reinterpret_cast<const __m256i*>(x+i)),key);
                                    _mm256_storeu_si256(reinterpret_cast<__m256i*>(out+i),value);
                                }
                            }
#endif
                            for(;i<end;++i)for(auto j=column;j<last;++j){
                                auto&out=work[j].sums[c][p][i];const auto x=work[j].transformed[p][r*P::n+i];
                                const auto value=out+arithmetic.multiplyLazy(x,y[i]);out=value>=2*arithmetic.modulus()?value-2*arithmetic.modulus():value;
                            }
                        }
                    }
#pragma omp for collapse(3) schedule(dynamic,1)
                for(std::size_t j=0;j<size;++j)for(std::size_t c=0;c<chunks*2;++c)for(std::size_t p=0;p<count;++p){
                    if(!work[j].nonzero)continue;
                    auto&v=work[j].sums[c][p];for(auto&x:v)x=plans[p].arithmetic().canonical(x);
                    plans[p].inverseBitReversed(v);
                }
                std::vector<Small> recovered(512);
#pragma omp for collapse(3) schedule(dynamic,1)
                for(std::size_t j=0;j<size;++j)for(std::size_t c=0;c<2;++c)for(std::size_t begin=0;begin<P::n;begin+=512){
                    const auto length=std::min(std::size_t{512},std::size_t(P::n)-begin);
                    for(std::size_t i=begin;i<begin+length;++i)output[offset+j][c][i]={};
                    if(!work[j].nonzero)continue;
                    for(std::size_t chunk=0;chunk<chunks;++chunk){
                        std::array<std::span<std::uint32_t>,count> residues;
                        for(std::size_t p=0;p<count;++p)residues[p]=std::span(work[j].sums[chunk*2+c][p]).subspan(begin,length);
                        crt.reconstruct(residues,std::span(recovered).first(length));
                        for(std::size_t i=0;i<length;++i){
                            typename P::T value{};
                            if constexpr(Small::limbs<P::T::limbs)if(recovered[i].limb.back()>>63)value.limb.fill(~std::uint64_t{0});
                            std::copy_n(recovered[i].limb.begin(),std::min(Small::limbs,P::T::limbs),value.limb.begin());
                            output[offset+j][c][begin+i]=ckks_detail::reduceToLevel<P,Q>(output[offset+j][c][begin+i]+(value<<(chunk*chunk_bits)));
                        }
                    }
                }
            }
        }
    }
};
} // namespace ckks_detail
// Global phases share key vectors across Columns outputs. Window bounds the
// temporary batch workspace; no input-dependent data is retained after return.
template<std::size_t Columns=4,class P,std::uint32_t Q,std::uint32_t Bits,std::uint32_t ChunkBits>
inline void CKKSApplyGadgetBatchRNS32(
    std::type_identity_t<std::span<TRLWE<P>>> output,
    std::type_identity_t<std::span<const Polynomial<P>>> input,
    const CKKSPreparedGadgetRNS32<P,Q,Bits,ChunkBits>&key,std::size_t window=8)
{
    ckks_detail::GadgetBatchRNS32Access::template apply<Columns,P,Q,Bits,ChunkBits>(key,output,input,window);
}
} // namespace TFHEpp
#endif
