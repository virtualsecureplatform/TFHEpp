#pragma once
#include "ckks.hpp"
#include <modular_rns32.hpp>

namespace TFHEpp::narrow_detail {
// Exact inverse transforms and CRT, followed immediately by the existing
// table-scale rounding. Each CRT block owns disjoint spans of every prime.
template<class Schedule,std::size_t Count,class First,class Plans,class CRT>
void inverseRoundFactorBatchRNS32(
    std::vector<First>&first,
    std::vector<TRLWE<typename Schedule::Param>>&second,
    const std::vector<std::array<std::vector<std::uint32_t>,Count>>&sums,
    std::size_t columns,std::size_t stride,const Plans&plans,const CRT&crt)
{
    using P=typename Schedule::Param;
    if(sums.size()!=2*(P::k+1)||stride<columns)
        throw std::invalid_argument("invalid factor finish shape");
    for(const auto&row:sums)for(const auto&v:row)
        if(v.size()!=stride*P::n)throw std::invalid_argument("invalid factor finish spectrum length");
    first.resize(columns);second.resize(columns);
    if(!columns)return;
    const auto rows=2*(P::k+1)*columns;
    std::vector<std::array<std::vector<std::uint32_t>,Count>> spectra(rows);
    for(auto&row:spectra)for(auto&v:row)v.resize(P::n);
#pragma omp parallel
    {
#pragma omp for collapse(2) schedule(dynamic,1)
        for(std::size_t row=0;row<rows;++row)for(std::size_t p=0;p<Count;++p){
            auto&v=spectra[row][p];const auto&sum=sums[row/columns][p];
            const auto column=row%columns;
            for(std::size_t i=0;i<P::n;++i)v[i]=sum[i*stride+column];
            plans[p].inverse(v);
        }
        std::vector<typename P::T> recovered(512);
#pragma omp for collapse(2) schedule(dynamic,1)
        for(std::size_t row=0;row<rows;++row)for(std::size_t begin=0;begin<P::n;begin+=512){
            const auto length=std::min(std::size_t{512},std::size_t(P::n)-begin);
            std::array<std::span<std::uint32_t>,Count> residues;
            for(std::size_t p=0;p<Count;++p)residues[p]=std::span(spectra[row][p]).subspan(begin,length);
            crt.reconstruct(residues,std::span(recovered).first(length));
            const auto t=row/columns/(P::k+1),component=row/columns%(P::k+1),column=row%columns;
            auto&first_ct=[&]() -> TRLWE<P>& {
                if constexpr(requires{first[column].ct;})return first[column].ct;
                else return first[column];
            }();
            auto&out=(t?second[column]:first_ct)[component];
            for(std::size_t i=0;i<length;++i)
                out[begin+i]=ckks_detail::roundedLevelRightShift<P,Schedule::selector_log_q,Schedule::plain_log_delta>(recovered[i]);
        }
    }
}
}
