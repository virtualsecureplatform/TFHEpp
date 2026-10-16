#pragma once
#include "narrow_frequency_basis.hpp"

#ifdef TFHEPP_HAS_EXTENDED_MULTILIMB_PARAMS
namespace TFHEpp {
// A first multiplication layer derives only the roots consumed by the final
// order-two basis. No coefficient weighting or rounding crosses a product.
class CKKSNarrowTwoLayerFrequencyBasis {
public:
    struct Expression { std::size_t first,second;bool conjugate_second; };
private:
    std::vector<std::uint32_t> seeds_;
    CKKSNarrowFrequencyBasis final_;
    std::vector<Expression> expressions_;
public:
    CKKSNarrowTwoLayerFrequencyBasis(std::uint32_t degree,
        std::vector<std::uint32_t> seeds,std::vector<std::uint32_t> intermediate)
        :seeds_(std::move(seeds)),final_(degree,std::move(intermediate),true)
    {
        if(seeds_.empty() || !seeds_.front() || seeds_.back()>degree ||
           !std::is_sorted(seeds_.begin(),seeds_.end()) ||
           std::adjacent_find(seeds_.begin(),seeds_.end())!=seeds_.end())
            throw std::invalid_argument("invalid two-layer seed frequencies");
        for(auto frequency:final_.frequencies()){
            const auto direct=std::lower_bound(seeds_.begin(),seeds_.end(),frequency);
            if(direct!=seeds_.end() && *direct==frequency){
                expressions_.push_back({std::size_t(direct-seeds_.begin()),seeds_.size(),false});continue;
            }
            bool found=false;
            for(std::size_t a=0;a<seeds_.size()&&!found;++a)
                for(std::size_t b=a;b<seeds_.size();++b)if(seeds_[a]+seeds_[b]==frequency){
                    expressions_.push_back({a,b,false});found=true;break;
                }
            for(std::size_t a=0;a<seeds_.size()&&!found;++a)
                for(std::size_t b=0;b<a;++b)if(seeds_[a]-seeds_[b]==frequency){
                    expressions_.push_back({a,b,true});found=true;break;
                }
            if(!found)throw std::invalid_argument("two-layer intermediate root is not covered");
        }
    }
    const auto& seeds()const{return seeds_;}
    const auto& finalPlan()const{return final_;}
    const auto& expressions()const{return expressions_;}
};

template<class P,std::uint32_t Q,std::uint32_t Delta,class ProductKey,class Conjugate>
inline void CKKSNarrowDeriveFrequencyRoots(
    std::vector<CKKSCiphertext<P,Q-Delta,Delta>>&result,
    std::span<const CKKSCiphertext<P,Q,Delta>> seeds,
    const CKKSNarrowTwoLayerFrequencyBasis&plan,const ProductKey&key,Conjugate conjugate)
{
    static_assert(ProductKey::log_q==Q && Q>2*Delta);
    if(seeds.size()!=plan.seeds().size())throw std::invalid_argument("two-layer seed count mismatch");
    using Root=CKKSCiphertext<P,Q,Delta>;
    std::vector<std::unique_ptr<Root>> negative(seeds.size());
    result.resize(plan.expressions().size());
    for(std::size_t i=0;i<result.size();++i){
        const auto e=plan.expressions()[i];
        if(e.second==seeds.size()){CKKSLevelReduce<P>(result[i],seeds[e.first]);continue;}
        const auto*right=&seeds[e.second];
        if(e.conjugate_second){
            auto&cached=negative[e.second];
            if(!cached){cached=std::make_unique<Root>();conjugate(cached->ct,right->ct);}
            right=cached.get();
        }
        CKKSMultRelinThenRescale<P>(result[i],seeds[e.first],*right,key);
    }
}
} // namespace TFHEpp
#endif
