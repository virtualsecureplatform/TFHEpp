#include "ckks_narrow_basis_benchmark.hpp"
#include <ckks/prepared_compact_linear.hpp>
#include <iostream>
#include <random>

using namespace TFHEpp;
using P=narrow_benchmark::BenchmarkParam<9,7>;
using Plan=CKKSCompactRealLinearTransformPlan<P,93,52,37>;
using Prepared=CKKSPreparedCompactRealLinearTransformFFT<P,93,52,37>;

int main()
{
    Plan plan; plan.has_direct=true; plan.has_conjugate=true;
    for (auto *part:{&plan.direct,&plan.conjugate}) {
        part->k_step=2;
        for (int g=0;g<2;++g) {
            typename decltype(part->groups)::value_type group;
            group.giant_step=g;
            for (int j=0;j<2;++j) {
                typename decltype(group.terms)::value_type term;
                term.baby_step=j;
                for (std::size_t i=0;i<P::n;++i) {
                    const std::array<std::int64_t,5> values{
                        0,-1,std::numeric_limits<std::int64_t>::min(),
                        std::numeric_limits<std::int64_t>::max(),(std::int64_t{1}<<39)+3};
                    term.coefficients[i]=values[(i+j+g)%values.size()];
                }
                group.terms.push_back(term);
            }
            part->groups.push_back(std::move(group));
        }
    }
    bool rejected=false;
    try { Prepared bad(plan,0); } catch (const std::invalid_argument&) { rejected=true; }
    if (!rejected) throw std::runtime_error("insufficient spectrum budget accepted");
    const auto budget=std::uint64_t{8}*sizeof(Prepared::PlainFFT);
    auto input_key=std::make_unique<CKKSGaloisKey<P,93>>();
    auto output_key=std::make_unique<CKKSGaloisKey<P,56>>();
    CKKSRetainedPreparedGaloisKeyFFT<P,93> retained_input(*input_key);
    CKKSRetainedPreparedGaloisKeyFFT<P,56> retained_output(*output_key);
    const auto builds=retained_input.prepared_keys+retained_output.prepared_keys;
    auto input=std::make_unique<CKKSCiphertext<P,93,52>>();
    std::mt19937_64 random(20261015);
    for (auto &component:input->ct) for (auto &value:component) {
        for (auto &word:value.limb) word=random();
        value=ckks_detail::reduceToLevel<P,93>(value);
    }
    for (int mode=0;mode<5;++mode) {
        plan.has_direct=mode!=1; plan.has_conjugate=mode!=0;
        if (mode>=3)
            for (auto *part:{&plan.direct,&plan.conjugate})
                for (auto &group:part->groups) for (auto &term:group.terms)
                    for (std::size_t i=0;i<P::n;++i)
                        term.coefficients[i]=mode==4?0:(i%2?-17:23);
        Prepared prepared(plan,budget);
        if (mode>=3 && prepared.spectrum_bytes!=8*sizeof(PolynomialInFD<P>))
            throw std::runtime_error("zero high plaintext digits retained");
        auto reference=std::make_unique<CKKSPlainMulResult<P,93,52,37>>();
        CKKSLazyPreparedGaloisKeyFFT<P,93> transient_input(*input_key);
        CKKSLazyPreparedGaloisKeyFFT<P,56> transient_output(*output_key);
        CKKSRealLinearTransform<P,93,52,37>(*reference,*input,plan,transient_input,transient_output);
        for (int repeat=0;repeat<2;++repeat) {
            auto actual=std::make_unique<CKKSPlainMulResult<P,93,52,37>>();
            CKKSRealLinearTransform<P,93,52,37>(*actual,*input,prepared,retained_input,retained_output);
            if (actual->ct!=reference->ct) throw std::runtime_error("prepared compact transform differs from reference");
        }
    }
    if (retained_input.prepared_keys+retained_output.prepared_keys!=builds)
        throw std::runtime_error("public Galois spectra rebuilt online");
    std::cout<<"direct/conjugate/both transforms byte-identical; signed int64 extremes and omitted zero digits; repeated public cache reuse; budget guard passed\n";
}
