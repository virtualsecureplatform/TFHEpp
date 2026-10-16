#include "ckks_narrow_basis_benchmark.hpp"
#include <ckks/prepared_compact_linear_rns32.hpp>
#include <iostream>
#include <random>
#include <boost/multiprecision/cpp_int.hpp>

using namespace TFHEpp;
using P=narrow_benchmark::BenchmarkParam<9,7>;
using Plan=CKKSCompactRealLinearTransformPlan<P,93,52,37>;
using Prepared=CKKSPreparedCompactRealLinearTransformRNS32<P,93,52,37>;

template<std::uint32_t Q>void checkWiderRotationArithmetic()
{
    using BP=narrow_benchmark::BenchmarkParam<9,7>;
    using Integer=boost::multiprecision::cpp_int;
    const auto integer=[](const BP::T&value){Integer out=0;for(std::size_t j=BP::T::limbs;j-->0;){out<<=64;out+=value.limb[j];}return out;};
    const Integer modulus=Integer{1}<<Q,mask=modulus-1,half=modulus/2;
    auto raw=std::make_unique<CKKSGaloisKey<BP,Q>>();std::mt19937_64 random(20261023+Q);
    for(std::size_t bit:{std::size_t{0},std::size_t{BP::nbit}})
        for(auto&source:(*raw)[bit].data)for(auto&row:source)for(auto&component:row)for(auto&v:component){for(auto&w:v.limb)w=random();v=ckks_detail::reduceToLevel<BP,Q>(v);}
    CKKSRetainedPreparedGaloisKeyFFT<BP,Q> original(*raw);
    CKKSRetainedPreparedGaloisKeyFFT<BP,Q,8> wider(*raw);
    if(wider.spectrum_bytes()*2!=original.spectrum_bytes())throw std::runtime_error("wider rotation cache did not halve spectra");
    for(int trial=0;trial<4;++trial){
        auto input=std::make_unique<TRLWE<BP>>();
        for(auto&component:*input){
            for(auto&v:component)for(auto&w:v.limb)w=random();
            component[0]={};component[1]=BP::T{1};component[2]=BP::T{}-BP::T{1};
            component[3]=BP::T{1}<<(Q-1);component[4]=(BP::T{1}<<Q)+BP::T{3};
        }
        // The FFT backend starts at N512. Keep eight nonzero mask
        // coefficients so a direct convolution oracle stays inexpensive;
        // higher storage bits still exercise noncanonical input handling.
        for(std::size_t i=8;i<BP::n;++i)(*input)[0][i]=((*input)[0][i]>>Q)<<Q;
        for(std::size_t bit:{std::size_t{0},std::size_t{BP::nbit}}){
            const auto exponent=bit?2*BP::n-1:5;
            auto actual=std::make_unique<TRLWE<BP>>(),permuted=std::make_unique<TRLWE<BP>>();
            CKKSEvalAuto<BP,Q>(*actual,*input,exponent,wider.get(bit));
            for(std::size_t c=0;c<2;++c)Automorphism<BP>((*permuted)[c],(*input)[c],exponent);
            std::array<std::vector<Integer>,2> expected;
            for(auto&component:expected)component.resize(BP::n);
            for(std::size_t i=0;i<BP::n;++i)expected[1][i]=integer((*permuted)[1][i])&mask;
            for(std::size_t i=0;i<BP::n;++i){
                Integer value=integer((*permuted)[0][i])&mask;if(value>=half)value-=modulus;
                for(std::size_t row=0;row<(Q+7)/8;++row){
                    Integer residue=(value+128)%256;if(residue<0)residue+=256;
                    const Integer digit=residue-128;value=(value-digit)/256;
                    if(digit==0)continue;
                    // Match the desired public power2^(8*row) to its original
                    // four-bit encrypted row; the oracle does not use FFTs.
                    const auto old_row=CKKSAutoKey<BP,Q>::key_switch_rows-1-2*row;
                    for(std::size_t c=0;c<2;++c)for(std::size_t j=0;j<BP::n;++j){
                        const auto destination=(i+j)%BP::n;
                        const Integer product=digit*integer((*raw)[bit][0][old_row][c][j]);
                        expected[c][destination]+=(i+j<BP::n?-product:product);
                    }
                }
            }
            for(std::size_t c=0;c<2;++c)for(std::size_t i=0;i<BP::n;++i)
                if(integer((*actual)[c][i])!=(expected[c][i]&mask))throw std::runtime_error("wider prepared rotation differs from integer convolution oracle");
        }
    }
    std::cout<<"wider S2C rotation Q="<<Q<<" primary8 matches independent integer oracle; half spectrum bytes\n";
}

int main()
{
    checkWiderRotationArithmetic<93>();checkWiderRotationArithmetic<56>();
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
    const auto budget=std::uint64_t{8}*Prepared::term_spectrum_bytes;
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
    {
        // Noncryptographic nonzero rows make FFT workspace corruption visible.
        auto fixture=std::make_unique<CKKSGaloisKey<P,93>>();
        for(std::size_t bit=0;bit<5;++bit)for(auto&source:(*fixture)[bit].data)for(auto&row:source)
            for(auto&component:row)for(auto&v:component){for(auto&w:v.limb)w=random();v=ckks_detail::reduceToLevel<P,93>(v);}
        CKKSRetainedPreparedGaloisKeyFFT<P,93> cache(*fixture);
        const auto prepared_count=cache.prepared_keys;
#ifdef _OPENMP
        const auto active_levels=omp_get_max_active_levels();omp_set_max_active_levels(2);
#endif
        for(int mode=0;mode<4;++mode){
            std::vector<bool> used(32,false);
            if(mode==0)std::fill(used.begin(),used.end(),true);
            if(mode==1)for(auto i:{3,7,13,17,31})used[i]=true;
            if(mode==2)used[1]=true;
            ckks_detail::CKKSSparseBabyStepTable<P> serial(32);
            ckks_detail::CKKSBuildSparseBabyStepRotationTable<P,93>(serial,input->ct,used,cache);
            for(std::size_t threads:{2,8}){
                ckks_detail::CKKSSparseBabyStepTable<P> concurrent(32);
                std::size_t peak=0;
                ckks_detail::CKKSBuildSparseBabyStepRotationTable<P,93>(concurrent,input->ct,used,cache,threads,&peak);
#ifdef _OPENMP
                if(mode==0 && !omp_get_dynamic() && omp_get_max_threads()>1 && omp_get_thread_limit()>1 && peak<2)throw std::runtime_error("baby table did not execute concurrent rotations");
#endif
                for(std::size_t i=0;i<serial.size();++i)
                    if(bool(serial[i])!=bool(concurrent[i]) || (serial[i]&&*serial[i]!=*concurrent[i]))throw std::runtime_error("concurrent baby table differs from serial nonzero-row fixture");
            }
        }
#ifdef _OPENMP
        if(omp_get_max_active_levels()!=2)throw std::runtime_error("baby table changed caller OpenMP policy");
        omp_set_max_active_levels(active_levels);
#endif
        CKKSRetainedPreparedGaloisKeyFFT<P,93,8> wider_cache(*fixture);
        std::vector<bool> wider_used(32,true);
        ckks_detail::CKKSSparseBabyStepTable<P> wider_serial(32),wider_parallel(32);
        ckks_detail::CKKSBuildSparseBabyStepRotationTable<P,93>(wider_serial,input->ct,wider_used,wider_cache);
        ckks_detail::CKKSBuildSparseBabyStepRotationTable<P,93>(wider_parallel,input->ct,wider_used,wider_cache,8);
        for(std::size_t i=0;i<wider_serial.size();++i)
            if(bool(wider_serial[i])!=bool(wider_parallel[i])||(wider_serial[i]&&*wider_serial[i]!=*wider_parallel[i]))
                throw std::runtime_error("wider retained rotations changed with parallel lookup");
        std::vector<bool> used(8,true);bool rejected=false;
        ckks_detail::CKKSSparseBabyStepTable<P> failed(8);
        CKKSLazyPreparedGaloisKeyFFT<P,93> lazy(*fixture);
        try{ckks_detail::CKKSBuildSparseBabyStepRotationTable<P,93>(failed,input->ct,used,lazy,8);}
        catch(const std::invalid_argument&){rejected=true;}
        if(!rejected || lazy.prepared_keys)throw std::runtime_error("parallel baby table accepted a mutable lazy cache");
        auto saved=std::move(cache.keys[1]->spectra);cache.keys[1]->spectra.clear();rejected=false;
        try{ckks_detail::CKKSBuildSparseBabyStepRotationTable<P,93>(failed,input->ct,used,cache,8);}
        catch(const std::invalid_argument&error){rejected=std::string(error.what())=="invalid prepared auto key";}
        cache.keys[1]->spectra=std::move(saved);
        if(!rejected || cache.prepared_keys!=prepared_count)throw std::runtime_error("parallel baby table failed to propagate worker error or reused keys incorrectly");
        std::cout<<"parallel baby tables match serial with nonzero rows; dependency gaps, worker errors and caller policy passed\n";
    }
    for (int mode=0;mode<6;++mode) {
        plan.has_direct=mode!=1; plan.has_conjugate=mode!=0;
        if (mode>=3)
            for (auto *part:{&plan.direct,&plan.conjugate})
                for (auto &group:part->groups) for (auto &term:group.terms)
                    for (std::size_t i=0;i<P::n;++i)
                        term.coefficients[i]=mode==4?0:(i%2?-17:23);
        if (mode==5) {
            for(auto *part:{&plan.direct,&plan.conjugate})
                for(auto &group:part->groups)for(auto &term:group.terms)term.coefficients.fill(0);
            plan.direct.groups[1].terms[1].coefficients[0]=42;
        }
        Prepared prepared(plan,budget);
        if(mode==5 && (!prepared.has_direct || prepared.has_conjugate ||
            prepared.direct.groups.size()!=1 || prepared.direct.groups[0].terms.size()!=1 ||
            prepared.direct.groups[0].giant_step!=1 || prepared.direct.groups[0].terms[0].baby_step!=1 ||
            prepared.spectrum_bytes!=Prepared::term_spectrum_bytes))
            throw std::runtime_error("zero groups or diagonals were not pruned correctly");
        if ((mode==3 && prepared.spectrum_bytes!=8*Prepared::term_spectrum_bytes) ||
            (mode==4 && (prepared.spectrum_bytes || prepared.has_direct || prepared.has_conjugate)))
            throw std::runtime_error("zero high plaintext digits retained");
        auto reference=std::make_unique<CKKSPlainMulResult<P,93,52,37>>();
        CKKSLazyPreparedGaloisKeyFFT<P,93> transient_input(*input_key);
        CKKSLazyPreparedGaloisKeyFFT<P,56> transient_output(*output_key);
        CKKSRealLinearTransform<P,93,52,37>(*reference,*input,plan,transient_input,transient_output);
        for (int repeat=0;repeat<2;++repeat) {
            auto actual=std::make_unique<CKKSPlainMulResult<P,93,52,37>>();
            for(auto &component:actual->ct)component.fill(P::T{12345});
            CKKSRealLinearTransform<P,93,52,37>(*actual,*input,prepared,retained_input,retained_output);
            if (actual->ct!=reference->ct) throw std::runtime_error("prepared compact transform differs from reference");
            prepared.apply(*actual,*input,retained_input,retained_output,nullptr,true,true);
            if (actual->ct!=reference->ct) throw std::runtime_error("Horner compact transform differs from reference");
        }
    }
    for(const std::vector<int> giants:{std::vector<int>{0,1,2,3,4,5,6,7},std::vector<int>{5,1,4,5,3},std::vector<int>{2,6,3}})
        for(int mode=0;mode<3;++mode){
            Plan varied;varied.has_direct=mode!=1;varied.has_conjugate=mode!=0;
            for(auto*part:{&varied.direct,&varied.conjugate}){
                part->k_step=part==&varied.direct?2:3;
                for(const auto giant:giants){
                    decltype(part->groups)::value_type group;group.giant_step=giant;
                    for(int j=0;j<3;++j){
                        decltype(group.terms)::value_type term;term.baby_step=j%2;
                        for(std::size_t i=0;i<P::n;++i)term.coefficients[i]=((i+j+giant)%2?-1:1)*std::int64_t((std::uint64_t{1}<<39)+i+j+1);
                        group.terms.push_back(term);
                    }
                    part->groups.push_back(std::move(group));
                }
            }
            Prepared prepared(varied,6*giants.size()*Prepared::term_spectrum_bytes);
            auto reference=std::make_unique<CKKSPlainMulResult<P,93,52,37>>(),ordinary=std::make_unique<CKKSPlainMulResult<P,93,52,37>>(),horner=std::make_unique<CKKSPlainMulResult<P,93,52,37>>();
            CKKSLazyPreparedGaloisKeyFFT<P,93> transient_input(*input_key);
            CKKSLazyPreparedGaloisKeyFFT<P,56> transient_output(*output_key);
            CKKSRealLinearTransform<P,93,52,37>(*reference,*input,varied,transient_input,transient_output);
            CKKSLinearRNSTimings before,after;
            prepared.apply(*ordinary,*input,retained_input,retained_output,&before,true,false);
            prepared.apply(*horner,*input,retained_input,retained_output,&after,true,true);
            if(horner->ct!=ordinary->ct || horner->ct!=reference->ct)throw std::runtime_error("dense, sparse or repeated giant steps changed the Horner result");
            for(bool use_horner:{false,true})for(std::size_t workers:{2,8}){
                CKKSLinearRNSTimings concurrent;
                auto parallel=std::make_unique<CKKSPlainMulResult<P,93,52,37>>();
#ifdef _OPENMP
                const auto old_levels=omp_get_max_active_levels();omp_set_max_active_levels(2);
#endif
                prepared.apply(*parallel,*input,retained_input,retained_output,&concurrent,true,use_horner,0,workers);
                if(parallel->ct!=reference->ct)throw std::runtime_error("parallel giant ranges changed exact zero-key result");
#ifdef _OPENMP
                if(omp_get_max_active_levels()!=2)throw std::runtime_error("giant workers changed caller OpenMP policy");
                omp_set_max_active_levels(old_levels);
                if(omp_get_max_threads()>1 && !omp_get_dynamic() && omp_get_thread_limit()>1 && concurrent.group_peak_workers<2)
                    throw std::runtime_error("giant ranges did not execute concurrently");
#endif
            }
            if(mode==0 && giants.size()==8){
                auto fixture=std::make_unique<CKKSGaloisKey<P,56>>();
                for(std::size_t bit=0;bit<5;++bit)for(auto&source:(*fixture)[bit].data)for(auto&row:source)
                    for(auto&component:row)for(auto&v:component){for(auto&w:v.limb)w=random();v=ckks_detail::reduceToLevel<P,56>(v);}
                CKKSRetainedPreparedGaloisKeyFFT<P,56> nonzero(*fixture);
                auto serial=std::make_unique<CKKSPlainMulResult<P,93,52,37>>(),parallel=std::make_unique<CKKSPlainMulResult<P,93,52,37>>();
                prepared.apply(*serial,*input,retained_input,nonzero,nullptr,true,false);
                for(std::size_t workers:{2,8}){
                    prepared.apply(*parallel,*input,retained_input,nonzero,nullptr,true,false,0,workers);
                    if(parallel->ct!=serial->ct)throw std::runtime_error("parallel independent giant rotations differ with nonzero key rows");
                }
#ifdef _OPENMP
                if(omp_get_max_threads()>1){
                    bool rejected=false;
                    CKKSLazyPreparedGaloisKeyFFT<P,56> lazy(*output_key);
                    try{prepared.apply(*horner,*input,retained_input,lazy,nullptr,true,true,0,8);}
                    catch(const std::invalid_argument&){rejected=true;}
                    if(!rejected || lazy.prepared_keys)throw std::runtime_error("giant workers accepted a lazy mutable output cache");
                    auto saved=std::move(retained_output.keys[1]->spectra);retained_output.keys[1]->spectra.clear();rejected=false;
                    try{prepared.apply(*horner,*input,retained_input,retained_output,nullptr,true,true,0,8);}
                    catch(const std::invalid_argument&e){rejected=std::string(e.what())=="invalid prepared auto key";}
                    retained_output.keys[1]->spectra=std::move(saved);
                    if(!rejected)throw std::runtime_error("giant worker error did not propagate");
                }
#endif
            }
            if(mode==0 && giants.size()==8 && (before.giant_rotation_bits!=12 || after.giant_rotation_bits!=7))throw std::runtime_error("dense Horner binary rotation counts differ");
        }
    for(std::size_t terms:{3,4,5,7,8,9}){
        auto varied=plan;varied.has_direct=true;varied.has_conjugate=true;
        for(auto*part:{&varied.direct,&varied.conjugate})for(auto&group:part->groups){
            group.terms.clear();
            for(std::size_t j=0;j<terms;++j){
                decltype(group.terms)::value_type term;term.baby_step=j%2;
                for(std::size_t i=0;i<P::n;++i)term.coefficients[i]=((i+j)%2?-1:1)*std::int64_t((std::uint64_t{1}<<39)+j+1);
                group.terms.push_back(term);
            }
        }
        Prepared prepared(varied,4*terms*Prepared::term_spectrum_bytes);
        auto reference=std::make_unique<CKKSPlainMulResult<P,93,52,37>>(),fused=std::make_unique<CKKSPlainMulResult<P,93,52,37>>(),unfused=std::make_unique<CKKSPlainMulResult<P,93,52,37>>();
        CKKSLazyPreparedGaloisKeyFFT<P,93> transient_input(*input_key);
        CKKSLazyPreparedGaloisKeyFFT<P,56> transient_output(*output_key);
        CKKSRealLinearTransform<P,93,52,37>(*reference,*input,varied,transient_input,transient_output);
        prepared.apply(*fused,*input,retained_input,retained_output,nullptr,true);
        prepared.apply(*unfused,*input,retained_input,retained_output,nullptr,false);
        if(fused->ct!=unfused->ct || fused->ct!=reference->ct)throw std::runtime_error("four-term S2C matrix batch or remainder changed exact output");
    }
    bool bound_rejected=false;
    try { CKKSPreparedCompactRealLinearTransformRNS32<P,93,52,37,2> insufficient(plan,budget); }
    catch(const std::invalid_argument&){bound_rejected=true;}
    if(!bound_rejected)throw std::runtime_error("insufficient RNS product accepted");
    // More terms than the N512 fixed FFT batch: stress signed digit maxima
    // and verify the complete transform against a uint128 modular oracle.
    Plan batched;batched.has_direct=true;batched.has_conjugate=false;
    batched.direct.k_step=1;
    decltype(batched.direct.groups)::value_type group;group.giant_step=0;
    constexpr std::uint64_t coefficient=(std::uint64_t{1}<<22)-1;
    for(std::size_t j=0;j<512;++j){
        decltype(group.terms)::value_type term;term.baby_step=0;term.coefficients.fill(coefficient);
        group.terms.push_back(std::move(term));
    }
    batched.direct.groups.push_back(std::move(group));
    P::T value{};
    for(std::size_t j=0;j<6;++j)value+=P::T::from_signed_i64(-32768)<<(16*j);
    value=ckks_detail::reduceToLevel<P,93>(value);
    for(auto &component:input->ct)component.fill(value);
    const unsigned __int128 mask93=(static_cast<unsigned __int128>(1)<<93)-1;
    const auto raw=(static_cast<unsigned __int128>(value.limb[1])<<64)|value.limb[0];
    for(bool adaptive:{true}){
        Prepared prepared(batched,std::uint64_t{512}*Prepared::term_spectrum_bytes);
        auto actual=std::make_unique<CKKSPlainMulResult<P,93,52,37>>();
        CKKSRealLinearTransform<P,93,52,37>(*actual,*input,prepared,retained_input,retained_output);
        bool exact=true;
        for(std::size_t i=0;i<P::n;++i){
            const auto multiplicity=std::int64_t(2*i+2)-P::n;
            const auto factor=static_cast<unsigned __int128>(static_cast<__int128>(multiplicity)*512*coefficient);
            const auto product=(raw*factor)&mask93;
            const bool negative=(product>>(93-1))!=0;
            const auto magnitude=negative?(static_cast<unsigned __int128>(1)<<93)-product:product;
            const auto rounded_magnitude=(magnitude+(static_cast<unsigned __int128>(1)<<36))>>37;
            const auto rounded=std::uint64_t(negative?static_cast<unsigned __int128>(0)-rounded_magnitude:rounded_magnitude)&((std::uint64_t{1}<<56)-1);
            for(const auto &component:actual->ct)
                if(component[i]!=P::T{rounded})exact=false;
        }
        if(adaptive && !exact)throw std::runtime_error("exact RNS batch differs from independent constant-polynomial convolution");
        std::cout<<"constant-polynomial stress adaptive="<<adaptive<<" exact="<<exact<<'\n';
    }
    if (retained_input.prepared_keys+retained_output.prepared_keys!=builds)
        throw std::runtime_error("public Galois spectra rebuilt online");
    std::cout<<"exact RNS direct/conjugate/both transforms match FFT; signed int64, zero operators, budget/bound guards and integer stress oracle passed\n";
}
