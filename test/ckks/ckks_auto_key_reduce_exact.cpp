#include "ckks_narrow_basis_benchmark.hpp"
#include <ckks/key_projection.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <iomanip>

using namespace TFHEpp;
using boost::multiprecision::cpp_int;
using P = narrow_benchmark::LeanBasisParam<9>;
using N = narrow_benchmark::Selector56BasisSchedule<narrow_benchmark::Secondary11BasisSchedule<
    CKKSNarrowBootstrapSchedule<P,48,518,47,39,31,true,44,true>>>;
using R = CKKSNarrowRefreshSchedule<N,44,37,0,true>;
using B = CKKSNarrowBasisRefreshSchedule<R,31>;

cpp_int integer(const P::T &v)
{
    cpp_int out=0;
    for (std::size_t j=v.limb.size();j-->0;) {out<<=64;out+=v.limb[j];}
    return out;
}
cpp_int residue(cpp_int x,const cpp_int &modulus)
{
    x%=modulus;if(x<0)x+=modulus;return x;
}

template <std::uint32_t Q,std::uint32_t SourceQ>
void check(const CKKSAutoKey<P,SourceQ> &source,const Key<P> &secret,
           const CKKSSlotVector<P> &slots,std::ostream &report,bool &first)
{
    auto reduced=std::make_unique<CKKSAutoKey<P,Q>>();
    CKKSReduceAutoKey<P,Q>(*reduced,source);
    const cpp_int modulus=cpp_int{1}<<Q,mask=modulus-1;
    std::array<int,P::n> s{};
    for(std::size_t i=0;i<P::n;++i) {
        if(secret[i]==P::T{})s[i]=0;
        else if(secret[i]==P::T{1})s[i]=1;
        else if(secret[i]==P::T{}-P::T{1})s[i]=-1;
        else throw std::runtime_error("retained secret is not ternary");
    }
    std::uint64_t coefficient_checks=0,phase_checks=0;
    cpp_int largest_noise=0;
    // Independent row selection by equality of encrypted gadget powers;
    // independent cpp_int negacyclic phase arithmetic, no production multiplier.
    constexpr auto width=std::numeric_limits<P::T>::digits;
    for(std::size_t r=0;r<reduced->key_switch_rows;++r) {
        const auto shift=width-(CKKSKeySwitchFirstRowForLevel<P,Q>()+r+1)*P::B̅gbit;
        std::size_t selected=source.key_switch_rows;
        for(std::size_t j=0;j<source.key_switch_rows;++j)
            if(width-(CKKSKeySwitchFirstRowForLevel<P,SourceQ>()+j+1)*P::B̅gbit==shift)selected=j;
        if(selected==source.key_switch_rows)throw std::runtime_error("gadget row absent");
        std::vector<cpp_int> phase(P::n),a(P::n);
        for(std::size_t c=0;c<=P::k;++c)
            for(std::size_t i=0;i<P::n;++i) {
                if(integer((*reduced)[0][r][c][i])!=(integer(source[0][selected][c][i])&mask))
                    throw std::runtime_error("public row reduction differs from integer mask");
                ++coefficient_checks;
            }
        for(std::size_t i=0;i<P::n;++i) {
            a[i]=integer((*reduced)[0][r][0][i]);
            phase[i]=integer((*reduced)[0][r][1][i]);
        }
        for(std::size_t i=0;i<P::n;++i)for(std::size_t j=0;j<P::n;++j)if(s[j]) {
            const auto index=i+j;
            if(index<P::n)phase[index]-=a[i]*s[j];
            else phase[index-P::n]+=a[i]*s[j];
        }
        for(std::size_t i=0;i<P::n;++i) {
            const int transformed=i==0?s[0]:-s[P::n-i];
            cpp_int error=residue(phase[i]-cpp_int{transformed}*(cpp_int{1}<<shift),modulus);
            if(error>=modulus/2)error-=modulus;
            if(error<0)error=-error;
            if(error>largest_noise)largest_noise=error;
            ++phase_checks;
        }
    }
    if(largest_noise>=64)throw std::runtime_error("projected encrypted gadget phase failed noise bound");
    using Cipher=CKKSCiphertext<P,Q,47>;
    auto input=std::make_unique<Cipher>(),raw=std::make_unique<Cipher>(),fft=std::make_unique<Cipher>();
    ckksSlotEncrypt<P>(*input,slots,secret,{P::α,0});
    CKKSEvalAuto<P,Q>(raw->ct,input->ct,2*P::n-1,*reduced);
    auto prepared=makeCKKSPreparedAutoKeyFFT<P,Q>(*reduced);
    CKKSEvalAuto<P,Q>(fft->ct,input->ct,2*P::n-1,*prepared);
    if(raw->ct!=fft->ct)throw std::runtime_error("raw/prepared reduced key switch differs");
    CKKSSlotVector<P> decoded{};ckksSlotDecrypt<P>(decoded,*fft,secret);
    double error=0;
    for(std::size_t i=0;i<P::n/2;++i)
        error=narrow_benchmark::finiteMaximum(error,std::abs(decoded[i]-std::conj(slots[i])));
    std::cout << "Q=" << Q << " scale=47 conjugation error=" << error << std::endl;
    if(error>=std::ldexp(1.0,-23))throw std::runtime_error("reduced key conjugation precision failed");
    if(!first)report<<',';first=false;
    report<<"{\"source_log_q\":"<<SourceQ<<",\"target_log_q\":"<<Q
          <<",\"ciphertext_log_delta\":47,\"coefficient_checks\":"<<coefficient_checks<<",\"independent_phase_checks\":"<<phase_checks
          <<",\"max_integer_gadget_noise\":"<<largest_noise<<",\"conjugation_max_error\":"<<error
          <<",\"raw_prepared_identical\":true}";
}

int main(int argc,char **argv)try {
    if(argc!=3)throw std::runtime_error("expected retained N512 checkpoint.bin and report path");
    auto dense=std::make_unique<Key<P>>(),sparse=std::make_unique<Key<P>>();
    auto slots=std::make_unique<CKKSSlotVector<P>>();auto input=std::make_unique<R::Input>();
    auto key=std::make_unique<CKKSNarrowRefreshKey<R>>();
    auto basis=std::make_unique<CKKSNarrowFrequencyBasisKey<B>>();
    auto retained=std::make_unique<CKKSCompactDDRelinKey<P,B::output_log_q>>();
    std::array<std::uint64_t,12> metadata{};
    {std::ifstream f(argv[1],std::ios::binary);if(!f)throw std::runtime_error("checkpoint missing");
     cereal::PortableBinaryInputArchive ar(f);ar(metadata);
     if(metadata!=std::array<std::uint64_t,12>{0x4e41525242433031ULL,P::n,576,518,47,31,85,44,18,8,44,37})
         throw std::runtime_error("incompatible retained checkpoint");
     ar(*dense,*sparse,*slots,*input,*key,*basis,*retained);}
    std::ofstream report(argv[2]);if(!report)throw std::runtime_error("cannot open report");
    report<<std::setprecision(17)<<"{\"n\":"<<P::n<<",\"fresh_keys_generated\":false,\"cases\":[";
    bool first=true;
    check<101>(basis->conjugation,*dense,*slots,report,first);
    check<100>(basis->conjugation,*dense,*slots,report,first);
    check<197>(basis->conjugation,*dense,*slots,report,first);
    check<345>(key->narrow.factor_galois[P::nbit],*dense,*slots,report,first);
    // All public Galois entries, including the unused zero entry, and aliased
    // equal-level destinations must retain exactly the same coefficients.
    auto galois=std::make_unique<CKKSGaloisKey<P,R::evaluation_input_log_q>>();
    CKKSReduceGaloisKey<P,R::evaluation_input_log_q>(*galois,key->s2c_input_galois);
    for(std::size_t i=0;i<galois->size();++i)
        if((*galois)[i].data!=key->s2c_input_galois[i].data)throw std::runtime_error("equal-level Galois copy differs");
    CKKSReduceGaloisKey<P,R::evaluation_input_log_q>(*galois,*galois);
    for(std::size_t i=0;i<galois->size();++i)
        if((*galois)[i].data!=key->s2c_input_galois[i].data)throw std::runtime_error("aliased Galois copy differs");
    report<<"],\"all_checks_passed\":true}";
    std::cout<<"Retained automorphism modulus reduction passed; no keys generated.\n";
} catch(const std::exception &e) {std::cerr<<e.what()<<'\n';return 1;}
