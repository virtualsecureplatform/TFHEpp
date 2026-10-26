#pragma once
#include "ckks_ram_fixture.hpp"
#include <ckks/subring_switch.hpp>
#include <iomanip>
#ifdef __linux__
#include <sys/resource.h>
#endif
#ifdef USE_BLAKE3
#include <blake3.h>
#endif

namespace ckks_ram {
using RingKey=TFHEpp::CKKSToTFHESubringSwitchKey32<P,L>;
using RingPrepared=TFHEpp::PreparedCKKSToTFHESubringSwitch32<P,L>;
using OutputInput=TFHEpp::CKKSCiphertext<P,32,29>;
inline constexpr const char *ring_file="to-tfhe32-subring-sigma256-v2.bin";
inline constexpr const char *direct256_file="to-tfhe32-direct-sigma256-v1.bin";
inline std::string fingerprint(const std::filesystem::path &path)
{
#ifdef USE_BLAKE3
    std::ifstream input(path,std::ios::binary);if(!input)throw std::runtime_error("cannot hash "+path.string());
    blake3_hasher hasher;blake3_hasher_init(&hasher);std::array<char,65536> buffer;
    while(input){input.read(buffer.data(),buffer.size());blake3_hasher_update(&hasher,buffer.data(),input.gcount());}
    if(!input.eof())throw std::runtime_error("hash read failed");
    std::array<std::uint8_t,32> digest;blake3_hasher_finalize(&hasher,digest.data(),digest.size());
    std::ostringstream result;result<<std::hex<<std::setfill('0');for(auto b:digest)result<<std::setw(2)<<unsigned(b);return result.str();
#else
    (void)path;throw std::runtime_error("RAM prototype fingerprints require USE_BLAKE3");
#endif
}
inline std::string optionalFingerprint(const std::filesystem::path &path){
#ifdef USE_BLAKE3
    return fingerprint(path);
#else
    (void)path;return {};
#endif
}
inline auto loadDenseCheckpoint(const std::filesystem::path &directory)
{
    auto dense=std::make_unique<TFHEpp::Key<P>>(),sparse=std::make_unique<TFHEpp::Key<P>>();
    std::ifstream file(directory/"checkpoint.bin",std::ios::binary);if(!file)throw std::runtime_error("missing checkpoint");
    cereal::PortableBinaryInputArchive ar(file);std::array<std::uint64_t,9> meta;std::string function;ar(meta,function);
    if(meta!=std::array<std::uint64_t,9>{0x4e415252504b3031ULL,P::n,85,41,379,49,34,100,448}||function!="tanh")
        throw std::runtime_error("RAM source checkpoint mismatch");
    ar(*dense,*sparse);return dense;
}
template<class KeyType> void readBridgeKey(const std::filesystem::path&path,KeyType&key,std::uint64_t kind,const std::string&source,const std::string&target)
{
    std::ifstream file(path,std::ios::binary);if(!file)throw std::runtime_error("missing public bridge key: "+path.string());
    TFHEpp::CKKSReadSubringBridgeKey<P,L>(file,key,kind,source,target);
}
template<class KeyType> void publishBridgeKey(const std::filesystem::path&path,const KeyType&key,std::uint64_t kind,const std::string&source,const std::string&target)
{
    const auto temporary=path.string()+".partial";
    const auto lock=path.string()+".lock";
    if(!std::filesystem::create_directory(lock))throw std::runtime_error("bridge publication already active");
    bool created=false;
    try {
        if(std::filesystem::exists(path)||std::filesystem::exists(temporary))throw std::runtime_error("bridge destination exists");
        {std::ofstream file(temporary,std::ios::binary);if(!file)throw std::runtime_error("cannot create public bridge archive");created=true;
         TFHEpp::CKKSWriteSubringBridgeKey<P,L>(file,key,kind,source,target);file.close();if(!file)throw std::runtime_error("bridge archive close failed");}
        // Publication is exclusive; the directory guard also protects the temporary file.
        std::filesystem::create_hard_link(temporary,path);std::filesystem::remove(temporary);std::filesystem::remove(lock);
    } catch(...) {std::error_code ignored;if(created)std::filesystem::remove(temporary,ignored);std::filesystem::remove(lock,ignored);throw;}

}
struct KeyNoiseAudit {
    std::size_t count=0;long double sum_squares=0;std::uint64_t l1=0,max_absolute=0;
    void add(std::uint32_t residue){const auto e=std::int64_t(std::bit_cast<std::int32_t>(residue));const auto a=std::uint64_t(e<0?-e:e);++count;sum_squares+=static_cast<long double>(e)*e;l1+=a;max_absolute=std::max(max_absolute,a);}
    double rms()const{return std::sqrt(double(sum_squares/count));}
    void check(double sigma)const{if(!count||rms()<.95*sigma||rms()>1.05*sigma)throw std::runtime_error("bridge key noise audit failed");}
    void json(std::ostream&out)const{out<<"{\"coefficients\":"<<count<<",\"rms\":"<<rms()<<",\"l1\":"<<l1<<",\"maximum_absolute\":"<<max_absolute<<'}';}
};
// Independent schoolbook small-ring phase audit. No production NTT/decomposition
// helper participates, and only aggregate residual statistics are published.
inline KeyNoiseAudit auditRingKey(const RingKey&key,const TFHEpp::Key<P>&source,const TFHEpp::Key<L>&target)
{
    KeyNoiseAudit audit;
    for(std::size_t row=0;row<RingKey::rows;++row)for(std::size_t r=0;r<RingKey::stride;++r)for(std::size_t i=0;i<L::n;++i){
        auto phase=key.ciphertexts[row][1][r+RingKey::stride*i];
        for(std::size_t j=0;j<L::n;++j){
            const auto index=i>=j?i-j:L::n+i-j;
            const auto term=key.ciphertexts[row][0][r+RingKey::stride*index]*target[j];
            if(i>=j)phase-=term;else phase+=term;
        }
        const auto expected=std::uint32_t(source[r+RingKey::stride*i].limb[0])<<(30-2*row);
        audit.add(phase-expected);
    }
    audit.check(256);
    if(std::ldexp(2.*double(audit.l1),-32)>.025746934137548232)
        throw std::runtime_error("ring key exceeds assessed uniform bridge-error bound");
    return audit;
}
inline KeyNoiseAudit auditDirectKey(const OutKey&key,const TFHEpp::Key<P>&source,const TFHEpp::Key<L>&target,double sigma)
{
    KeyNoiseAudit audit;
    // Distributed deterministic sample across the complete public key.
    constexpr std::size_t total=P::n*16*2,samples=4096;
    for(std::size_t sample=0;sample<samples;++sample){
        const auto index=sample*(total-1)/(samples-1),i=index/32,level=(index/2)%16,digit=index%2;
        auto phase=TFHEpp::tlweSymPhase<L>(key[i][level][digit],target);
        const auto expected=std::uint32_t(source[i].limb[0])*std::uint32_t(digit+1)<<(30-2*level);
        audit.add(phase-expected);
    }
    audit.check(sigma);return audit;
}
inline int addRingSwitchKeys(const std::filesystem::path &checkpoint,const std::filesystem::path &root)
{
    const auto source_id=fingerprint(checkpoint/"checkpoint.bin"),target_id=fingerprint(root/"owner.bin");
    auto dense=loadDenseCheckpoint(checkpoint);auto owner=std::make_unique<SavedOwner>();load(root/"owner.bin",*owner);
    auto ring=std::make_unique<RingKey>();const bool ring_reused=std::filesystem::exists(root/ring_file);
    if(ring_reused)readBridgeKey(root/ring_file,*ring,1,source_id,target_id);
    else TFHEpp::CKKSToTFHESubringSwitchKeyGen32(*ring,*dense,owner->secret());
    auto ring_audit=auditRingKey(*ring,*dense,owner->secret());
    if(!ring_reused)publishBridgeKey(root/ring_file,*ring,1,source_id,target_id);
    ring.reset();
    auto direct=std::make_unique<OutKey>();const bool direct_reused=std::filesystem::exists(root/direct256_file);
    if(direct_reused)readBridgeKey(root/direct256_file,*direct,2,source_id,target_id);
    else TFHEpp::CKKSToTFHEKeyGen<P,L,2,16>(*direct,*dense,owner->secret(),std::ldexp(256.,-32));
    auto direct_audit=auditDirectKey(*direct,*dense,owner->secret(),256);
    if(!direct_reused)publishBridgeKey(root/direct256_file,*direct,2,source_id,target_id);
    direct.reset();
    std::ofstream report(root/"ring-switch-key-audit.json");report<<std::setprecision(17)
        <<"{\"source_fingerprint\":\""<<source_id<<"\",\"target_fingerprint\":\""<<target_id
        <<"\",\"fingerprint_algorithm\":\"BLAKE3-256\",\"sigma\":256,\"fresh_private_keys_generated\":false,\"new_evaluation_keys_generated\":"<<(!ring_reused+!direct_reused)
        <<",\"ring_key_fingerprint\":\""<<fingerprint(root/ring_file)<<"\",\"direct256_key_fingerprint\":\""<<fingerprint(root/direct256_file)<<"\",\"ring\":";
    ring_audit.json(report);report<<",\"actual_uniform_ring_error_bound\":"<<std::ldexp(2.*double(ring_audit.l1),-32)<<",\"direct_sample\":";direct_audit.json(report);
    report<<",\"security_validated\":false}\n";if(!report)throw std::runtime_error("key audit write failed");
    std::cout<<"Ring and matched direct keys validated; new evaluation keys="<<(!ring_reused+!direct_reused)<<std::endl;return 0;
}
inline std::uint64_t peakRSSBytes(){
#ifdef __linux__
    rusage usage{};if(getrusage(RUSAGE_SELF,&usage)==0)return std::uint64_t(usage.ru_maxrss)*1024;
#endif
    return 0;
}
inline std::filesystem::path reportDirectory(const std::filesystem::path&root)
{
    const auto env=std::getenv("CKKS_RAM_REPORT_DIRECTORY");const auto result=env?std::filesystem::path(env):root;
    std::filesystem::create_directories(result);return result;
}
inline std::string outputSwitchMode()
{
    const auto env=std::getenv("CKKS_RAM_OUTPUT_SWITCH");std::string mode=env?env:"direct";
    if(mode!="direct"&&mode!="direct256"&&mode!="ring")throw std::invalid_argument("unknown RAM output switch");
    return mode;
}
class OutputBridge {
    std::unique_ptr<OutKey> direct_;
    std::unique_ptr<RingPrepared> ring_;
public:
    const std::string mode;
    std::string key_fingerprint;
    double preparation_seconds=0;
    std::size_t key_bytes=0,cache_bytes=0;
    TFHEpp::CKKSSubringSwitchTimings stages;
    OutputBridge(const std::filesystem::path&root,std::string selected,const std::string&source_id):mode(std::move(selected)){
        const auto started=std::chrono::steady_clock::now();
        const auto target_id=mode=="direct"?std::string{}:fingerprint(root/"owner.bin");
        const auto path=root/(mode=="ring"?ring_file:mode=="direct256"?direct256_file:"to-tfhe32.bin");
        if(mode=="ring"){
            auto raw=std::make_unique<RingKey>();readBridgeKey(path,*raw,1,source_id,target_id);key_bytes=sizeof(*raw);
            ring_=std::make_unique<RingPrepared>(*raw);cache_bytes=ring_->byteSize();
        }else if(mode=="direct"||mode=="direct256"){
            direct_=std::make_unique<OutKey>();if(mode=="direct")load(path,*direct_);else readBridgeKey(path,*direct_,2,source_id,target_id);key_bytes=sizeof(*direct_);
        }else throw std::invalid_argument("unknown output bridge");
        key_fingerprint=mode=="direct"?optionalFingerprint(path):fingerprint(path);preparation_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
    }
    void apply(std::vector<TFHEpp::TLWE<L>>&out,const OutputInput&input,std::size_t count){
        if(!count||count>P::n)throw std::invalid_argument("invalid output bridge count");
        if(ring_)ring_->apply(out,input,count,&stages);
        else {out.resize(count);
#pragma omp parallel for schedule(static)
            for(std::size_t i=0;i<count;++i)TFHEpp::CKKSCoeffToTFHE<P,L,32,29,2,16>(out[i],input,i,*direct_);
        }
    }
    void json(std::ostream&out)const{out<<"{\"mode\":\""<<mode<<"\",\"sigma\":"<<(mode=="direct"?std::ldexp(L::α,32):256)<<",\"key_fingerprint\":\""<<key_fingerprint<<"\",\"key_bytes\":"<<key_bytes<<",\"prepared_bytes\":"<<cache_bytes<<",\"setup_seconds\":"<<preparation_seconds<<",\"cumulative_kernel_seconds\":"<<stages.switching<<",\"cumulative_split_seconds\":"<<stages.splitting<<",\"cumulative_extraction_seconds\":"<<stages.extraction<<'}';}
};
inline void saveOutputInput(const std::filesystem::path&path,const OutputInput&input,const std::string&source_id,const std::string&target_id){
    if(std::filesystem::exists(path))throw std::runtime_error("saved bridge input exists");
    std::ofstream file(path,std::ios::binary);if(!file)throw std::runtime_error("cannot save bridge input");
    TFHEpp::CKKSWriteSubringBridgeKey<P,L>(file,input,3,source_id,target_id);
}
inline double auditInput(const OutputInput&input,const TFHEpp::Key<P>&source,std::size_t count=P::n/2){
    auto product=std::make_unique<TFHEpp::Polynomial<P>>();
    TFHEpp::ckks_detail::exactActiveTorusMulBySmall<P,32>(*product,input.ct[0],source);
    double error=0;
    for(std::size_t i=0;i<count;++i){
        const auto phase=std::uint32_t(input.ct[1][i].limb[0])-std::uint32_t((*product)[i].limb[0]);
        const double x=std::ldexp(double(std::bit_cast<std::int32_t>(phase)),-32);
        error=std::max(error,std::min(std::abs(x-.125),std::abs(x+.125)));
    }
    if(!std::isfinite(error)||error>1./16)throw std::runtime_error("pre-switch CKKS error exceeds assessed budget: "+std::to_string(error));
    return error;
}
} // namespace ckks_ram
