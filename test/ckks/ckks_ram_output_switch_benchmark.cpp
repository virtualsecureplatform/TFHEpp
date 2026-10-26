#include "ckks_ram_output_switch.hpp"
#include <sys/resource.h>

int main(int argc,char**argv)try {
    using namespace ckks_ram;using namespace TFHEpp;
    if(argc!=5){std::cerr<<"Usage: ckks_ram_output_switch_benchmark CHECKPOINT FIXTURE INPUT32 OUTPUT_DIRECTORY\n";return 2;}
    const std::filesystem::path checkpoint=argv[1],fixture=argv[2],input_path=argv[3],output=argv[4];
    std::filesystem::create_directories(output);
    if(std::filesystem::exists(output/"bridge.jsonl")||std::filesystem::exists(output/"bridge.json"))throw std::runtime_error("benchmark receipts exist");
    const auto source_id=fingerprint(checkpoint/"checkpoint.bin"),target_id=fingerprint(fixture/"owner.bin");
    auto input=std::make_unique<OutputInput>();{std::ifstream file(input_path,std::ios::binary);if(!file)throw std::runtime_error("missing derived input");CKKSReadSubringBridgeKey<P,L>(file,*input,3,source_id,target_id);}
    auto dense=loadDenseCheckpoint(checkpoint);auto owner=std::make_unique<SavedOwner>();load(fixture/"owner.bin",*owner);
    auto product=std::make_unique<Polynomial<P>>();ckks_detail::exactActiveTorusMulBySmall<P,32>(*product,input->ct[0],*dense);
    std::vector<std::uint32_t> phases(P::n/2);
    for(std::size_t i=0;i<phases.size();++i)phases[i]=std::uint32_t(input->ct[1][i].limb[0])-std::uint32_t((*product)[i].limb[0]);
    const auto pre_error=auditInput(*input,*dense);
    std::array<std::unique_ptr<OutputBridge>,3> bridges;
    const std::array<std::string,3> names{"direct","direct256","ring"};
    for(std::size_t b=0;b<3;++b)bridges[b]=std::make_unique<OutputBridge>(fixture,names[b],source_id);
    const std::array<std::array<int,3>,6> permutations{{{0,1,2},{2,1,0},{1,0,2},{2,0,1},{0,2,1},{1,2,0}}};
    std::ofstream journal(output/"bridge.jsonl");journal<<std::setprecision(17);
    for(auto count:{256u,512u,1024u,8192u}){
        for(int trial=-1;trial<6;++trial)for(int order=0;order<3;++order){
            const auto b=trial<0?order:permutations[trial][order];auto &bridge=*bridges[b];
            auto previous=bridge.stages;std::vector<TLWE<L>> out;
            const auto before=std::chrono::steady_clock::now();bridge.apply(out,*input,count);
            const auto seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-before).count();
            double max_error=0;
            for(std::size_t i=0;i<count;++i){
                const auto actual=tlweSymPhase<L>(out[i],owner->secret());
                const auto expected=std::bit_cast<std::int32_t>(phases[i]);
                if(tlweSymDecrypt<L>(out[i],owner->secret())!=(expected>0))throw std::runtime_error("bridge benchmark bit mismatch");
                max_error=std::max(max_error,std::abs(std::ldexp(double(std::bit_cast<std::int32_t>(actual-phases[i])),-32)));
            }
            if(b==2&&max_error>.025746934137548232)throw std::runtime_error("bridge exceeds assessed error bound");
            journal<<"{\"bits\":"<<count<<",\"backend\":\""<<names[b]<<"\",\"trial\":"<<trial<<",\"order\":"<<order<<",\"warmup\":"<<(trial<0?"true":"false")<<",\"seconds\":"<<seconds
                <<",\"kernel_seconds\":"<<bridge.stages.switching-previous.switching<<",\"split_seconds\":"<<bridge.stages.splitting-previous.splitting<<",\"extraction_seconds\":"<<bridge.stages.extraction-previous.extraction
                <<",\"maximum_added_phase_error\":"<<max_error<<",\"correct\":true}\n";journal.flush();if(!journal)throw std::runtime_error("benchmark journal write failed");
            std::cout<<"bridge bits="<<count<<" trial="<<trial<<" backend="<<names[b]<<" seconds="<<seconds<<" passed"<<std::endl;
        }
    }
    rusage usage{};getrusage(RUSAGE_SELF,&usage);
    std::ofstream report(output/"bridge.json");report<<std::setprecision(17)<<"{\"source_fingerprint\":\""<<source_id<<"\",\"target_fingerprint\":\""<<target_id<<"\",\"input_fingerprint\":\""<<fingerprint(input_path)
        <<"\",\"pre_switch_error\":"<<pre_error<<",\"combined_peak_rss_bytes\":"<<std::uint64_t(usage.ru_maxrss)*1024<<",\"all_checks_passed\":true,\"security_validated\":false,\"backends\":[";
    for(std::size_t b=0;b<3;++b){if(b)report<<',';bridges[b]->json(report);}report<<"]}\n";
    if(!report)throw std::runtime_error("benchmark summary write failed");
    return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}
