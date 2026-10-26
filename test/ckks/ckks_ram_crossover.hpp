#pragma once
#include "ckks_ram_output_switch.hpp"

namespace ckks_ram {
// Pair identical encrypted requests, alternating execution order. Each RAM
// retains its own state across requests; correctness checks are outside timing.
template<std::size_t A,std::size_t W,class TFHERefresh,class CKKSRefresh,class Validate,class Describe>
int runPaired(const std::filesystem::path &root,const Fixture<A,W> &fixture,
              const TFHEpp::Key<L> &secret,const TFHEpp::EvalKey &ek,
              TFHERefresh tfhe_refresh,CKKSRefresh ckks_refresh,
              double &packing,double &bootstrap,double &decoding,double &switching,
              Validate validate,Describe describe)
{
    using RAM=typename Fixture<A,W>::RAM;
    constexpr std::size_t count=std::size_t{1}<<A;
    const auto requests=trace<A,W>();
    if(requests.size()<5||fixture.addresses.size()!=requests.size()||fixture.enables.size()!=requests.size()||fixture.data.size()!=requests.size())
        throw std::runtime_error("paired RAM fixture length mismatch");
    std::array<std::unique_ptr<RAM>,2> memories;
    std::array<std::array<std::uint32_t,count>,2> expected{};
    for(std::size_t backend=0;backend<2;++backend){
        memories[backend]=std::make_unique<RAM>();memories[backend]->memory=fixture.initial;
        for(std::size_t i=0;i<count;++i)expected[backend][i]=i&((1u<<W)-1);
    }
    const auto prefix="crossover-"+std::to_string(count)+"x"+std::to_string(W);
    std::ofstream records(root/(prefix+".jsonl"));
    if(!records)throw std::runtime_error("cannot create paired RAM receipt");
    records<<std::setprecision(17);
    std::array<std::vector<double>,2> times,refresh_times;
    std::array<double,2> phase_error{};
    double maximum_input_error=0;
    std::array<double,4> stage_start{};
    for(std::size_t step=0;step<5;++step){
        if(step==1)stage_start={packing,bootstrap,decoding,switching};
        for(std::size_t order=0;order<2;++order){
            const std::size_t backend=order^(step&1);
            typename RAM::Timing timing;
            const auto refresh=[&](const auto &in){return backend?ckks_refresh(in):tfhe_refresh(in);};
            const auto before=std::chrono::steady_clock::now();
            const auto output=memories[backend]->access(fixture.addresses[step],fixture.enables[step],fixture.data[step],ek,refresh,&timing);
            const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-before).count();
            const double input_error=backend?validate(count*W):0;
            maximum_input_error=std::max(maximum_input_error,input_error);
            for(std::size_t w=0;w<W;++w)
                if(TFHEpp::tlweSymDecrypt<L>(output[w],secret)!=bool((expected[backend][requests[step].address]>>w)&1))
                    throw std::runtime_error("paired RAM read mismatch");
            if(requests[step].enable)expected[backend][requests[step].address]=requests[step].data;
            double error=0;
            for(std::size_t w=0;w<W;++w)for(std::size_t i=0;i<count;++i){
                TFHEpp::TLWE<L> cell;TFHEpp::SampleExtractIndex<L>(cell,memories[backend]->memory[w][i],0);
                const bool bit=(expected[backend][i]>>w)&1;
                if(TFHEpp::tlweSymDecrypt<L>(cell,secret)!=bit)
                    throw std::runtime_error("paired RAM state mismatch at cell "+std::to_string(i));
                const double phase=std::ldexp(double(std::bit_cast<std::int32_t>(TFHEpp::tlweSymPhase<L>(cell,secret))),-32);
                error=std::max(error,std::abs(phase-(bit?0.125:-0.125)));
            }
            phase_error[backend]=std::max(phase_error[backend],error);
            if(step){times[backend].push_back(seconds);refresh_times[backend].push_back(timing.refresh);}
            records<<"{\"backend\":\""<<(backend?"ckks":"tfhe-serial")<<"\",\"bits\":"<<count*W
                   <<",\"step\":"<<step<<",\"order\":"<<order<<",\"warmup\":"<<(step?"false":"true")
                   <<",\"seconds\":"<<seconds<<",\"refresh_seconds\":"<<timing.refresh
                   <<",\"address_seconds\":"<<timing.address<<",\"read_seconds\":"<<timing.read<<",\"write_seconds\":"<<timing.write
                   <<",\"max_phase_error\":"<<error<<",\"pre_switch_error\":"<<input_error<<",\"correct\":true}\n";records.flush();
            std::cout<<"paired RAM bits="<<count*W<<" step="<<step<<" backend="<<(backend?"CKKS":"TFHE")
                     <<" seconds="<<seconds<<" all cells passed"<<std::endl;
        }
    }
    const auto values=[](std::ostream &out,const std::vector<double> &v){out<<'[';for(std::size_t i=0;i<v.size();++i){if(i)out<<',';out<<v[i];}out<<']';};
    std::ofstream summary(root/(prefix+".json"));summary<<std::setprecision(17);
    summary<<"{\"address_bits\":"<<A<<",\"word_bits\":"<<W<<",\"warmups_per_backend\":1,\"measured_accesses_per_backend\":4,\"alternating_order\":true,\"all_cells_checked\":true,\"backends\":[";
    for(std::size_t backend=0;backend<2;++backend){
        if(backend)summary<<',';
        summary<<"{\"name\":\""<<(backend?"ckks":"tfhe-serial")<<"\",\"access_seconds\":";
        values(summary,times[backend]);summary<<",\"refresh_seconds\":";values(summary,refresh_times[backend]);
        summary<<",\"max_phase_error\":"<<phase_error[backend]<<'}';
    }
    summary<<"],\"measured_ckks_stages\":{\"packing\":"<<packing-stage_start[0]<<",\"bootstrap\":"<<bootstrap-stage_start[1]
           <<",\"decoding\":"<<decoding-stage_start[2]<<",\"switching\":"<<switching-stage_start[3]<<"},\"lvl1_sigma\":"<<std::ldexp(L::α,32)
           <<",\"pre_switch_max_error\":"<<maximum_input_error<<",\"security_validated\":false,\"output_bridge\":";
    describe(summary);summary<<"}\n";
    if(!summary)throw std::runtime_error("paired RAM summary write failed");
    return 0;
}
} // namespace ckks_ram

namespace ckks_ram {
struct BridgeStageTimes {
    double packing=0,bootstrap=0,decoding=0,switching=0;
    void json(std::ostream&out)const{out<<"{\"packing\":"<<packing<<",\"bootstrap\":"<<bootstrap<<",\"decoding\":"<<decoding<<",\"switching\":"<<switching<<'}';}
};
// Four independent RAM states; each measured backend occupies every execution
// position once. Validation, source-phase audit and receipt I/O are untimed.
template<std::size_t A,class Refresh,class Validate,class Describe>
int runRingComparison(const std::filesystem::path &directory,const Fixture<A,1>&fixture,
    const TFHEpp::Key<L>&secret,const TFHEpp::EvalKey&ek,
    Refresh refresh,Validate validate,Describe describe,std::array<BridgeStageTimes,3>&stages)
{
    using RAM=typename Fixture<A,1>::RAM;constexpr std::size_t count=std::size_t{1}<<A;
    const auto requests=trace<A,1>();
    if(requests.size()<5 || fixture.addresses.size()!=requests.size() || fixture.enables.size()!=requests.size() || fixture.data.size()!=requests.size())
        throw std::runtime_error("ring comparison fixture length mismatch");
    const std::array<std::string,4> names{"tfhe-serial","direct","direct256","ring"};
    std::array<std::unique_ptr<RAM>,4> memories;
    std::array<std::array<std::uint32_t,count>,4> expected{};
    for(std::size_t b=0;b<4;++b){memories[b]=std::make_unique<RAM>();memories[b]->memory=fixture.initial;for(std::size_t i=0;i<count;++i)expected[b][i]=i&1;}
    const auto prefix=directory/("ring-comparison-"+std::to_string(count)+"x1");
    if(std::filesystem::exists(prefix.string()+".jsonl")||std::filesystem::exists(prefix.string()+".json"))throw std::runtime_error("ring comparison receipt exists");
    std::ofstream journal(prefix.string()+".jsonl");journal<<std::setprecision(17);
    std::array<std::vector<double>,4> times;
    std::array<BridgeStageTimes,3> start{};
    for(std::size_t step=0;step<5;++step){
        if(step==1)start=stages;
        for(std::size_t order=0;order<4;++order){
            const auto b=(order+step)%4;
            typename RAM::Timing timing;
            const auto before_stages=b?stages[b-1]:BridgeStageTimes{};
            const auto before=std::chrono::steady_clock::now();
            const auto output=memories[b]->access(fixture.addresses[step],fixture.enables[step],fixture.data[step],ek,
                [&](const auto&input){return refresh(b,input);},&timing);
            const auto elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-before).count();
            if(TFHEpp::tlweSymDecrypt<L>(output[0],secret)!=bool(expected[b][requests[step].address]))throw std::runtime_error("ring RAM read mismatch");
            if(requests[step].enable)expected[b][requests[step].address]=requests[step].data;
            double error=0;
            for(std::size_t i=0;i<count;++i){
                TFHEpp::TLWE<L> cell;TFHEpp::SampleExtractIndex<L>(cell,memories[b]->memory[0][i],0);
                if(TFHEpp::tlweSymDecrypt<L>(cell,secret)!=bool(expected[b][i]))throw std::runtime_error("ring RAM cell mismatch");
                const double phase=std::ldexp(double(std::bit_cast<std::int32_t>(TFHEpp::tlweSymPhase<L>(cell,secret))),-32);
                error=std::max(error,std::abs(phase-(expected[b][i]?.125:-.125)));
            }
            const double input_error=b?validate(b,count,step):0;
            if(step)times[b].push_back(elapsed);
            journal<<"{\"bits\":"<<count<<",\"backend\":\""<<names[b]<<"\",\"step\":"<<step<<",\"order\":"<<order<<",\"warmup\":"<<(step?"false":"true")
                <<",\"seconds\":"<<elapsed<<",\"address_seconds\":"<<timing.address<<",\"read_seconds\":"<<timing.read<<",\"write_seconds\":"<<timing.write<<",\"refresh_seconds\":"<<timing.refresh
                <<",\"max_phase_error\":"<<error<<",\"pre_switch_error\":"<<input_error<<",\"all_cells_checked\":true,\"correct\":true";
            if(b){const BridgeStageTimes delta{stages[b-1].packing-before_stages.packing,stages[b-1].bootstrap-before_stages.bootstrap,stages[b-1].decoding-before_stages.decoding,stages[b-1].switching-before_stages.switching};journal<<",\"stages\":";delta.json(journal);}
            journal<<"}\n";journal.flush();
            if(!journal)throw std::runtime_error("ring RAM journal write failed");
            std::cout<<"ring RAM bits="<<count<<" step="<<step<<" backend="<<names[b]<<" seconds="<<elapsed<<" all cells passed"<<std::endl;
        }
    }
    std::ofstream summary(prefix.string()+".json");summary<<std::setprecision(17)<<"{\"bits\":"<<count<<",\"warmups_per_backend\":1,\"measured_accesses_per_backend\":4,\"all_cells_checked\":true,\"security_validated\":false,\"backends\":[";
    for(std::size_t b=0;b<4;++b){if(b)summary<<',';summary<<"{\"name\":\""<<names[b]<<"\",\"access_seconds\":[";for(std::size_t i=0;i<times[b].size();++i){if(i)summary<<',';summary<<times[b][i];}summary<<"]";
        if(b){const BridgeStageTimes measured{stages[b-1].packing-start[b-1].packing,stages[b-1].bootstrap-start[b-1].bootstrap,stages[b-1].decoding-start[b-1].decoding,stages[b-1].switching-start[b-1].switching};summary<<",\"measured_stage_totals\":";measured.json(summary);summary<<",\"bridge\":";describe(b,summary);}summary<<'}';}
    summary<<"],\"combined_peak_rss_bytes\":"<<peakRSSBytes()<<"}\n";if(!summary)throw std::runtime_error("ring RAM summary write failed");return 0;
}
} // namespace ckks_ram
