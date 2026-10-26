#pragma once
#include "ckks_ram_output_switch.hpp"
#include "ckks_ram_crossover.hpp"

namespace ckks_ram {
struct RAMRunOptions {
    std::string output_switch;
    std::filesystem::path reports;
    bool capture_only=false;
};
template<class Narrow,std::size_t A=2,std::size_t W=2,class Transform>
int run(const std::filesystem::path &root,const TFHEpp::Key<P> &sparse,
        const TFHEpp::CKKSNarrowBootstrapKey<Narrow> &bootstrap,
        Transform slots_to_coefficients,bool paired=false,
        const std::filesystem::path &checkpoint={},bool ring_comparison=false,
        const RAMRunOptions *options=nullptr)
{
    using namespace TFHEpp;
    using Clock=std::chrono::steady_clock;
    auto started=Clock::now();
    auto owner=std::make_unique<SavedOwner>();load(root/"owner.bin",*owner);
    auto ek=std::make_unique<EvalKey>();load(root/"tfhe-eval.bin",*ek);
    auto fixture=std::make_unique<Fixture<A,W>>();load(root/("trace-"+std::to_string(std::size_t{1}<<A)+"x"+std::to_string(W)+".bin"),*fixture);
    auto packing=std::make_unique<PackKey>();load(root/"pack48.bin",*packing);
    PreparedTLWECoefficientPacking48<P,L> prepared(*packing);packing.reset();
    const auto reports=options&&!options->reports.empty()?options->reports:reportDirectory(root);
    std::filesystem::create_directories(reports);
    const auto source_id=checkpoint.empty()?std::string{}:optionalFingerprint(checkpoint/"checkpoint.bin");
    const auto target_id=optionalFingerprint(root/"owner.bin");
    std::array<std::unique_ptr<OutputBridge>,3> bridges;
    if(ring_comparison){
        bridges[0]=std::make_unique<OutputBridge>(root,"direct",source_id);
        bridges[1]=std::make_unique<OutputBridge>(root,"direct256",source_id);
        bridges[2]=std::make_unique<OutputBridge>(root,"ring",source_id);
    }else bridges[0]=std::make_unique<OutputBridge>(root,options&&!options->output_switch.empty()?options->output_switch:outputSwitchMode(),source_id);
    auto dense=checkpoint.empty()?std::unique_ptr<Key<P>>{}:loadDenseCheckpoint(checkpoint);
    std::unique_ptr<OutputInput> last_input;
    double pre_switch_error=0;
    std::array<BridgeStageTimes,3> stages{};
    auto &totals=stages[0];

    const double bridge_setup=std::chrono::duration<double>(Clock::now()-started).count();
    std::cout<<"RAM bridge keys loaded/prepared in "<<bridge_setup<<" seconds"<<std::endl;
    const auto &secret=owner->secret();
    // Audit the public packing evaluator against independent source phases.
    // This check is outside the timed, secret-free refresh callback.
    std::vector<TLWE<L>> audit;
    for(const auto &address:fixture->addresses)for(const auto &bit:address)audit.push_back(bit);
    auto packed=std::make_unique<CKKSCiphertext<P,48,44>>();prepared.apply(*packed,audit);
    auto polynomial=std::make_unique<std::array<double,P::n>>();
    ckksDecrypt<P,48,44>(*polynomial,*packed,sparse);
    double packing_error=0;
    for(std::size_t i=0;i<audit.size();++i){
        const auto phase=std::bit_cast<std::int32_t>(tlweSymPhase<L>(audit[i],secret));
        packing_error=std::max(packing_error,std::abs((*polynomial)[i]-std::ldexp(double(phase),-28)));
    }
    if(!std::isfinite(packing_error)||packing_error>1e-4)throw std::runtime_error("TLWE packing phase audit failed");
    std::cout<<"packing phase error="<<packing_error<<std::endl;
    auto &pack_seconds=totals.packing,&bootstrap_seconds=totals.bootstrap,&decode_seconds=totals.decoding,&switch_seconds=totals.switching;
    std::size_t refreshes=0;
    // This callable only uses public evaluation material and ciphertexts.
    const auto ckks_refresh_for=[&](const std::vector<TLWE<L>> &in,OutputBridge &bridge,BridgeStageTimes &stage){
        auto before=Clock::now();
        auto coefficients=std::make_unique<CKKSCiphertext<P,48,44>>();prepared.apply(*coefficients,in);
        stage.packing+=std::chrono::duration<double>(Clock::now()-before).count();before=Clock::now();
        // A TLWE Boolean phase +/-1/8 becomes x=+/-2 at scale 2^44.
        // sin(4*pi*x/16) restores +/-1 with zero first derivative there.
        const std::array<CKKSNarrowBlindRotationRequest<Narrow,44>,1> requests{{{coefficients.get(),2,{0,-0.5}}}};
        std::vector<typename Narrow::Output> roots;
        CKKSNarrowBlindRotateBatch<Narrow,44>(roots,requests,bootstrap,1);
        auto conjugate=std::make_unique<typename Narrow::Output>();
        bootstrap.automorphism_evaluator(conjugate->ct,roots[0].ct,Narrow::output_log_q,2*P::n-1);
        CKKSAddInPlace<P>(roots[0],*conjugate);
        auto slots=std::make_unique<CKKSCiphertext<P,93,44>>();
        static_assert(Narrow::output_log_q>=93 && Narrow::log_delta==44);
        for(std::size_t c=0;c<2;++c)for(std::size_t i=0;i<P::n;++i)
            slots->ct[c][i]=ckks_detail::reduceToLevel<P,93>(roots[0].ct[c][i]);
        stage.bootstrap+=std::chrono::duration<double>(Clock::now()-before).count();before=Clock::now();
        auto decoded=std::make_unique<CKKSCiphertext<P,48,44>>();slots_to_coefficients(*decoded,*slots);
        auto reduced=std::make_unique<CKKSCiphertext<P,33,29>>();CKKSScaleDown<P,48,44,15>(*reduced,*decoded);
        auto extracted=std::make_unique<CKKSCiphertext<P,32,29>>();
        for(std::size_t c=0;c<2;++c)for(std::size_t i=0;i<P::n;++i)
            extracted->ct[c][i]=ckks_detail::reduceToLevel<P,32>(reduced->ct[c][i]);
        stage.decoding+=std::chrono::duration<double>(Clock::now()-before).count();before=Clock::now();
        std::vector<TLWE<L>> out;bridge.apply(out,*extracted,in.size());
        stage.switching+=std::chrono::duration<double>(Clock::now()-before).count();++refreshes;
        last_input=std::move(extracted);
        return out;
    };
    const auto ckks_refresh=[&](const auto &in){return ckks_refresh_for(in,*bridges[0],stages[0]);};
    const auto tfhe_refresh=[&](const std::vector<TLWE<L>> &in){
        std::vector<TLWE<L>> out(in.size());
        // Match the existing RAM's serial per-cell refresh loop.
        for(std::size_t i=0;i<in.size();++i){
            TLWE<typename KS::targetP> switched;TRLWE<L> cell;
            IdentityKeySwitch<KS>(switched,in[i],ek->getiksk<KS>());
            BlindRotate<BR>(cell,switched,ek->getbkfft<BR>(),μpolygen<L,L::μ>());
            SampleExtractIndex<L>(out[i],cell,0);
        }
        return out;
    };
    const auto validate_input=[&](std::size_t count){
        if(!dense||!last_input)throw std::runtime_error("missing retained source for bridge audit");
        if(count==P::n/2&&!std::filesystem::exists(reports/"output-input32.bin"))saveOutputInput(reports/"output-input32.bin",*last_input,source_id,target_id);
        const auto error=auditInput(*last_input,*dense,count);pre_switch_error=std::max(pre_switch_error,error);
        return error;
    };
    const auto capture=std::getenv("CKKS_RAM_CAPTURE_ONLY");
    if((options&&options->capture_only)||(capture&&std::string(capture)=="1")){
        if constexpr((std::size_t{1}<<A)*W==P::n/2){
            std::vector<TLWE<L>> cells(P::n/2);
            for(std::size_t w=0;w<W;++w)for(std::size_t i=0;i<(std::size_t{1}<<A);++i)
                SampleExtractIndex<L>(cells[w*(std::size_t{1}<<A)+i],fixture->initial[w][i],0);
            const auto out=ckks_refresh(cells);
            for(std::size_t i=0;i<out.size();++i)if(tlweSymDecrypt<L>(out[i],secret)!=tlweSymDecrypt<L>(cells[i],secret))throw std::runtime_error("capture bridge bit mismatch");
            validate_input(cells.size());
            std::ofstream capture_report(reports/"capture.json");capture_report<<std::setprecision(17)<<"{\"security_validated\":false,\"input_reused\":true,\"fresh_private_keys_generated\":false,\"pre_switch_error\":"<<pre_switch_error<<",\"stages\":";totals.json(capture_report);capture_report<<",\"output_bridge\":";bridges[0]->json(capture_report);capture_report<<",\"combined_peak_rss_bytes\":"<<peakRSSBytes()<<",\"correct\":true}\n";
            if(!capture_report)throw std::runtime_error("capture receipt write failed");
            std::cout<<"Retained initial memory restored and q32 input captured"<<std::endl;return 0;
        }else throw std::invalid_argument("capture requires a full batch");
    }
    if(ring_comparison){
        if constexpr(W==1){return runRingComparison<A>(reports,*fixture,secret,*ek,
            [&](std::size_t b,const auto &in){return b?ckks_refresh_for(in,*bridges[b-1],stages[b-1]):tfhe_refresh(in);},
            [&](std::size_t,std::size_t count,std::size_t){return validate_input(count);},
            [&](std::size_t b,std::ostream&out){bridges[b-1]->json(out);},stages);}
        else throw std::invalid_argument("ring comparison requires one-bit words");
    }
    if(paired)return runPaired<A,W>(reports,*fixture,secret,*ek,tfhe_refresh,ckks_refresh,
        pack_seconds,bootstrap_seconds,decode_seconds,switch_seconds,
        [&](std::size_t count){return dense?validate_input(count):0.;},
        [&](std::ostream&out){bridges[0]->json(out);});
    const auto plain=trace<A,W>();
    std::ofstream report(reports/(A==2?"comparison.json":"comparison-256x1.json"));report<<std::setprecision(17);
    report<<"{\"address_bits\":"<<A<<",\"word_bits\":"<<W<<",\"operations\":"<<plain.size()
          <<",\"read_before_write\":true,\"fixture_reused\":true,\"keys_generated_during_evaluation\":0"
          <<",\"bridge_setup_seconds\":"<<bridge_setup<<",\"packing_phase_max_error\":"<<packing_error
          <<",\"security_validated\":false,\"backends\":[";
    for(int backend=0;backend<2;++backend){
        auto ram=std::make_unique<typename Fixture<A,W>::RAM>();ram->memory=fixture->initial;
        std::array<std::uint32_t,std::size_t{1}<<A> expected{};
        for(std::size_t i=0;i<expected.size();++i)expected[i]=i&((std::uint32_t{1}<<W)-1);
        double total=0,address=0,read=0,write=0,refresh=0,max_phase_error=0;
        typename Fixture<A,W>::RAM::Word previous_output{};
        for(std::size_t step=0;step<plain.size();++step){
            typename Fixture<A,W>::RAM::Timing timing;
            const auto before=Clock::now();
            const auto refreshed=[&](const auto &in){return backend?ckks_refresh(in):tfhe_refresh(in);};
            // Step 8 reads the old value 1 before writing 2; feed that
            // encrypted output straight back as step 9's write data.
            const auto &write_data=A==2&&step==9?previous_output:fixture->data[step];
            const auto output=ram->access(fixture->addresses[step],fixture->enables[step],write_data,*ek,refreshed,&timing);
            previous_output=output;
            const double elapsed=std::chrono::duration<double>(Clock::now()-before).count();total+=elapsed;
            if(backend && dense)validate_input((std::size_t{1}<<A)*W);
            address+=timing.address;read+=timing.read;write+=timing.write;refresh+=timing.refresh;
            for(std::size_t w=0;w<W;++w)if(tlweSymDecrypt<L>(output[w],secret)!=bool((expected[plain[step].address]>>w)&1))throw std::runtime_error("RAM read mismatch at step "+std::to_string(step));
            if(plain[step].enable)expected[plain[step].address]=plain[step].data;
            // Validate every cell after every access, including untouched cells.
            for(std::size_t w=0;w<W;++w)for(std::size_t i=0;i<expected.size();++i){
                TLWE<L> cell;SampleExtractIndex<L>(cell,ram->memory[w][i],0);
                const bool bit=(expected[i]>>w)&1;
                if(tlweSymDecrypt<L>(cell,secret)!=bit)throw std::runtime_error("RAM state mismatch at step "+std::to_string(step)+" cell "+std::to_string(i));
                const double phase=std::ldexp(double(std::bit_cast<std::int32_t>(tlweSymPhase<L>(cell,secret))),-32);
                max_phase_error=std::max(max_phase_error,std::abs(phase-(bit?0.125:-0.125)));
            }
            std::cout<<(backend?"CKKS":"TFHE")<<" RAM step="<<step<<" seconds="<<elapsed<<" all cells passed"<<std::endl;
        }
        if(backend)report<<',';
        report<<"{\"name\":\""<<(backend?"ckks-batch":"tfhe-cell")<<"\",\"correct\":true,\"seconds\":"<<total
              <<",\"address_seconds\":"<<address<<",\"read_seconds\":"<<read<<",\"write_seconds\":"<<write
              <<",\"refresh_seconds\":"<<refresh<<",\"state_max_phase_error\":"<<max_phase_error<<'}';report.flush();
    }
    report<<"],\"ckks_refreshes\":"<<refreshes<<",\"pack_seconds\":"<<pack_seconds<<",\"bootstrap_seconds\":"<<bootstrap_seconds
          <<",\"slots_to_coefficients_seconds\":"<<decode_seconds<<",\"to_tlwe_seconds\":"<<switch_seconds
          <<",\"combined_peak_rss_bytes\":"<<peakRSSBytes()<<",\"pre_switch_max_error\":"<<pre_switch_error<<",\"output_bridge\":";bridges[0]->json(report);report<<"}\n";
    return 0;
}
} // namespace ckks_ram
