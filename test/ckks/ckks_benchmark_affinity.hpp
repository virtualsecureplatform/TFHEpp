#pragma once
#include <algorithm>
#include <array>
#include <atomic>
#include <cerrno>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <vector>
#ifdef __linux__
#include <sched.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

namespace narrow_benchmark {
using CPUProfile=std::vector<int>;
// Record entry CPUs without serializing the measured independent tree work.
class CPURecorder {
#ifdef __linux__
    std::array<std::atomic<bool>,CPU_SETSIZE> seen_{};
#endif
public:
    void record(){
#ifdef __linux__
        const auto cpu=sched_getcpu();
        if(cpu>=0 && cpu<CPU_SETSIZE)seen_[cpu].store(true,std::memory_order_relaxed);
#endif
    }
    CPUProfile cpus()const{
        CPUProfile result;
#ifdef __linux__
        for(int cpu=0;cpu<CPU_SETSIZE;++cpu)if(seen_[cpu].load(std::memory_order_relaxed))result.push_back(cpu);
#endif
        return result;
    }
};
inline int benchmarkThreads(){
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}
// Use one logical CPU per physical core first, then its SMT siblings. This
// respects the caller's allowed CPU set and never changes system-wide policy.
inline std::vector<CPUProfile> coreSweepProfiles(){
#if defined(__linux__) && defined(_OPENMP)
    cpu_set_t allowed;CPU_ZERO(&allowed);
    if(omp_get_proc_bind()==omp_proc_bind_close){
        // Native binding may pin the initial thread before main. Recover the
        // complete permitted set from runtime places, rather than that one pin.
        for(int place=0;place<omp_get_num_places();++place){
            std::vector<int> cpus(omp_get_place_num_procs(place));
            omp_get_place_proc_ids(place,cpus.data());
            for(int cpu:cpus)if(cpu>=0&&cpu<CPU_SETSIZE)CPU_SET(cpu,&allowed);
        }
    }else if(sched_getaffinity(0,sizeof(allowed),&allowed))throw std::runtime_error("cannot read benchmark CPU affinity");
    std::map<std::pair<int,int>,CPUProfile> cores;
    for(int cpu=0;cpu<CPU_SETSIZE;++cpu)if(CPU_ISSET(cpu,&allowed)){
        const auto path="/sys/devices/system/cpu/cpu"+std::to_string(cpu)+"/topology/";
        int socket=0,core=0;
        if(!(std::ifstream(path+"physical_package_id")>>socket) || !(std::ifstream(path+"core_id")>>core))
            throw std::runtime_error("cannot read benchmark physical CPU topology");
        cores[{socket,core}].push_back(cpu);
    }
    if(cores.empty())throw std::runtime_error("no benchmark CPUs available");
    CPUProfile physical,all;
    for(const auto &[id,cpus]:cores)physical.push_back(cpus.front());
    all=physical;
    for(const auto &[id,cpus]:cores)all.insert(all.end(),cpus.begin()+1,cpus.end());
    CPUProfile baseline(physical.begin(),physical.begin()+std::min(std::size_t{8},physical.size()));
    return {baseline,physical,all,all,physical,baseline};
#else
    throw std::invalid_argument("thread sweep requires Linux and OpenMP");
#endif
}
inline void applyCPUProfile(const CPUProfile&cpus){
#if defined(__linux__) && defined(_OPENMP)
    const auto binding=omp_get_proc_bind();
    if(cpus.empty() || (binding!=omp_proc_bind_false && binding!=omp_proc_bind_close))
        throw std::invalid_argument("core replay requires CPUs and OMP_PROC_BIND=false or close");
    omp_set_dynamic(0);omp_set_num_threads(cpus.size());
    std::vector<int> errors(cpus.size());int actual=0;
#pragma omp parallel
    {
        const auto tid=omp_get_thread_num();
        if(tid==0)actual=omp_get_num_threads();
        cpu_set_t mask;CPU_ZERO(&mask);CPU_SET(cpus[tid],&mask);
        if(binding==omp_proc_bind_false){
            if(sched_setaffinity(0,sizeof(mask),&mask))errors[tid]=errno;
        }else{
            cpu_set_t observed;CPU_ZERO(&observed);
            if(sched_getaffinity(0,sizeof(observed),&observed))errors[tid]=errno;
            else if(CPU_COUNT(&observed)!=1||!CPU_ISSET(cpus[tid],&observed))errors[tid]=EINVAL;
        }
    }
    if(actual!=static_cast<int>(cpus.size()) || std::any_of(errors.begin(),errors.end(),[](int e){return e!=0;}))
        throw std::runtime_error("cannot apply CPU profile; close binding requires ordered singleton CPU places");
#else
    throw std::invalid_argument("thread sweep requires Linux and OpenMP");
#endif
}
inline std::string CPUProfileJSON(const CPUProfile&cpus){
    std::ostringstream out;out<<'[';
    for(std::size_t i=0;i<cpus.size();++i){if(i)out<<',';out<<cpus[i];}
    out<<']';return out.str();
}
} // namespace narrow_benchmark
