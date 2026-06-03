#ifdef USE_PERF
#include <gperftools/profiler.h>
#endif

#include <cereal/archives/portable_binary.hpp>
#include <cereal/types/vector.hpp>
#include <chrono>
#include <fstream>
#include <iostream>
#include <memory>
#include <random>
#include <tfhe++.hpp>

int main()
{
    // generate a random key
    using brP = TFHEpp::lvl02param;
    using xorbrP = TFHEpp::lvl01param;
    using briksP = TFHEpp::lvl10param;
    using cbiksP = TFHEpp::lvl21param;
    using privksP = TFHEpp::lvl21param;

    // To see performance
    std::chrono::system_clock::time_point start, end;
#ifdef USE_PERF
    ProfilerStart("evalkeyen.prof");
#endif
    start = std::chrono::system_clock::now();
    std::unique_ptr<TFHEpp::SecretKey> sk(new TFHEpp::SecretKey);
    TFHEpp::EvalKey ek;
    ek.emplaceiksk<briksP>(*sk);
    ek.emplaceiksk<cbiksP>(*sk);
    ek.emplacebkfft<brP>(*sk);
    ek.emplacebkfft<xorbrP>(*sk);
    ek.emplaceprivksk4cb<privksP>(*sk);
#ifdef USE_PERF
    ProfilerStop();
#endif

    end = std::chrono::system_clock::now();
    double elapsed =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
            .count();
    std::cout << elapsed << "ms" << std::endl;
}