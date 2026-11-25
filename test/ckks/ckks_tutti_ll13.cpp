// Full-ring TUTTI refresh at the LL13 schedule: precision by the paper's
// metric (inputs 2i/n - 1, -log2 of the largest slot error) and latency.
// Needs about 5 GB. Pass a repetition count as the first argument.
#include <ckks/tutti.hpp>

#include <cstdlib>
#include <iostream>

int main(int argc, char **argv)
{
#ifdef TFHEPP_HAS_EXTENDED_MULTILIMB_PARAMS
    using namespace TFHEpp;
    using Config = CKKSTuttiLL13Config;
    using Schedule = CKKSTuttiSchedule<Config>;
    using P = Schedule::Param;
    using Clock = std::chrono::steady_clock;
    const auto seconds = [](Clock::time_point from) {
        return std::chrono::duration<double>(Clock::now() - from).count();
    };
    const int repetitions = argc > 1 ? std::atoi(argv[1]) : 1;

    auto start = Clock::now();
    auto dense = std::make_unique<Key<P>>();
    auto sparse = std::make_unique<Key<P>>();
    CKKSTuttiDenseKeyGen<P>(*dense, Config::dense_weight);
    CKKSTuttiSparseKeyGen<Schedule>(*sparse);
    auto key = std::make_unique<CKKSTuttiKey<Schedule>>();
    CKKSTuttiKeyGen(*key, *dense, *sparse);
    std::cout << "N = " << P::n << ", blocks " << Schedule::blocks << ", stride "
              << Schedule::stride << ", indicators " << Schedule::indicator_count
              << " at 2^" << Schedule::selector_log_q << ", input 2^" << Schedule::input_log_q
              << ", output 2^" << Schedule::output_log_q << "\nkey generation "
              << seconds(start) << " s\n";

    auto slots = std::make_unique<CKKSSlotVector<P>>();
    auto decoded = std::make_unique<CKKSSlotVector<P>>();
    const double count = static_cast<double>(slots->size());
    double worst_overall = 0;
    for (int run = 0; run < repetitions; ++run) {
        for (std::size_t i = 0; i < slots->size(); ++i)
            (*slots)[i] = Config::message_bound * (2 * static_cast<double>(i) / count - 1);
        auto input = std::make_unique<Schedule::Input>();
        ckksSlotEncrypt<P>(*input, *slots, *dense);
        auto output = std::make_unique<Schedule::Output>();
        CKKSTuttiTimings timings;
        start = Clock::now();
        CKKSTuttiBootstrap(*output, *input, *key, &timings);
        const double total = seconds(start);
        ckksSlotDecrypt<P>(*decoded, *output, *dense);
        double worst = 0, squares = 0;
        for (std::size_t i = 0; i < slots->size(); ++i) {
            const double error = std::abs((*decoded)[i] - (*slots)[i]);
            worst = std::max(worst, error);
            squares += error * error;
        }
        worst_overall = std::max(worst_overall, worst);
        std::cout << "refresh " << total << " s (linear " << timings.linear_seconds
                  << ", accmul " << timings.accmul_seconds << ", tree " << timings.tree_seconds
                  << ", finish " << timings.finish_seconds << "), max error " << worst
                  << " = " << -std::log2(worst) << " bits, rms "
                  << -std::log2(std::sqrt(squares / count)) << " bits\n";

        // One multiplication returns to the refresh input level.
        constexpr auto q = Schedule::output_log_q;
        constexpr auto scale = Schedule::output_log_delta;
        auto relin = std::make_unique<CKKSDDRelinKeyChain<P, q + scale, scale, 1>>();
        CKKSDDRelinKeyChainGen<P, q + scale, scale, 1>(*relin, *dense);
        auto squared = std::make_unique<CKKSCiphertext<P, q - scale, scale>>();
        CKKSMultRelinThenRescale<P>(*squared, *output, *output, relin->template get<0>());
        static_assert(q - scale == Schedule::input_log_q);
        ckksSlotDecrypt<P>(*decoded, *squared, *dense);
        worst = 0;
        for (std::size_t i = 0; i < slots->size(); ++i)
            worst = std::max(worst, std::abs((*decoded)[i] - (*slots)[i] * (*slots)[i]));
        std::cout << "square after refresh: max error " << worst << " = "
                  << -std::log2(worst) << " bits\n";
        worst_overall = std::max(worst_overall, worst);
    }
    if (!(worst_overall < 0.125)) {
        std::cerr << "TUTTI LL13 refresh lost its message\n";
        return 1;
    }
    std::cout << "TUTTI LL13 refresh passed\n";
#else
    (void)argc; (void)argv;
#endif
}
