// Encrypted RAM whose cells are refreshed by one TUTTI Boolean bootstrap per
// 4096 cells, against the per-cell TFHE refresh on the same host.
//   ckks_tutti_ram [address bits: 8, 10, 12 or 13] [measured accesses] [baseline]
// baseline: serial (default), parallel, both or none. All keys are generated
// here; decryption happens only in checks outside the timed accesses.
// 8192 cells need a large stack (ulimit -s 131072).
#include <ckks/tutti_boolean.hpp>
#include <tfhe/refreshable_ram.hpp>

#include <cstdlib>
#include <iomanip>
#include <iostream>

#ifdef TFHEPP_HAS_EXTENDED_MULTILIMB_PARAMS
namespace {
using namespace TFHEpp;
using Config = CKKSTuttiBooleanLL13Config;
using CP = Config::Param;
using KS = independentlvl10param;
using CB = cblvl02param;
using PKS = cblvl21param;
using BR = independentlvl01param;
using L = PKS::targetP;
using RefreshKey = CKKSTuttiBooleanRefreshKey<Config, L>;
using Clock = std::chrono::steady_clock;

double since(Clock::time_point start)
{
    return std::chrono::duration<double>(Clock::now() - start).count();
}

struct Request { std::uint32_t address; bool enable, data; };

struct Context {
    const Key<L> &secret;
    const Key<CP> &dense;
    const EvalKey &ek;
    const RefreshKey &key;
    int measured;
    std::string baseline;
};

// Largest distance of a refreshed coefficient from +-1/8 before the switch.
double audit(const CKKSCiphertext<CP, 32, 29> &ct, const Key<CP> &dense, std::size_t count)
{
    auto values = std::make_unique<std::array<double, CP::n>>();
    ckksDecrypt<CP, 32, 29>(*values, ct, dense);
    double worst = 0;
    for (std::size_t i = 0; i < count; ++i)
        worst = std::max(worst, std::abs(std::abs((*values)[i]) - 1.0) / 8);
    return worst;
}

template <std::size_t A>
void run(const Context &context)
{
    using RAM = RefreshableRAM<KS, CB, PKS, BR, A, 1>;
    constexpr std::size_t cells = RAM::cells;
    const auto &secret = context.secret;
    std::vector<Request> requests{{0, true, true}, {cells - 1, true, false},
        {cells / 2 - 1, false, false}, {0, false, false}, {cells - 1, false, true},
        {cells / 2, true, true}, {cells / 2, false, false}, {1, true, false}};
    requests.resize(std::min<std::size_t>(requests.size(), context.measured + 1));

    auto initial = std::make_unique<typename RAM::Memory>();
    for (std::size_t i = 0; i < cells; ++i) {
        Polynomial<L> message{};
        message[0] = i & 1 ? L::μ : -L::μ;
        trlweSymEncrypt<L>((*initial)[0][i], message, secret);
    }
    std::vector<typename RAM::Address> addresses(requests.size());
    std::vector<TLWE<L>> enables(requests.size());
    std::vector<typename RAM::Word> data(requests.size());
    for (std::size_t r = 0; r < requests.size(); ++r) {
        for (std::size_t i = 0; i < A; ++i)
            tlweSymEncrypt<L>(addresses[r][i], (requests[r].address >> i) & 1 ? L::μ : -L::μ, secret);
        tlweSymEncrypt<L>(data[r][0], requests[r].data ? L::μ : -L::μ, secret);
        tlweSymEncrypt<L>(enables[r], requests[r].enable ? L::μ : -L::μ, secret);
    }

    CKKSTuttiBooleanTimings stages;
    std::size_t batches = 0;
    double pre_switch_error = 0;
    std::vector<std::unique_ptr<CKKSCiphertext<CP, 32, 29>>> captured;
    // Public evaluation material and ciphertexts only.
    const auto tutti_refresh = [&](const std::vector<TLWE<L>> &in) {
        std::vector<TLWE<L>> out;
        out.reserve(in.size());
        captured.clear();
        for (std::size_t offset = 0; offset < in.size(); offset += RefreshKey::batch) {
            const auto count = std::min<std::size_t>(RefreshKey::batch, in.size() - offset);
            const std::vector<TLWE<L>> part(in.begin() + offset, in.begin() + offset + count);
            std::vector<TLWE<L>> refreshed;
            captured.push_back(std::make_unique<CKKSCiphertext<CP, 32, 29>>());
            CKKSTuttiBooleanRefresh(refreshed, part, context.key, &stages, captured.back().get());
            out.insert(out.end(), refreshed.begin(), refreshed.end());
            ++batches;
        }
        return out;
    };
    const auto tfhe_cell = [&](TLWE<L> &out, const TLWE<L> &in) {
        TLWE<typename KS::targetP> switched;
        TRLWE<L> cell;
        IdentityKeySwitch<KS>(switched, in, context.ek.template getiksk<KS>());
        BlindRotate<BR>(cell, switched, context.ek.template getbkfft<BR>(), μpolygen<L, L::μ>());
        SampleExtractIndex<L>(out, cell, 0);
    };
    const auto tfhe_serial = [&](const std::vector<TLWE<L>> &in) {
        std::vector<TLWE<L>> out(in.size());
        for (std::size_t i = 0; i < in.size(); ++i) tfhe_cell(out[i], in[i]);
        return out;
    };
    const auto tfhe_parallel = [&](const std::vector<TLWE<L>> &in) {
        std::vector<TLWE<L>> out(in.size());
#pragma omp parallel for schedule(static)
        for (std::size_t i = 0; i < in.size(); ++i) tfhe_cell(out[i], in[i]);
        return out;
    };

    std::vector<std::string> backends{"tutti"};
    if (context.baseline == "serial" || context.baseline == "both") backends.push_back("tfhe-serial");
    if (context.baseline == "parallel" || context.baseline == "both") backends.push_back("tfhe-parallel");
    std::cout << std::setprecision(4) << "RAM " << cells << " x 1 bit, " << requests.size() - 1
              << " measured accesses after one warm-up\n";
    for (const auto &backend : backends) {
        auto ram = std::make_unique<RAM>();
        ram->memory = *initial;
        std::vector<bool> expected(cells);
        for (std::size_t i = 0; i < cells; ++i) expected[i] = i & 1;
        double total = 0, address = 0, read = 0, write = 0, refresh = 0, phase_error = 0;
        for (std::size_t step = 0; step < requests.size(); ++step) {
            if (step == 1) { stages = {}; batches = 0; }
            typename RAM::Timing timing;
            const auto start = Clock::now();
            const auto output = ram->access(addresses[step], enables[step], data[step], context.ek,
                [&](const auto &in) {
                    return backend == "tutti" ? tutti_refresh(in) :
                        backend == "tfhe-serial" ? tfhe_serial(in) : tfhe_parallel(in);
                }, &timing);
            const double seconds = since(start);
            if (tlweSymDecrypt<L>(output[0], secret) != expected[requests[step].address])
                throw std::runtime_error("RAM read mismatch");
            if (requests[step].enable) expected[requests[step].address] = requests[step].data;
            for (std::size_t i = 0; i < cells; ++i) {
                TLWE<L> cell;
                SampleExtractIndex<L>(cell, ram->memory[0][i], 0);
                if (tlweSymDecrypt<L>(cell, secret) != expected[i])
                    throw std::runtime_error("RAM state mismatch at cell " + std::to_string(i));
                const double phase = std::ldexp(static_cast<double>(
                    std::bit_cast<std::int32_t>(tlweSymPhase<L>(cell, secret))), -32);
                phase_error = std::max(phase_error, std::abs(phase - (expected[i] ? 0.125 : -0.125)));
            }
            if (backend == "tutti")
                for (std::size_t b = 0; b < captured.size(); ++b)
                    pre_switch_error = std::max(pre_switch_error, audit(*captured[b], context.dense,
                        std::min<std::size_t>(RefreshKey::batch, cells - b * RefreshKey::batch)));
            std::cout << "  " << backend << " access " << step << (step ? "" : " (warm-up)") << ": "
                      << seconds << " s, refresh " << timing.refresh << " s" << std::endl;
            if (!step) continue;
            total += seconds; address += timing.address; read += timing.read;
            write += timing.write; refresh += timing.refresh;
        }
        const double n = static_cast<double>(requests.size() - 1);
        std::cout << backend << " mean access " << total / n << " s (address " << address / n
                  << ", read " << read / n << ", write " << write / n << ", refresh " << refresh / n
                  << "), cell phase error " << phase_error << '\n';
        if (backend == "tutti") {
            const double b = static_cast<double>(batches);
            std::cout << "  per batch (" << b / n << " per access): packing "
                      << stages.packing_seconds / b << ", bootstrap " << stages.bootstrap_seconds / b
                      << " (sparse switch " << stages.bootstrap.linear_seconds / b << ", accmul "
                      << stages.bootstrap.accmul_seconds / b << ", tree "
                      << stages.bootstrap.tree_seconds / b << ", finish "
                      << stages.bootstrap.finish_seconds / b << "), decoding "
                      << stages.decoding_seconds / b << ", switching "
                      << stages.switching_seconds / b << "\n  pre-switch error " << pre_switch_error
                      << " (limit 0.0625)\n";
        }
    }
}
} // namespace
#endif

int main(int argc, char **argv)
{
#ifdef TFHEPP_HAS_EXTENDED_MULTILIMB_PARAMS
    const int address_bits = argc > 1 ? std::atoi(argv[1]) : 8;
    const int measured = argc > 2 ? std::atoi(argv[2]) : 2;
    const std::string baseline = argc > 3 ? argv[3] : "serial";
    if (measured < 1 || measured > 7) { std::cerr << "measured accesses: 1 to 7\n"; return 2; }

    auto start = Clock::now();
    auto sk = std::make_unique<SecretKey>();
    auto ek = std::make_unique<EvalKey>();
    ek->emplaceiksk<KS>(*sk);
    ek->emplacebkfft<CB>(*sk);
    ek->emplacebkfft<BR>(*sk);
    ek->emplaceprivksk4cb<PKS>(*sk);
    const auto &secret = sk->key.getIndependent<L>();
    std::cout << "TFHE keys " << since(start) << " s" << std::endl;

    start = Clock::now();
    auto dense = std::make_unique<Key<CP>>();
    auto sparse = std::make_unique<Key<CP>>();
    CKKSTuttiDenseKeyGen<CP>(*dense, Config::dense_weight);
    CKKSTuttiSparseKeyGen<RefreshKey::Schedule>(*sparse);
    auto key = std::make_unique<RefreshKey>();
    CKKSTuttiBooleanKeyGen(*key, secret, *dense, *sparse);
    std::cout << "TUTTI Boolean refresh keys " << since(start) << " s" << std::endl;

    // A full batch of random bits, checked before any RAM access.
    {
        std::vector<TLWE<L>> in(RefreshKey::batch), out;
        std::vector<bool> bits(in.size());
        std::uniform_int_distribution<int> coin(0, 1);
        for (std::size_t i = 0; i < in.size(); ++i) {
            bits[i] = coin(generator);
            tlweSymEncrypt<L>(in[i], bits[i] ? L::μ : -L::μ, secret);
        }
        CKKSTuttiBooleanTimings stages;
        auto captured = std::make_unique<CKKSCiphertext<CP, 32, 29>>();
        start = Clock::now();
        CKKSTuttiBooleanRefresh(out, in, *key, &stages, captured.get());
        const double seconds = since(start);
        double worst = 0;
        for (std::size_t i = 0; i < in.size(); ++i) {
            if (tlweSymDecrypt<L>(out[i], secret) != bits[i]) {
                std::cerr << "TUTTI Boolean refresh flipped bit " << i << '\n';
                return 1;
            }
            const double phase = std::ldexp(static_cast<double>(
                std::bit_cast<std::int32_t>(tlweSymPhase<L>(out[i], secret))), -32);
            worst = std::max(worst, std::abs(phase - (bits[i] ? 0.125 : -0.125)));
        }
        std::cout << "batch of " << in.size() << " random bits refreshed in " << seconds
                  << " s, pre-switch error " << audit(*captured, *dense, in.size())
                  << ", output phase error " << worst << std::endl;
    }

    const Context context{secret, *dense, *ek, *key, measured, baseline};
    switch (address_bits) {
    case 8: run<8>(context); break;
    case 10: run<10>(context); break;
    case 12: run<12>(context); break;
    case 13: run<13>(context); break;
    default: std::cerr << "address bits: 8, 10, 12 or 13\n"; return 2;
    }
    std::cout << "TUTTI RAM passed\n";
#else
    (void)argc; (void)argv;
#endif
}
