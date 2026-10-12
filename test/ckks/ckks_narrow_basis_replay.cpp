#include "ckks_narrow_basis_benchmark.hpp"

int main(int argc, char **argv)
try
{
    if (argc == 3 && (std::string(argv[1]) == "--small-lean-basis-storage" ||
                      std::string(argv[1]) == "--lean-basis-storage")) {
        if (std::filesystem::exists(argv[2]))
            throw std::invalid_argument("storage report already exists");
        return std::string(argv[1]) == "--small-lean-basis-storage"
            ? narrow_benchmark::basisBenchmark<9, true>("", "", argv[2], true)
            : narrow_benchmark::basisBenchmark<15, true>("", "", argv[2], true);
    }
    const bool create = argc >= 2 && (std::string(argv[1]) == "--small-lean-basis-fbs-checkpoint" ||
                                     std::string(argv[1]) == "--lean-basis-fbs-checkpoint");
    if ((argc != 5 && argc != 6) || (create ? argc != 5 :
            (std::string(argv[1]) != "--small-lean-basis-fbs-replay" &&
             std::string(argv[1]) != "--lean-basis-fbs-replay" &&
             std::string(argv[1]) != "--small-lean-basis-fbs-replay-secondary11" &&
             std::string(argv[1]) != "--lean-basis-fbs-replay-secondary11" &&
             std::string(argv[1]) != "--small-lean-basis-fbs-replay-secondary11-selector56" &&
             std::string(argv[1]) != "--lean-basis-fbs-replay-secondary11-selector56"))) {
        std::cerr << "Usage: ckks_narrow_basis_replay --[small-]lean-basis-fbs-replay "
                     "sigmoid|gelu|tanh|all CHECKPOINT_DIRECTORY REPORT.json [PREPARED_SELECTOR_DIRECTORY]\n"
                     "       ckks_narrow_basis_replay --[small-]lean-basis-fbs-replay-secondary11 "
                     "sigmoid|gelu|tanh|all CHECKPOINT_DIRECTORY REPORT.json [PREPARED_SELECTOR_DIRECTORY]\n"
                     "       ckks_narrow_basis_replay --[small-]lean-basis-fbs-replay-secondary11-selector56 "
                     "sigmoid|gelu|tanh|all CHECKPOINT_DIRECTORY REPORT.json [PREPARED_SELECTOR_DIRECTORY]\n"
                     "       ckks_narrow_basis_replay --[small-]lean-basis-fbs-checkpoint "
                     "sigmoid|gelu|tanh|all NEW_CHECKPOINT_DIRECTORY REPORT.json\n"
                     "       ckks_narrow_basis_replay --[small-]lean-basis-storage REPORT.json\n";
        return 2;
    }
    const std::string function = argv[2];
    if (function != "all" && function != "sigmoid" && function != "gelu" && function != "tanh")
        throw std::invalid_argument("unsupported retained basis function");
    if (std::string(argv[3]).empty() || std::string(argv[4]).empty())
        throw std::invalid_argument("checkpoint and report paths must be nonempty");
    if (std::filesystem::exists(argv[4]))
        throw std::invalid_argument("replay report already exists");
    const auto parent = std::filesystem::path(argv[4]).parent_path();
    const auto suffix = create ? "-initial.bin" : "-replay.bin";
    for (const auto &name : {"sigmoid", "gelu", "tanh"})
        if ((function == "all" || function == name) &&
            std::filesystem::exists(parent / (std::string(name) + suffix)))
            throw std::invalid_argument("basis ciphertext already exists");
    const bool small = std::string(argv[1]).starts_with("--small-");
    if (create) {
        if (std::filesystem::exists(argv[3]))
            throw std::invalid_argument("checkpoint directory already exists; use replay to reuse its keys");
        const auto selectors = (std::filesystem::path(argv[3]) / "selectors").string();
        return small
            ? narrow_benchmark::basisBenchmark<9, true>(function, selectors, argv[4], false, argv[3])
            : narrow_benchmark::basisBenchmark<15, true>(function, selectors, argv[4], false, argv[3]);
    }
    // Check only the public profile before allocating the full evaluation key.
    std::ifstream checkpoint(std::filesystem::path(argv[3]) / "checkpoint.bin", std::ios::binary);
    if (!checkpoint) throw std::invalid_argument("cannot open retained basis checkpoint");
    cereal::PortableBinaryInputArchive archive(checkpoint);
    std::array<std::uint64_t, 12> metadata{};
    archive(metadata);
    const std::array<std::uint64_t, 12> expected{
        0x4e41525242433031ULL, small ? 512ULL : 32768ULL, 576, 518, 47, 31, 85, 44, 18, 8, 44, 37};
    if (metadata != expected) throw std::invalid_argument("retained basis profile mismatch");
    checkpoint.close();
    const auto selectors = (std::filesystem::path(argv[3]) / "selectors").string();
    if (std::string(argv[1]).ends_with("-selector56")) {
        if (small)
            return narrow_benchmark::basisBenchmark<9, true, true, true>(argv[2], selectors, argv[4], false, argv[3], true, argc == 6 ? argv[5] : "");
        return narrow_benchmark::basisBenchmark<15, true, true, true>(argv[2], selectors, argv[4], false, argv[3], true, argc == 6 ? argv[5] : "");
    }
    if (std::string(argv[1]).ends_with("-secondary11")) {
        if (small)
            return narrow_benchmark::basisBenchmark<9, true, true>(argv[2], selectors, argv[4], false, argv[3], true, argc == 6 ? argv[5] : "");
        return narrow_benchmark::basisBenchmark<15, true, true>(argv[2], selectors, argv[4], false, argv[3], true, argc == 6 ? argv[5] : "");
    }
    if (small)
        return narrow_benchmark::basisBenchmark<9, true>(argv[2], selectors, argv[4], false, argv[3], true, argc == 6 ? argv[5] : "");
    return narrow_benchmark::basisBenchmark<15, true>(argv[2], selectors, argv[4], false, argv[3], true, argc == 6 ? argv[5] : "");
}
catch (const std::exception &error)
{
    std::cerr << "retained basis replay: " << error.what() << '\n';
    return 2;
}
