#include "scheme_switch_roundtrip.hpp"
#include "scheme_switch_test_params.hpp"

namespace {
// Deliberately small test parameters, not a security claim. Keep enough
// precision for both the rectangular transform and periodic reduction.
struct QuickCKKS : TFHEpp::lvl5param {
    static constexpr std::uint32_t nbit = 7, n = 1 << nbit;
    static constexpr double ckks_min_noise_stddev = 3.2;
    static constexpr std::uint32_t B̅gbit = 8, B̅gₐbit = 8;
    static constexpr std::uint32_t l̅ = 56, l̅ₐ = 56;
};
struct QuickFullPrecision : scheme_switch_test::FullCKKS {
    static constexpr std::uint32_t nbit = 7, n = 1 << nbit;
};
struct QuickTFHE : TFHEpp::lvlhalfparam {
    static constexpr std::uint32_t n = 8;
};
}

int main(int argc, char **argv)
try {
    if (argc == 2 && std::string_view(argv[1]) == "--full-precision")
        return scheme_switch_test::run<QuickFullPrecision, QuickTFHE,
                                       1000, 58, 60, 31, 9>();
    if (argc != 1) throw std::invalid_argument("usage: ckks_tfhe_roundtrip_quick [--full-precision]");
    return scheme_switch_test::run<QuickCKKS, QuickTFHE, 448, 40, 40, 15, 4>();
}
catch (const std::exception &e) {
    std::cerr << e.what() << std::endl;
    return 1;
}
