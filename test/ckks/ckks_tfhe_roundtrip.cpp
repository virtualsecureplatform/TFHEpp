#include "scheme_switch_roundtrip.hpp"
#include "scheme_switch_test_params.hpp"

namespace {
void schedule_check()
{
    using Large = TFHEpp::lvl6param;
    TFHEpp::CKKSRotationKeyIndexSet<Large> input{}, output{};
    TFHEpp::CKKSSlotsToCoeffsRotationKeyIndices<Large, 880, 50, 35>(
        input, output, 128);
    // Seven baby bits, seven giant bits, plus input conjugation.
    if (TFHEpp::CKKSRotationKeyIndexSetCount<Large>(input) != 8 ||
        TFHEpp::CKKSRotationKeyIndexSetCount<Large>(output) != 7)
        throw std::runtime_error("unexpected full-ring rotation schedule");
    TFHEpp::CKKSDirectRotationKeyIndexSet<Large> input_direct{}, output_direct{};
    TFHEpp::CKKSSlotsToCoeffsHybridRotationKeyIndices<Large, 880, 50, 35>(
        input, input_direct, output, output_direct, 128);
    std::cout << "lvl6 schedule collected without the 112 GiB encoded plans\n";
}

}  // namespace

int main(int argc, char **argv)
try {
    schedule_check();
    if (argc == 1 || (argc == 2 && std::string_view(argv[1]) == "--schedule-only")) {
        std::cout << "Use --full for the multi-GiB encrypted roundtrip." << std::endl;
        return 0;
    }
    if (argc != 2 || std::string_view(argv[1]) != "--full")
        throw std::invalid_argument("usage: ckks_tfhe_roundtrip [--schedule-only | --full]");
    std::cout << "Functional prototype only: this profile FAILS the 128-bit "
                 "security screen; see docs/SchemeSwitchingCalibration.md."
              << std::endl;
    return scheme_switch_test::run<scheme_switch_test::FullCKKS, TFHEpp::lvlhalfparam,
                                   1000, 58, 60, 31, 9>();
}
catch (const std::exception &e) {
    std::cerr << e.what() << std::endl;
    return 1;
}
