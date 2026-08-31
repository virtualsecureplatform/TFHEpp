#include <params.hpp>

int main()
{
    using namespace TFHEpp::shallowboot;
    if (!isStructuredOneHotCompatible(structured_binary_gate_std128) ||
        !isStructuredOneHotCompatible(structured_binary_gate_medium) ||
        !isStructuredOneHotCompatible(structured_binary_function_std128))
        return 1;
    if (general_binary_gate_std128.lwe_hamming_weight != 43 ||
        binary_ntt_paper_std128.rlwe_modulus_log2 != 105 ||
        binary_ntt_paper_std128.key_switch_digits != 4 ||
        binary_ntt_source_screened.lwe_hamming_weight != 37)
        return 1;
}
