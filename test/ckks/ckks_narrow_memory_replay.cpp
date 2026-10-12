#include <ckks/narrow_refresh.hpp>
#include <ckks/narrow_selector_file.hpp>
#include <cereal/archives/portable_binary.hpp>
#include <fstream>
#include <iostream>

namespace {
// One fixed small plain37 profile keeps sanitizer builds independent of the
// benchmark's many full-ring template instantiations.
struct P {
    static constexpr int32_t key_value_max = 1, key_value_min = -1;
    static constexpr std::uint32_t nbit = 9, n = 1 << nbit, k = 1;
    static constexpr std::uint32_t l = 1, l_a = 1, lₐ = 1;
    static constexpr std::uint32_t Bgbit = 4, Bg_abit = 4, Bgₐbit = 4;
    using T = TFHEpp::MultiLimbUInt<7>;
    static constexpr T Bg = T{1} << Bgbit, Bg_a = T{1} << Bg_abit, Bgₐ = Bg_a;
    static constexpr auto errordist = TFHEpp::ErrorDistribution::ModularGaussian;
    static const inline double α = std::ldexp(1.0, -448);
    static constexpr double ckks_min_noise_stddev = 3.2;
    static constexpr T μ = T{1} << 253;
    static constexpr std::uint32_t plain_modulusbit = 20;
    static constexpr T plain_modulus = T{786433};
    static constexpr double Δ = 0.0;
    static constexpr std::uint32_t l̅ = 112, l̅ₐ = l̅, B̅gbit = 4, B̅gₐbit = 4;
};
using Narrow = TFHEpp::CKKSNarrowBootstrapSchedule<P, 48, 385, 48, 37, 31, true, 28, true>;
using Refresh = TFHEpp::CKKSNarrowRefreshSchedule<Narrow, 44, 37, 0, true>;

void replay(const std::filesystem::path &directory)
{
    using namespace TFHEpp;
    auto key = std::make_unique<CKKSNarrowRefreshKey<Refresh>>();
    auto input = std::make_unique<Refresh::Input>();
    Key<P> dense{}, sparse{};
    CKKSSlotVector<P> slots{};
    std::ifstream file(directory / "checkpoint.bin", std::ios::binary);
    if (!file) throw std::runtime_error("cannot open small plain37 checkpoint");
    {
        cereal::PortableBinaryInputArchive archive(file);
        std::array<std::uint64_t, 9> metadata{};
        std::string function;
        archive(metadata, function);
        const std::array<std::uint64_t, 9> expected{0x4e41525243503031ULL, P::n,
            Refresh::input_log_q, Refresh::input_log_delta, Narrow::selector_log_q,
            Narrow::log_delta, Narrow::plain_log_delta, Narrow::output_log_q, P::T::digits};
        if (metadata != expected || function != "sigmoid")
            throw std::runtime_error("memory replay requires the small plain37 sigmoid profile");
        archive(dense, sparse, slots, *input, *key);
    }
    if (file.peek() != std::char_traits<char>::eof())
        throw std::runtime_error("trailing small checkpoint data");
    key->narrow.selector_source = std::make_shared<CKKSNarrowSelectorFile<Narrow>>(
        directory / "selectors", typename CKKSNarrowSelectorFile<Narrow>::Existing{});
    CKKSNarrowUseStagedFFT<Narrow>(key->narrow, 7);
    std::vector<std::complex<double>> coefficients(13);
    coefficients[0] = 0.125;
    for (std::size_t k = 1; k < coefficients.size(); ++k)
        coefficients[k] = {0.005 * (static_cast<int>(k % 3) - 1) / k,
                           (k % 2 ? -0.003 : 0.003) / k};
    auto first = std::make_unique<Refresh::Output>();
    auto second = std::make_unique<Refresh::Output>();
    CKKSNarrowFunctionalRefresh<Refresh>(*first, *input, *key, coefficients, true, 12);
    CKKSNarrowFunctionalRefresh<Refresh>(*second, *input, *key, coefficients, true, 1);
    if (first->ct != second->ct)
        throw std::runtime_error("changing tiles changed the retained-state ciphertext");
    CKKSSlotVector<P> decoded{};
    ckksSlotDecrypt<P>(decoded, *first, dense);
    const auto target = [&](double x) {
        double result = coefficients[0].real();
        for (std::size_t k = 1; k < coefficients.size(); ++k)
            result += 2 * (coefficients[k] * std::polar(1.0,
                2 * std::acos(-1.0) * k * x / Refresh::gap())).real();
        return result;
    };
    double error = 0;
    for (std::size_t i = 0; i < decoded.size(); ++i) {
        const double candidate = std::abs(decoded[i] -
            std::complex<double>{target(slots[i].real()), target(slots[i].imag())});
        if (!std::isfinite(candidate)) throw std::runtime_error("nonfinite memory replay result");
        error = std::max(error, candidate);
    }
    if (error >= std::ldexp(1.0, -23)) throw std::runtime_error("memory replay missed 23-bit gate");
    if (key->narrow.relin_fft || key->narrow.factor_conjugation_fft ||
        key->narrow.output_conjugation_fft || key->narrow.tail_rotations_fft)
        throw std::runtime_error("staged evaluation retained a temporary cache");
    std::cout << "retained small plain37 memory replay: tiles12/1 exact, slots="
              << decoded.size() << " max_error=" << error << '\n';
}
} // namespace

int main(int argc, char **argv)
{
    if (argc == 1) {
        std::cout << "opt-in memory replay: ckks_narrow_memory_replay CHECKPOINT_DIRECTORY\n";
        return 0;
    }
    try {
        if (argc != 2) throw std::invalid_argument("expected one checkpoint directory");
        replay(argv[1]);
        return 0;
    }
    catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
