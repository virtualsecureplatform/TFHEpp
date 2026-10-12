#include <ckks/narrow_refresh.hpp>
#include <ckks/narrow_selector_file.hpp>
#include <cereal/archives/portable_binary.hpp>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <random>

namespace {
template <std::uint32_t Nbit>
struct ArchiveParam {
    static constexpr std::int32_t key_value_max = 1, key_value_min = -1;
    static constexpr std::uint32_t nbit = Nbit, n = 1 << nbit, k = 1;
    static constexpr std::uint32_t l = 1, l_a = 1, lₐ = 1;
    static constexpr std::uint32_t Bgbit = 4, Bg_abit = 4, Bgₐbit = 4;
    using T = TFHEpp::MultiLimbUInt<7>;
    static constexpr T Bg = T{1} << Bgbit, Bg_a = Bg, Bgₐ = Bg;
    static constexpr TFHEpp::ErrorDistribution errordist = TFHEpp::ErrorDistribution::ModularGaussian;
    static const inline double α = std::ldexp(1.0, -448);
    static constexpr double ckks_min_noise_stddev = 3.2;
    static constexpr T μ = T{1} << 253;
    static constexpr std::uint32_t plain_modulusbit = 20;
    static constexpr T plain_modulus = T{786433};
    static constexpr double Δ = 0.0;
    static constexpr std::uint32_t l̅ = 112, l̅ₐ = l̅, B̅gbit = 4, B̅gₐbit = 4;
};


template <std::uint32_t Nbit>
struct WideParam : ArchiveParam<Nbit> {
    static constexpr std::uint32_t l̅ = 56, l̅ₐ = l̅, B̅gbit = 8, B̅gₐbit = 8;
};
template <std::uint32_t Nbit, bool Auxiliary8 = false>
int run(const std::filesystem::path &state, const std::filesystem::path &destination)
{
    using namespace TFHEpp;
    using P = ArchiveParam<Nbit>;
    using W = WideParam<Nbit>;
    using T = typename P::T;
    using Raw = CKKSAutoKey<P, 356>;
    using WideRaw = CKKSAutoKey<W, 356>;
    constexpr std::uint32_t q = Auxiliary8 ? 348 : 356;
    constexpr std::uint32_t delta = Auxiliary8 ? 48 : 49;
    using Cipher = CKKSCiphertext<P, q, delta>;
    using WideCipher = CKKSCiphertext<W, q, delta>;
    std::ifstream file(state, std::ios::binary);
    if (!file) throw std::runtime_error("cannot open retained auxiliary8 state");
    cereal::PortableBinaryInputArchive archive(file);
    std::array<std::uint64_t, 10> metadata{};
    std::string function; bool finished{};
    archive(metadata, function, finished);
    if (metadata != std::array<std::uint64_t, 10>{0x4155583850523031ULL,P::n,85,44,385,48,37,108,8,448} ||
        function != "sigmoid" || !finished) throw std::runtime_error("incompatible retained state");
    auto dense = std::make_unique<Key<P>>(); archive(*dense);
    for (const auto &value : *dense)
        if (value != T{} && value != T{1} && value != T{} - T{1})
            throw std::runtime_error("retained dense key is not ternary");
    const auto file_bytes = std::filesystem::file_size(state);
    auto skip = [&](std::uint64_t bytes) {
        const auto position = static_cast<std::uint64_t>(file.tellg());
        if (position > file_bytes || bytes > file_bytes - position) throw std::runtime_error("truncated retained state");
        file.seekg(static_cast<std::streamoff>(bytes), std::ios::cur);
        if (!file) throw std::runtime_error("cannot seek retained state");
    };
    // Skip sparse key, slots, input and saved output: fixed-array payloads.
    skip(sizeof(Key<P>) + (P::n / 2) * 2 * sizeof(double) + 2 * sizeof(TRLWE<P>));
    bool direct{}, conjugate{}; archive(direct, conjugate);
    for (unsigned branch = 0; branch < 2; ++branch) {
        std::array<std::uint64_t, 7> plan{}; int step{}; std::uint64_t groups{};
        archive(plan, step, groups);
        if (plan != std::array<std::uint64_t, 7>{0x434b4b53434c5031ULL,P::n,1,93,52,37,64} || groups > P::n / 2)
            throw std::runtime_error("unexpected compact S2C plan");
        std::uint64_t total_terms = 0;
        for (std::uint64_t g = 0; g < groups; ++g) {
            int giant{}; std::uint64_t terms{}; archive(giant, terms);
            if (terms > P::n / 2 || total_terms + terms > P::n / 2) throw std::runtime_error("unexpected S2C term count");
            total_terms += terms; skip(terms * (sizeof(int) + P::n * sizeof(std::int64_t)));
        }
    }
    skip(sizeof(CKKSGaloisKey<P, 93>) + sizeof(CKKSGaloisKey<P, 56>) + sizeof(CKKSSecretKeySwitchKey<P, 56>));
    std::uint64_t head{}, tail{}; archive(head, tail);
    if (head || tail) throw std::runtime_error("requires file-backed retained selectors");
    skip(P::nbit * sizeof(Raw));
    const auto raw_offset = static_cast<std::uint64_t>(file.tellg());
    auto raw = std::make_unique<Raw>(); archive(*raw);
    auto projected = std::make_unique<WideRaw>();
    constexpr auto first = CKKSKeySwitchFirstRowForLevel<P, 356>();
    constexpr auto wide_first = CKKSKeySwitchFirstRowForLevel<W, 356>();
    for (std::size_t r = 0; r < WideRaw::key_switch_rows; ++r) {
        const auto source = 2 * (wide_first + r + 1) - 1 - first;
        if (source >= Raw::key_switch_rows) throw std::runtime_error("projected row out of bounds");
        (*projected)[0][r] = (*raw)[0][source];
    }
    auto old_fft = makeCKKSPreparedAutoKeyFFT<P,356>(*raw);
    auto native_fft = makeCKKSPreparedAutoKeyFFT<P,356,16,8>(*raw);
    raw.reset();
    auto new_fft = makeCKKSPreparedAutoKeyFFT<W,356>(*projected);
    projected.reset();
    // Diagnostic ciphertexts with exact zero phase and public deterministic
    // masks under the retained secret. No secret/evaluation key is generated.
    std::mt19937_64 random(20262119);
    auto input = std::make_unique<Cipher>();
    auto old = std::make_unique<Cipher>();
    auto actual = std::make_unique<WideCipher>();
    auto native = std::make_unique<Cipher>();
    Key<W> wide_dense = *dense;
    std::array<std::complex<double>,P::n / 2> slots{};
    double old_error = 0, new_error = 0, max_difference = 0, input_error = 0;
    double old_seconds = 0, new_seconds = 0;
    const auto update = [](double &maximum, const auto &value) {
        const double error = std::abs(value);
        if (!std::isfinite(error)) throw std::runtime_error("nonfinite decoded phase");
        maximum = std::max(maximum, error);
    };
    std::ofstream outputs(destination, std::ios::binary);
    if (!outputs) throw std::runtime_error("cannot create public ciphertext diagnostic");
    cereal::PortableBinaryOutputArchive out(outputs);
    out(std::array<std::uint64_t,5>{Auxiliary8 ? 0x4155544f41383031ULL : 0x4155544f50383031ULL,P::n,q,delta,3});
    if constexpr (Auxiliary8) out(std::uint64_t{8});
    for (std::size_t sample = 0; sample < 3; ++sample) {
        CKKSSlotVector<P> target{};
        if constexpr (Auxiliary8) {
            std::uniform_real_distribution<double> distribution(-1, 1);
            for (auto &value : target) value = {distribution(random), distribution(random)};
            ckksSlotEncrypt<P>(*input, target, *dense, {P::α, 0});
            ckksSlotDecrypt<P,q,delta>(slots,*input,*dense);
            for (std::size_t i=0;i<slots.size();++i) update(input_error,slots[i]-target[i]);
            out(target);
        }
        else {
            for (auto &value : input->ct[0]) {
                for (auto &word : value.limb) word = random();
                value &= (T{1} << q) - T{1};
            }
            auto secret = std::make_unique<Polynomial<P>>();
            std::copy(dense->begin(), dense->end(), secret->begin());
            PolyMulTorusByDigit<P,q>(input->ct[1], input->ct[0], *secret);
            ckks_detail::reducePolynomialToLevel<P,q>(input->ct[1]);
        }
        // Keep input/output ciphertexts for an independent integer phase check.
        auto start = std::chrono::steady_clock::now();
        if constexpr (Auxiliary8)
            CKKSEvalAutoWithAuxiliaryScale<P,q,delta,8>(*old,*input,2 * P::n - 1,*old_fft);
        else CKKSEvalAuto<P,356>(old->ct,input->ct,2 * P::n - 1,*old_fft);
        old_seconds += std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        start = std::chrono::steady_clock::now();
        if constexpr (Auxiliary8) {
            auto wide_input = std::make_unique<WideCipher>(); wide_input->ct = input->ct;
            CKKSEvalAutoWithAuxiliaryScale<W,q,delta,8>(*actual,*wide_input,2 * P::n - 1,*new_fft);
        }
        else CKKSEvalAuto<W,356>(actual->ct,input->ct,2 * P::n - 1,*new_fft);
        new_seconds += std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        if constexpr (Auxiliary8)
            CKKSEvalAutoWithAuxiliaryScale<P,q,delta,8>(*native,*input,2 * P::n - 1,*native_fft);
        else CKKSEvalAuto<P,356>(native->ct,input->ct,2 * P::n - 1,*native_fft);
        if (native->ct != actual->ct) throw std::runtime_error("native primary projection differs from explicit retained-row projection");
        ckksSlotDecrypt<P,q,delta>(slots,*old,*dense);
        for (std::size_t i=0;i<slots.size();++i) update(old_error,slots[i]-std::conj(target[i]));
        ckksSlotDecrypt<W,q,delta>(slots,*actual,wide_dense);
        for (std::size_t i=0;i<slots.size();++i) update(new_error,slots[i]-std::conj(target[i]));
        auto difference = std::make_unique<WideCipher>();
        for (std::size_t c=0;c<2;++c)
            for (std::size_t i=0;i<P::n;++i) difference->ct[c][i]=actual->ct[c][i]-old->ct[c][i];
        ckksSlotDecrypt<W,q,delta>(slots,*difference,wide_dense);
        for (const auto &value : slots) update(max_difference,value);
        out(*input,*old,*actual);
    }
    outputs.flush(); if (!outputs) throw std::runtime_error("cannot flush diagnostic ciphertexts");
    const bool passed = std::isfinite(new_error) && new_error < std::ldexp(1.0,-23) && old_error < std::ldexp(1.0,-23);
    std::cout << std::setprecision(17) << "{\"n\":" << P::n << ",\"log_q\":" << q << ",\"log_delta\":" << delta << ",\"samples\":3"
        << ",\"auxiliary_bits\":" << (Auxiliary8 ? 8 : 0) << ",\"input_max_slot_error\":" << input_error
        << ",\"native_primary_projection_identical\":true"
        << ",\"raw_key_offset\":" << raw_offset << ",\"raw_key_bytes\":" << sizeof(Raw)
        << ",\"old_primary_rows\":" << Raw::key_switch_rows << ",\"new_primary_rows\":" << WideRaw::key_switch_rows
        << ",\"old_max_slot_error\":" << old_error << ",\"new_max_slot_error\":" << new_error
        << ",\"max_slot_difference\":" << max_difference << ",\"old_switch_seconds\":" << old_seconds
        << ",\"new_switch_seconds\":" << new_seconds << ",\"fresh_keys_generated\":false,\"accuracy_passed\":" << (passed ? "true" : "false") << "}\n";
    return passed ? 0 : 2;
}
} // namespace
int main(int argc,char **argv)
{
    try {
        if (argc != 4 && argc != 5) throw std::runtime_error("Usage: ckks_narrow_auto_primary_projection STATE N PUBLIC_OUTPUT [--auxiliary8]");
        if (argc == 5 && std::string(argv[4]) != "--auxiliary8") throw std::runtime_error("unknown mode");
        if (std::filesystem::exists(argv[3])) throw std::runtime_error("output must be new");
        if (argc == 5) {
            if (std::string(argv[2]) == "512") return run<9,true>(argv[1],argv[3]);
            if (std::string(argv[2]) == "16384") return run<14,true>(argv[1],argv[3]);
        }
        if (std::string(argv[2]) == "512") return run<9>(argv[1],argv[3]);
        if (std::string(argv[2]) == "16384") return run<14>(argv[1],argv[3]);
        throw std::runtime_error("unsupported ring");
    } catch (const std::exception &e) { std::cerr << e.what() << '\n'; return 2; }
}
