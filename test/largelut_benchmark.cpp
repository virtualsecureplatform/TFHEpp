#include <array>
#include <chrono>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numeric>
#include <tfhe++.hpp>
#include <vector>

namespace {

using Clock = std::chrono::steady_clock;

template <class F>
double time_ms(F &&f)
{
    const auto start = Clock::now();
    f();
    const auto end = Clock::now();
    return std::chrono::duration<double, std::milli>(end - start).count();
}

double average(const std::vector<double> &values)
{
    return std::accumulate(values.begin(), values.end(), 0.0) / values.size();
}

template <class P>
uint64_t checksum(const TFHEpp::TRLWE<P> &trlwe)
{
    uint64_t acc = 0;
    for (const auto &poly : trlwe)
        for (const auto coeff : poly) acc += static_cast<uint64_t>(coeff);
    return acc;
}

}  // namespace

int main()
{
    using P = TFHEpp::lvl1param;
    using r2rP = TFHEpp::lvl10param;

    constexpr uint32_t W = 2;
    constexpr uint32_t K = 3;
    constexpr uint32_t step = 1;
    constexpr uint32_t repeats = 5;
    constexpr uint32_t r2r_t = r2rP::t;
    constexpr uint32_t r2r_basebit = r2rP::basebit;
    constexpr uint32_t table_tlwes = 4;

    auto sk = std::make_unique<TFHEpp::SecretKey>();

    auto ahk = std::make_unique<TFHEpp::AnnihilateKey<P>>();
    const double ah_keygen_ms =
        time_ms([&] { TFHEpp::annihilatekeygen<P>(*ahk, *sk); });

    const auto positions = TFHEpp::LargeLUTStepR2RPositions<P, W, K>(step);
    const uint32_t R = TFHEpp::LargeLUTStepBlockSize<P, W, K>(step);

    auto r2rk = std::make_unique<TFHEpp::R2RKey<P, r2r_t, r2r_basebit>>();
    const double r2r_keygen_ms = time_ms([&] {
        TFHEpp::R2RKeyGen<P, r2r_t, r2r_basebit>(
            *r2rk, positions, R, sk->key.get<P>());
    });

    TFHEpp::TRLWE<P> input;
    TFHEpp::trlweSymEncryptZero<P>(input, sk->key.get<P>());

    auto extract_and_duplicate = [&] {
        std::vector<TFHEpp::TLWE<P>> packed;
        packed.reserve(positions.size() * static_cast<std::size_t>(R));
        for (const uint32_t position : positions) {
            TFHEpp::TLWE<P> extracted;
            TFHEpp::SampleExtractIndex<P>(extracted, input, position);
            for (uint32_t r = 0; r < R; r++) packed.push_back(extracted);
        }
        return packed;
    };

    auto extract_positions_array = [&] {
        std::array<TFHEpp::TLWE<P>, table_tlwes> extracted_positions{};
        for (uint32_t p = 0; p < positions.size(); p++)
            TFHEpp::SampleExtractIndex<P>(extracted_positions[p], input,
                                          positions[p]);
        return extracted_positions;
    };

    auto packed = extract_and_duplicate();
    auto extracted_positions = extract_positions_array();

    TFHEpp::TRLWE<P> ah_out;
    TFHEpp::TRLWE<P> ah_specialized_out;
    TFHEpp::TRLWE<P> r2r_out;

    TFHEpp::TLWE2TRLWEPacking<P>(ah_out, packed, *ahk);
    TFHEpp::TLWE2TablePacking<P, table_tlwes>(ah_specialized_out,
                                              extracted_positions, *ahk);
    TFHEpp::R2RPKS<P, r2r_t, r2r_basebit>(r2r_out, input, positions, R,
                                          *r2rk);

    std::vector<double> extract_ms;
    std::vector<double> extract_positions_ms;
    std::vector<double> ah_pack_ms;
    std::vector<double> ah_total_ms;
    std::vector<double> ah_specialized_pack_ms;
    std::vector<double> ah_specialized_total_ms;
    std::vector<double> r2r_ms;
    uint64_t sink =
        checksum<P>(ah_out) + checksum<P>(ah_specialized_out) + checksum<P>(r2r_out);

    for (uint32_t repeat = 0; repeat < repeats; repeat++) {
        std::vector<TFHEpp::TLWE<P>> local_packed;
        std::array<TFHEpp::TLWE<P>, table_tlwes> local_extracted_positions;

        extract_ms.push_back(time_ms([&] { local_packed = extract_and_duplicate(); }));
        extract_positions_ms.push_back(
            time_ms([&] { local_extracted_positions = extract_positions_array(); }));

        ah_pack_ms.push_back(time_ms([&] {
            TFHEpp::TLWE2TRLWEPacking<P>(ah_out, packed, *ahk);
            sink += checksum<P>(ah_out);
        }));

        ah_total_ms.push_back(time_ms([&] {
            local_packed = extract_and_duplicate();
            TFHEpp::TLWE2TRLWEPacking<P>(ah_out, local_packed, *ahk);
            sink += checksum<P>(ah_out);
        }));

        ah_specialized_pack_ms.push_back(time_ms([&] {
            auto local = extracted_positions;
            TFHEpp::TLWE2TablePacking<P, table_tlwes>(ah_specialized_out, local,
                                                      *ahk);
            sink += checksum<P>(ah_specialized_out);
        }));

        ah_specialized_total_ms.push_back(time_ms([&] {
            local_extracted_positions = extract_positions_array();
            TFHEpp::TLWE2TablePacking<P, table_tlwes>(
                ah_specialized_out, local_extracted_positions, *ahk);
            sink += checksum<P>(ah_specialized_out);
        }));

        r2r_ms.push_back(time_ms([&] {
            TFHEpp::R2RPKS<P, r2r_t, r2r_basebit>(r2r_out, input, positions, R,
                                                  *r2rk);
            sink += checksum<P>(r2r_out);
        }));
    }

    const std::size_t trlwe_bytes = sizeof(TFHEpp::TRLWE<P>);
    const std::size_t r2rk_bytes =
        P::n * r2r_t * (std::size_t{1} << (r2r_basebit - 1)) * trlwe_bytes;
    const std::size_t packed_bytes = packed.size() * sizeof(TFHEpp::TLWE<P>);

    std::cout << std::fixed << std::setprecision(3);
    std::cout << "LargeLUT packing benchmark, P=lvl1param, W=" << W
              << ", K=" << K << ", step=" << step << "\n";
    std::cout << "positions=" << positions.size() << ", R=" << R
              << ", packed TLWEs=" << packed.size() << "\n";
    std::cout << "AnnihilateKey size MiB: "
              << static_cast<double>(sizeof(TFHEpp::AnnihilateKey<P>)) /
                     (1024.0 * 1024.0)
              << "\n";
    std::cout << "R2RKey size MiB: "
              << static_cast<double>(r2rk_bytes) / (1024.0 * 1024.0)
              << "\n";
    std::cout << "Intermediate TLWE vector MiB: "
              << static_cast<double>(packed_bytes) / (1024.0 * 1024.0)
              << "\n";
    std::cout << "AnnihilateKey gen ms: " << ah_keygen_ms << "\n";
    std::cout << "R2RKey gen ms: " << r2r_keygen_ms << "\n";
    std::cout << "SampleExtract+duplicate ms: " << average(extract_ms) << "\n";
    std::cout << "SampleExtract positions only ms: "
              << average(extract_positions_ms) << "\n";
    std::cout << "TLWE2TRLWEPacking only ms: " << average(ah_pack_ms) << "\n";
    std::cout << "SampleExtract+TLWE2TRLWEPacking ms: " << average(ah_total_ms)
              << "\n";
    std::cout << "Specialized AH table packing only ms: "
              << average(ah_specialized_pack_ms) << "\n";
    std::cout << "SampleExtract+specialized AH table packing ms: "
              << average(ah_specialized_total_ms) << "\n";
    std::cout << "R2RPKS ms: " << average(r2r_ms) << "\n";
    std::cout << "runtime speedup vs total AH path: "
              << average(ah_total_ms) / average(r2r_ms) << "x\n";
    std::cout << "runtime speedup vs specialized AH path: "
              << average(ah_specialized_total_ms) / average(r2r_ms) << "x\n";
    std::cout << "checksum: " << sink << "\n";
}
