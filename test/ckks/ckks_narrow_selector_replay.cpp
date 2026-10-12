#include <ckks/narrow_selector_file.hpp>
#include <chrono>
#include <fstream>
#include <iterator>
#include <iostream>

namespace {
void require(bool value, const char *message) { if (!value) throw std::runtime_error(message); }
struct P : TFHEpp::lvl5param { static constexpr std::uint32_t nbit = 4, n = 16; };
// Only the file source's layout fields are exercised here. Complete plain37
// checkpoint/replay is separately checked by the encrypted benchmark.
struct Schedule {
    using Param = P;
    using Selector = TFHEpp::CKKSCiphertext<P, 53, 23>;
    static constexpr std::size_t blocks = 2, tail_packs = 1, block_width = 3;
    static constexpr std::uint32_t selector_log_q = 53, log_delta = 23, plain_log_delta = 7;
};
using File = TFHEpp::CKKSNarrowSelectorFile<Schedule>;
}
int main()
{
    const auto directory = std::filesystem::temp_directory_path() /
        ("ckks-selector-replay-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    std::vector<Schedule::Selector> input(Schedule::block_width), restored;
    for (std::size_t col = 0; col < input.size(); col++)
        for (std::size_t component = 0; component < 2; component++)
            for (std::size_t point = 0; point < P::n; point++)
                for (std::size_t limb = 0; limb < P::T::limbs; limb++)
                    input[col].ct[component][point].limb[limb] = 1000 * col + 100 * component + 10 * point + limb;
    {
        File file(directory);
        bool rejected = false;
        try { file.preserve(); } catch (const std::logic_error &) { rejected = true; }
        require(rejected, "incomplete cache was retained");
        for (bool tail : {false, true})
            for (std::size_t group = 0; group < (tail ? Schedule::tail_packs : Schedule::blocks); group++)
                file.write(tail, group, input);
        file.finish(); file.preserve();
    }
    require(std::filesystem::exists(directory / "head-0.bin"), "completed cache was not retained");
    {
        File reopened(directory, File::Existing{});
        for (bool tail : {false, true})
            for (std::size_t group = 0; group < (tail ? Schedule::tail_packs : Schedule::blocks); group++) {
                reopened.load(tail, group, restored);
                require(restored.size() == input.size(), "reopened selector count changed");
                for (std::size_t i = 0; i < input.size(); i++)
                    require(restored[i].ct == input[i].ct, "reopened native selector words changed");
            }
        bool rejected = false;
        try { reopened.write(false, 0, input); } catch (const std::invalid_argument &) { rejected = true; }
        require(rejected, "reopened source accepted a write");
    }
    require(std::filesystem::exists(directory / "tail-0.bin"), "reader deleted retained files");
    const auto file = directory / "head-0.bin";
    std::ifstream stream(file, std::ios::binary);
    const std::string original{std::istreambuf_iterator<char>(stream), std::istreambuf_iterator<char>()};
    stream.close();
    for (int kind = 0; kind < 4; kind++) {
        auto changed = original;
        if (kind == 0) changed[9] ^= 1; // N metadata, after endianness and magic
        if (kind == 1) changed.pop_back();
        if (kind == 2) changed.push_back('x');
        if (kind == 3) std::filesystem::rename(file, directory / "saved.bin");
        else { std::ofstream out(file, std::ios::binary | std::ios::trunc); out.write(changed.data(), changed.size()); }
        bool rejected = false;
        try { File bad(directory, File::Existing{}); } catch (const std::exception &) { rejected = true; }
        require(rejected, "malformed retained selector source accepted");
        require(std::filesystem::exists(directory / "tail-0.bin"), "failed reader deleted retained files");
        if (kind == 3) std::filesystem::rename(directory / "saved.bin", file);
        else { std::ofstream out(file, std::ios::binary | std::ios::trunc); out.write(original.data(), original.size()); }
    }
    bool rejected = false;
    try { File duplicate(directory); } catch (const std::runtime_error &) { rejected = true; }
    require(rejected && std::filesystem::exists(file), "new-cache constructor overwrote retained source");
    for (const auto &entry : std::filesystem::directory_iterator(directory)) std::filesystem::remove(entry.path());
    std::filesystem::remove(directory);
    {
        File owned(directory);
        owned.write(false, 0, input);
    }
    require(!std::filesystem::exists(directory), "normal owned source cleanup changed");
    std::cout << "selector retain/reopen checks passed: three exact group roundtrips, four malformed source rejections, read-only lifetime and normal cleanup\n";
}
