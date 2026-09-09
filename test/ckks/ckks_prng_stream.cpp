#include <array>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <vector>
#include <BLAKE3PRNG.hpp>

template <class Word>
bool check_stream()
{
    std::array<std::uint8_t, 32> seed{};
    for (std::size_t i = 0; i < seed.size(); i++) seed[i] = static_cast<std::uint8_t>(i);
    BLAKE3PRNG::BLAKE3PRNG<Word> stream(seed);
    blake3_hasher reference;
    blake3_hasher_init_keyed(&reference, seed.data());
    // Cross multiple refill boundaries; an incorrect byte/word stride can
    // pass a test that reads only the initial buffer.
    std::vector<std::uint8_t> expected(4*4096);
    blake3_hasher_finalize_seek(&reference, 0, expected.data(), expected.size());
    for (std::size_t offset = 0; offset < expected.size(); offset += sizeof(Word)) {
        Word word;
        std::memcpy(&word, expected.data()+offset, sizeof(Word));
        if (stream() != word) {
            std::cerr << "BLAKE3 stream mismatch at byte " << offset
                      << " for " << 8*sizeof(Word) << "-bit words\n";
            return false;
        }
    }
    return true;
}

int main()
{
    if (!check_stream<std::uint64_t>() || !check_stream<std::uint32_t>()) return 1;
    std::cout << "Passed: four BLAKE3 buffers match the reference XOF\n";
}
