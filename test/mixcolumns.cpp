#include <cassert>
#include <chrono>
#include <iostream>
#include <memory>
#include <random>
#include <tfhe++.hpp>

void MixColumns(unsigned char state[4][4]) {
  unsigned char temp_state[4][4];

  for (size_t i = 0; i < 4; ++i) {
    memset(temp_state[i], 0, 4);
  }

  for (size_t i = 0; i < 4; ++i) {
    for (size_t k = 0; k < 4; ++k) {
      for (size_t j = 0; j < 4; ++j) {
        if (CMDS[i][k] == 1)
          temp_state[i][j] ^= state[k][j];
        else
          temp_state[i][j] ^= GF_MUL_TABLE[CMDS[i][k]][state[k][j]];
      }
    }
  }

  for (size_t i = 0; i < 4; ++i) {
    memcpy(state[i], temp_state[i], 4);
  }
}

int main()
{
    using P = TFHEpp::lvl2param;
    std::random_device seed_gen;
    std::default_random_engine engine(seed_gen());
    std::uniform_int_distribution<uint32_t> binary(0, 1);

    std::unique_ptr<TFHEpp::SecretKey> sk(new TFHEpp::SecretKey());
    constexpr uint num_test = 1000;
    std::vector<std::array<TFHEpp::TLWE<P>,128>> cstate(num_test);
    std::vector<std::array<uint8_t, 128>> plaintext(num_test);

    for (int i = 0; i < num_test; i++) {
        for (int j = 0; j < 128; j++){
            plaintext[i][j] = binary(engine);
            cstate[i][j] = TFHEpp::tlweSymEncrypt<P>(plaintext[i][j]?1ULL << (std::numeric_limits<typename P::T>::digits - 2):-(1ULL << (std::numeric_limits<typename P::T>::digits - 2)), *sk);
        }
    }

    std::chrono::system_clock::time_point start, end;
    start = std::chrono::system_clock::now();
    for (int test = 0; test < num_test; test++) {
        // std::cout << "test: " << test << std::endl;
        TFHEpp::MixColumns<P>(cstate[test]);
    }

    end = std::chrono::system_clock::now();
    double elapsed =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
            .count();
    std::cout << elapsed / num_test << "ms" << std::endl;

    for (int i = 0; i < num_test; i++) {
        std::array<uint8_t, 128> pres;
        for (int j = 0; j < 128; j++)
            pres[j] = TFHEpp::tlweSymDecrypt<P>(cstate[i][j], *sk);
        unsigned char state[4][4];
        for (int j = 0; j < 4; j++)
            for (int k = 0; k < 4; k++){
                uint8_t byte = 0;
                for (int l = 0; l < 8; l++)
                    byte |= plaintext[i][j*32 + k * 8 + l] << l;
                state[j][k] = byte;
            }
        MixColumns(state);
        for (int j = 0; j < 4; j++)
            for (int k = 0; k < 4; k++){
                uint8_t byte = 0;
                for (int l = 0; l < 8; l++)
                    byte |= pres[j*32 + k * 8 + l] << l;
                // std::cout <<"j: " << j << " k: " << k << std::endl;
                // std::cout << (int)state[j][k] << " " << (int)byte << std::endl;
                assert(state[j][k] == byte);
            }
    }
  std::cout << "PASS" << std::endl;
}