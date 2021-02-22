#include <cassert>
#include <chrono>
#include <iostream>
#include <random>
#include <tfhe++.hpp>

using namespace std;
using namespace TFHEpp;

int main()
{
    const uint32_t num_test = 10;
    random_device seed_gen;
    default_random_engine engine(seed_gen());
    uniform_int_distribution<uint32_t> binary(0, 1);

    SecretKey *sk = new SecretKey;
    CircuitKeylvl22 *ck = new CircuitKeylvl22(*sk);
    vector<array<uint8_t, TFHEpp::DEF_nbar>> pa(num_test);
    vector<array<uint64_t, TFHEpp::DEF_nbar>> pmu(num_test);
    vector<uint8_t> pones(num_test);
    array<bool, TFHEpp::DEF_nbar> pres;
    for (array<uint8_t, TFHEpp::DEF_nbar> &i : pa)
        for (uint8_t &p : i) p = binary(engine);
    for (int i = 0; i < num_test; i++)
        for (int j = 0; j < TFHEpp::DEF_nbar; j++)
            pmu[i][j] = pa[i][j] ? DEF_μbar : -DEF_μbar;
    for (int i = 0; i < num_test; i++) pones[i] = true;
    vector<TRLWElvl2> ca(num_test);
    vector<TLWElvl0> cones(num_test);
    vector<TRGSWFFTlvl2> bootedTGSW(num_test);

    for (int i = 0; i < num_test; i++)
        ca[i] = trlweSymEncryptlvl2(pmu[i], TFHEpp::DEF_αbklvl02, sk->key.lvl2);
    cones = bootsSymEncrypt(pones, *sk);

    chrono::system_clock::time_point start, end;
    start = chrono::system_clock::now();
    for (int test = 0; test < num_test; test++) {
        CircuitBootstrappingFFTlvl22(bootedTGSW[test], cones[test], *ck);
    }
    end = chrono::system_clock::now();

    for (int test = 0; test < num_test; test++) {
        trgswfftExternalProductlvl2(ca[test], ca[test], bootedTGSW[test]);
        pres = trlweSymDecryptlvl2(ca[test], sk->key.lvl2);
        for (int i = 0; i < TFHEpp::DEF_nbar; i++)
            assert(pres[i] == pa[test][i]);
    }
    cout << "Passed" << endl;
    double elapsed =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
            .count();
    cout << elapsed / num_test << "ms" << endl;
}
