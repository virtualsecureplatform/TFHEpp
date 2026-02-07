#include <cassert>
#include <chrono>
#include <iostream>
#include <random>
#include <tfhe++.hpp>

int main()
{
#if defined(USE_CONCRETE) || defined(USE_CONCRETE_FFT)
    // Skip this test for CONCRETE builds - lvl3param uses nbit=13 which exceeds
    // the FFT table sizes initialized for CONCRETE parameters
    std::cout << "Skipping DD gate bootstrapping test for CONCRETE build" << std::endl;
    return 0;
#else
    using namespace std;
    using namespace TFHEpp;

    constexpr uint32_t num_test = 10;  // Reduced for faster testing
    random_device seed_gen;
    default_random_engine engine(seed_gen());
    uniform_int_distribution<uint32_t> binary(0, 1);

    // Use lvl03param which bootstraps from lvl0param to lvl3param
    // lvl3param has nbit=12 (n=4096) and non-trivial DD: l=2, l̅=2, Bgbit=16, B̅gbit=16
    using bkP = lvl03param;

    cout << "=== Testing GateBootstrappingTLWE2TLWE with DD (lvl3param, nbit=12) ===" << endl;
    cout << "Domain: n=" << bkP::domainP::n << ", T=" << sizeof(typename bkP::domainP::T) * 8 << "-bit" << endl;
    cout << "Target: n=" << bkP::targetP::n << ", nbit=" << bkP::targetP::nbit << ", T=" << sizeof(typename bkP::targetP::T) * 8 << "-bit" << endl;
    cout << "Primary decomposition: l=" << bkP::targetP::l << ", Bgbit=" << bkP::targetP::Bgbit << endl;
    cout << "Auxiliary decomposition: l̅=" << bkP::targetP::l̅ << ", B̅gbit=" << bkP::targetP::B̅gbit << endl;
    cout << "Total decomposition levels: " << (bkP::targetP::l * bkP::targetP::l̅) << endl;
    cout << "Bits used: " << (bkP::targetP::l * bkP::targetP::Bgbit + bkP::targetP::l̅ * bkP::targetP::B̅gbit) << " / " << std::numeric_limits<typename bkP::targetP::T>::digits << endl;
    cout << endl;

    // Generate keys
    cout << "Generating secret keys..." << endl;
    array<typename bkP::domainP::T, bkP::domainP::n> domainKey;
    for (int i = 0; i < bkP::domainP::n; i++) {
        domainKey[i] = binary(engine);
    }

    array<typename bkP::targetP::T, bkP::targetP::n> targetKey;
    for (int i = 0; i < bkP::targetP::n; i++) {
        targetKey[i] = (binary(engine) == 1)
                           ? static_cast<typename bkP::targetP::T>(bkP::targetP::key_value_max)
                           : static_cast<typename bkP::targetP::T>(bkP::targetP::key_value_min);
    }

    // Generate bootstrapping key
    cout << "Generating bootstrapping key (this may take a while for n=4096)..." << endl;
    auto bkfft = make_unique<BootstrappingKeyFFT<bkP>>();
    bkfftgen<bkP>(*bkfft, domainKey, targetKey);

    // Test arrays - use heap allocation for large 128-bit TLWEs to avoid stack overflow
    auto tlwe = make_unique<array<TLWE<typename bkP::domainP>, num_test>>();
    auto bootedtlwe = make_unique<array<TLWE<typename bkP::targetP>, num_test>>();
    array<bool, num_test> p;

    // Encrypt test values
    cout << "Encrypting test values..." << endl;
    for (int i = 0; i < num_test; i++) {
        p[i] = binary(engine) > 0;
        tlweSymEncrypt<typename bkP::domainP>(
            (*tlwe)[i], p[i] ? bkP::domainP::μ : -bkP::domainP::μ, bkP::domainP::α,
            domainKey);
    }

    // Perform gate bootstrapping (automatically uses DD when l̅ > 1)
    cout << "Running GateBootstrappingTLWE2TLWE (with DD)..." << endl;
    chrono::system_clock::time_point start, end;
    start = chrono::system_clock::now();

    for (int test = 0; test < num_test; test++) {
        GateBootstrappingTLWE2TLWE<bkP>(
            (*bootedtlwe)[test], (*tlwe)[test], *bkfft,
            μpolygen<typename bkP::targetP, bkP::targetP::μ>());
    }

    end = chrono::system_clock::now();

    // Verify results
    cout << "Verifying results..." << endl;
    int errors = 0;
    for (int i = 0; i < num_test; i++) {
        bool p2 = tlweSymDecrypt<typename bkP::targetP>((*bootedtlwe)[i], targetKey);
        if (p[i] != p2) {
            cerr << "Error at index " << i << ": expected " << p[i] << " got " << p2 << endl;
            errors++;
        }
    }

    if (errors == 0) {
        cout << "All tests passed!" << endl;
    }
    else {
        cerr << errors << " out of " << num_test << " tests failed!" << endl;
        return 1;
    }

    double elapsed =
        chrono::duration_cast<chrono::milliseconds>(end - start).count();
    cout << "Average time per bootstrapping: " << elapsed / num_test << "ms" << endl;

    cout << "=== GateBootstrappingTLWE2TLWE (DD) test completed ===" << endl;
    return 0;
#endif
}
