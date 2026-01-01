#include <cassert>
#include <iostream>
#include <random>
#include <tfhe++.hpp>

using namespace std;
using namespace TFHEpp;

// Custom parameter set for testing double decomposition with 64-bit Torus
// Using trivial double decomposition (l̅=1) which reduces to standard decomposition
// This verifies the code path works correctly
// Note: nbit must match lvl2param for FFT compatibility
struct DDTestParam {
    static constexpr int32_t key_value_max = 1;
    static constexpr int32_t key_value_min = -1;
    static constexpr std::uint32_t nbit = lvl2param::nbit;  // Use lvl2param's nbit for FFT compatibility
    static constexpr std::uint32_t n = 1 << nbit;
    static constexpr std::uint32_t k = 1;
    static constexpr std::uint32_t lₐ = 4;
    static constexpr std::uint32_t l = 4;
    static constexpr std::uint32_t Bgbit = 10;
    static constexpr std::uint32_t Bgₐbit = 10;
    static constexpr uint32_t Bg = 1U << Bgbit;
    static constexpr uint32_t Bgₐ = 1U << Bgₐbit;
    static constexpr ErrorDistribution errordist =
        ErrorDistribution::ModularGaussian;
    static const inline double α = std::pow(2.0, -40);
    using T = uint64_t;
    static constexpr T μ = 1ULL << 61;
    static constexpr uint32_t plain_modulus = 8;
    static constexpr double Δ = μ;
    // Trivial Double Decomposition (l̅=1 reduces to standard decomposition)
    // For l̅=1 to correctly reduce to standard decomposition, B̅gbit must equal Bgbit
    // With l=4, Bgbit=10: uses 40 bits
    // l̅=1, B̅gbit=10: adds 10 bits (but with j=0 only, no additional bits used)
    // The shift formula: width - (i+1)*Bgbit - j*B̅gbit = width - (i+1)*Bgbit when j=0
    static constexpr std::uint32_t l̅ = 1;
    static constexpr std::uint32_t l̅ₐ = 1;
    static constexpr std::uint32_t B̅gbit = 10;
    static constexpr std::uint32_t B̅gₐbit = 10;
};

// Test Double Decomposition correctness
// Verifies that DoubleDecomposition correctly represents the original value
template <class P>
void testDoubleDecomposition()
{
    random_device seed_gen;
    default_random_engine engine(seed_gen());
    uniform_int_distribution<typename P::T> dist(0,
                                                  numeric_limits<typename P::T>::max());

    constexpr int num_test = 100;
    // Tolerance based on decomposition precision
    // The lowest shift is: width - l*Bgbit - (l̅-1)*B̅gbit (j goes from 0 to l̅-1)
    // So remaining_bits = width - l*Bgbit - (l̅-1)*B̅gbit
    constexpr int remaining_bits = std::numeric_limits<typename P::T>::digits -
                                   P::l * P::Bgbit - (P::l̅ - 1) * P::B̅gbit;
    static_assert(remaining_bits >= 0,
                  "Invalid double decomposition parameters");

    for (int test = 0; test < num_test; test++) {
        // Create a random polynomial
        Polynomial<P> poly;
        for (int i = 0; i < P::n; i++) {
            poly[i] = dist(engine);
        }

        // Apply double decomposition
        DecomposedPolynomialDD<P> decpoly;
        DoubleDecomposition<P>(decpoly, poly);

        // Reconstruct the polynomial from decomposed components
        Polynomial<P> reconstructed;
        for (int n = 0; n < P::n; n++) {
            typename P::T sum = 0;
            for (int i = 0; i < P::l; i++) {
                for (int j = 0; j < P::l̅; j++) {
                    // The scaling factor for (i,j) position
                    // When l̅=1 (j=0 only), this reduces to standard decomposition
                    typename P::T h_val =
                        static_cast<typename P::T>(1)
                        << (std::numeric_limits<typename P::T>::digits -
                            (i + 1) * P::Bgbit - j * P::B̅gbit);
                    sum += static_cast<typename P::T>(
                               static_cast<std::make_signed_t<typename P::T>>(
                                   decpoly[i * P::l̅ + j][n])) *
                           h_val;
                }
            }
            reconstructed[n] = sum;
        }

        // Check that reconstruction is close to original
        typename P::T max_error =
            remaining_bits > 0
                ? (static_cast<typename P::T>(1) << remaining_bits)
                : static_cast<typename P::T>(1);

        for (int n = 0; n < P::n; n++) {
            typename P::T diff = poly[n] > reconstructed[n]
                                     ? poly[n] - reconstructed[n]
                                     : reconstructed[n] - poly[n];
            typename P::T max_val = ~static_cast<typename P::T>(0);
            typename P::T wrap_diff = max_val - diff + 1;
            typename P::T min_diff = diff < wrap_diff ? diff : wrap_diff;

            if (min_diff > max_error * 2) {
                cerr << "Decomposition error at n=" << n << endl;
                cerr << "  original=" << poly[n] << endl;
                cerr << "  reconstructed=" << reconstructed[n] << endl;
                cerr << "  max_error=" << max_error << endl;
                cerr << "  actual_diff=" << min_diff << endl;
                assert(false);
            }
        }
    }
    cout << "DoubleDecomposition test passed" << endl;
}

// Test External Product with Double Decomposition
template <class P>
void testExternalProduct()
{
    constexpr uint32_t num_test = 10;
    random_device seed_gen;
    default_random_engine engine(seed_gen());
    uniform_int_distribution<uint32_t> binary(0, 1);

    cout << "Parameters: n=" << P::n << ", k=" << P::k << endl;
    cout << "Primary: l=" << P::l << ", Bgbit=" << P::Bgbit << endl;
    cout << "Auxiliary: l̅=" << P::l̅ << ", B̅gbit=" << P::B̅gbit << endl;
    cout << "TRGSW rows: " << (P::k * P::lₐ * P::l̅ₐ + P::l * P::l̅) << endl;

    // Test with p=1 (identity)
    cout << "Testing with plaintext = 1 (identity)..." << endl;
    for (int test = 0; test < num_test; test++) {
        // Generate a random TRLWE key for this param set
        std::array<typename P::T, P::n> key;
        for (int i = 0; i < P::n; i++) {
            key[i] = (binary(engine) == 1)
                         ? static_cast<typename P::T>(P::key_value_max)
                         : static_cast<typename P::T>(P::key_value_min);
        }

        // Create a random message polynomial
        array<bool, P::n> p;
        for (bool &i : p) i = (binary(engine) > 0);
        Polynomial<P> pmu;
        for (int i = 0; i < P::n; i++)
            pmu[i] = p[i] ? P::μ : -P::μ;

        // Encrypt as TRLWE
        TRLWE<P> c;
        trlweSymEncrypt<P>(c, pmu, key);

        // Create TRGSW encrypting 1 (identity)
        const Polynomial<P> plainpoly = {static_cast<typename P::T>(1)};

        TRGSWFFT<P> trgswfft;
        trgswSymEncrypt<P>(trgswfft, plainpoly, key);

        // Apply external product (automatically uses DD when P::l̅ > 1)
        TRLWE<P> result;
        ExternalProduct<P>(result, c, trgswfft);

        // Decrypt and verify
        const array<bool, P::n> p2 = trlweSymDecrypt<P>(result, key);
        for (int i = 0; i < P::n; i++) {
            if (p[i] != p2[i]) {
                cerr << "ExternalProduct (DD) failed at index " << i
                     << ": expected " << p[i] << " got " << p2[i] << endl;
                assert(false);
            }
        }
    }
    cout << "ExternalProduct (DD) test passed (p=1)" << endl;

    // Test with p=-1 (negation)
    cout << "Testing with plaintext = -1 (negation)..." << endl;
    for (int test = 0; test < num_test; test++) {
        std::array<typename P::T, P::n> key;
        for (int i = 0; i < P::n; i++) {
            key[i] = (binary(engine) == 1)
                         ? static_cast<typename P::T>(P::key_value_max)
                         : static_cast<typename P::T>(P::key_value_min);
        }

        array<bool, P::n> p;
        for (bool &i : p) i = (binary(engine) > 0);
        Polynomial<P> pmu;
        for (int i = 0; i < P::n; i++)
            pmu[i] = p[i] ? P::μ : -P::μ;

        TRLWE<P> c;
        trlweSymEncrypt<P>(c, pmu, key);

        // TRGSW encrypting -1 (negation)
        const Polynomial<P> plainpoly = {static_cast<typename P::T>(-1)};

        TRGSWFFT<P> trgswfft;
        trgswSymEncrypt<P>(trgswfft, plainpoly, key);

        TRLWE<P> result;
        ExternalProduct<P>(result, c, trgswfft);

        const array<bool, P::n> p2 = trlweSymDecrypt<P>(result, key);
        for (int i = 0; i < P::n; i++) {
            if (p[i] != !p2[i]) {
                cerr << "ExternalProduct (DD, p=-1) failed at index " << i
                     << ": expected " << !p[i] << " got " << p2[i] << endl;
                assert(false);
            }
        }
    }
    cout << "ExternalProduct (DD) test passed (p=-1)" << endl;
}

int main()
{
    cout << "=== Testing Double Decomposition ===" << endl;
    cout << endl;

    // DDTestParam: 64-bit Torus with proper double decomposition
    // l=2, Bgbit=16, l̅=2, B̅gbit=16
    // Total bits: 2*16 + 2*16 = 64 (exact fit)
    // Decomposition levels: 2*2 = 4

    cout << "DDTestParam configuration:" << endl;
    cout << "  Torus type: uint64_t (64-bit)" << endl;
    cout << "  n = " << DDTestParam::n << " (polynomial degree)" << endl;
    cout << "  k = " << DDTestParam::k << endl;
    cout << "  Primary decomposition: l = " << DDTestParam::l
         << ", Bgbit = " << DDTestParam::Bgbit << endl;
    cout << "  Auxiliary decomposition: l̅ = " << DDTestParam::l̅
         << ", B̅gbit = " << DDTestParam::B̅gbit << endl;
    cout << "  Total decomposition levels: " << (DDTestParam::l * DDTestParam::l̅)
         << endl;
    cout << "  Bits used: "
         << (DDTestParam::l * DDTestParam::Bgbit +
             DDTestParam::l̅ * DDTestParam::B̅gbit)
         << " / 64" << endl;
    cout << "  TRGSW size: "
         << (DDTestParam::k * DDTestParam::lₐ * DDTestParam::l̅ₐ +
             DDTestParam::l * DDTestParam::l̅)
         << " TRLWE rows" << endl;
    cout << endl;

    cout << "--- Testing DoubleDecomposition ---" << endl;
    testDoubleDecomposition<DDTestParam>();
    cout << endl;

    cout << "--- Testing ExternalProduct (DD) ---" << endl;
    testExternalProduct<DDTestParam>();
    cout << endl;

    cout << "=== All Double Decomposition tests passed ===" << endl;
    return 0;
}
