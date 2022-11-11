#ifdef ENABLE_LSHTEST
#include <gmp.h>
#include <gmpxx.h>
#endif

#include <cassert>
#include <cuhe++.hpp>
#include <iostream>
#include <random>

using namespace std;
using namespace cuHEpp;

int main()
{
    constexpr int numTest = 100000;
    random_device seed_gen;
    default_random_engine engine(seed_gen());
    uniform_int_distribution<uint64_t> dist(0, P);

    // Add Test
    for (int i = 0; i < numTest; i++) {
        __uint128_t a = dist(engine);
        __uint128_t b = dist(engine);
        INTorus A(a);
        INTorus B(b);
        assert((A + B).value == ((a + b) % P));
    }
    cout << "Add PASS" << endl;

    // Add and assignment Test
    for (int i = 0; i < numTest; i++) {
        __uint128_t a = dist(engine);
        __uint128_t b = dist(engine);
        INTorus A(a);
        INTorus B(b);
        A += B;
        assert(A.value == ((a + b) % P));
    }
    cout << "Add and assignment PASS" << endl;

    // Sub Test
    for (int i = 0; i < numTest; i++) {
        __uint128_t a = dist(engine);
        __uint128_t b = dist(engine);
        INTorus A(a);
        INTorus B(b);
        assert((A - B).value == ((a + (P - b)) % P));
    }
    cout << "Sub PASS" << endl;

    // Sub and assignment Test
    for (int i = 0; i < numTest; i++) {
        __uint128_t a = dist(engine);
        __uint128_t b = dist(engine);
        INTorus A(a);
        INTorus B(b);
        A -= B;
        // cout<<A.value<<":"<<static_cast<uint64_t>((a+(P-b))%P)<<endl;
        assert(A.value == ((a + (P - b)) % P));
    }
    cout << "Sub and assignment PASS" << endl;

    // Mul Test
    for (int i = 0; i < numTest; i++) {
        __uint128_t a = dist(engine);
        __uint128_t b = dist(engine);
        INTorus A(a);
        INTorus B(b);
        // cout<<(A*B).value<<":"<<static_cast<uint64_t>((a*b)%P)<<endl;
        assert((A * B).value == ((a * b) % P));
    }
    cout << "Mul PASS" << endl;

    // Mul and assignment Test
    for (int i = 0; i < numTest; i++) {
        __uint128_t a = dist(engine);
        __uint128_t b = dist(engine);
        INTorus A(a);
        INTorus B(b);
        A *= B;
        // cout<<A.value<<":"<<static_cast<uint64_t>((a*b)%P)<<endl;
        assert(A.value == ((a * b) % P));
    }
    cout << "Mul and assignment PASS" << endl;

#ifdef ENABLE_LSHTEST
    // Lsh Test
    for (int upper = 1; upper <= 6; upper++) {
        for (int l = 32 * (upper - 1); l < 32 * upper; l++) {
            // cout<<l<<endl;
            for (int i = 0; i < numTest; i++) {
                uint64_t temp = dist(engine);
                INTorus A(temp);
                mpz_class a = temp;
                // cout<<static_cast<uint64_t>(a)<<":";
                // if((A<<l).value!=(static_cast<uint64_t>((a<<l)%P))){
                if ((A << l).value != (a << l) % P) {
                    mpz_class res = (a << l) % P;
                    cout << (A << l).value << ":" << res.get_str() << endl;
                    std::cout << l << std::endl;
                    // cout<<(A<<l).value -
                    // static_cast<uint64_t>((a<<l)%P)<<endl;
                }
                // assert((A<<l).value==(static_cast<uint64_t>((a<<l)%P)));
                assert((A << l).value == (a << l) % P);
            }
            const __uint128_t temp = ((1UL << 32) - 1) << 32;
            INTorus A(temp);
            const mpz_class a = ((1UL << 32) - 1) << 32;
            mpz_class res = (a << l) % P;
            if ((A << l).value != ((a << l) % P)) {
                cout << (A << l).value << ":" << res.get_str() << endl;
                cout << "Here" << endl;
            }
            assert((A << l).value == ((a << l) % P));
        }
        cout << "Lsh" << upper * 32 << " PASS" << endl;
    }
#endif

    // InvPow2 Test
    for (int i = 1; i <= 31; i++) {
        assert((INTorus(1U << i, false) * cuHEpp::InvPow2(i)).value == 1);
    }
    cout << "InvPow2 PASS" << endl;
}