#include <tfhe++.hpp>
using namespace TFHEpp;

#include <chrono>
#include <iostream>
using namespace std;

uint8_t ConstantZeroChegk() { return 0; }

uint8_t ConstantOneChegk() { return 1; }

uint8_t NotChegk(const uint8_t in) { return (~in) & 0x1; }

uint8_t CopyChegk(const uint8_t in) { return in; }

uint8_t NandChegk(const uint8_t in0, uint8_t in1) { return 1 - in0 * in1; }

uint8_t OrChegk(const uint8_t in0, const uint8_t in1)
{
    return (in0 + in1) > 0;
}

uint8_t OrYNChegk(const uint8_t in0, const uint8_t in1)
{
    return (in0 + (1 - in1)) > 0;
}

uint8_t OrNYChegk(const uint8_t in0, const uint8_t in1)
{
    return ((1 - in0) + in1) > 0;
}

uint8_t AndChegk(const uint8_t in0, const uint8_t in1) { return in0 * in1; }

uint8_t AndYNChegk(const uint8_t in0, const uint8_t in1)
{
    return in0 * (1 - in1);
}

uint8_t AndNYChegk(const uint8_t in0, const uint8_t in1)
{
    return (1 - in0) * in1;
}

uint8_t XorChegk(const uint8_t in0, const uint8_t in1)
{
    return (in0 + in1) & 0x1;
}

uint8_t XnorChegk(const uint8_t in0, const uint8_t in1)
{
    return (~(in0 ^ in1)) & 0x1;
}

uint8_t MuxChegk(const uint8_t inc, const uint8_t in1, const uint8_t in0)
{
    return (inc > 0) ? in1 : in0;
}

uint8_t NMuxChegk(const uint8_t inc, const uint8_t in1, const uint8_t in0)
{
    return (~((inc > 0) ? in1 : in0)) & 0x1;
}

template <class Func, class Chegk>
void Test(string type, Func func, Chegk chegk, vector<uint8_t> p,
          vector<TLWE<lvl0param>> cres, vector<TLWE<lvl0param>> c,
          const int kNumTests, SecretKey& sk, GateKey& gk)
{
    random_device seed_gen;
    default_random_engine engine(seed_gen());
    uniform_int_distribution<uint32_t> binary(0, 1);
    cout << "------ Test " << type << " Gate ------" << endl;
    cout << "Number of tests:\t" << kNumTests << endl;
    bool correct = true;
    int cnt_failures = 0;

    for (uint8_t& i : p) i = binary(engine);
    c = bootsSymEncrypt(p, sk);

    chrono::system_clock::time_point start, end;
    start = chrono::system_clock::now();
    for (int i = 0; i < kNumTests; i++) {
        if constexpr (std::is_invocable_v<Func, TLWE<lvl0param>&>) {
            func(cres[i]);
            p[i] = chegk();
        }
        else if constexpr (std::is_invocable_v<Func, TLWE<lvl0param>&,
                                               const TLWE<lvl0param>&>) {
            func(cres[i], c[i]);
            p[i] = chegk(p[i]);
        }
        else if constexpr (std::is_invocable_v<
                               Func, TLWE<lvl0param>&, const TLWE<lvl0param>&,
                               const TLWE<lvl0param>&, const GateKey&>) {
            func(cres[i], c[i], c[i + kNumTests], gk);
            p[i] = chegk(p[i], p[i + kNumTests]);
        }
        else if constexpr (std::is_invocable_v<
                               Func, TLWE<lvl0param>&, const TLWE<lvl0param>&,
                               const TLWE<lvl0param>&, const TLWE<lvl0param>&,
                               const GateKey&>) {
            func(cres[i], c[i], c[i + kNumTests], c[i + kNumTests * 2], gk);
            p[i] = chegk(p[i], p[i + kNumTests], p[i + kNumTests * 2]);
        }
        else if constexpr (std::is_invocable_v<
                               Func, TLWE<lvl0param>&, const TLWE<lvl0param>&,
                               const TLWE<lvl0param>&, const TLWE<lvl0param>&,
                               const TLWE<lvl0param>&, const GateKey&>) {
            func(cres[i], c[i], c[i + kNumTests], c[i + kNumTests * 2],
                 c[i + kNumTests * 3], gk);
            p[i] = chegk(p[i], p[i + kNumTests], p[i + kNumTests * 2],
                         p[i + kNumTests * 3]);
        }
        else {
            std::cout << "Invalid Function" << std::endl;
        }
    }
    end = chrono::system_clock::now();
    vector<uint8_t> p2(cres.size());
    p2 = bootsSymDecrypt(cres, sk);
    for (int i = 0; i < kNumTests; i++) assert(p[i] == p2[i]);
    std::cout << "Passed" << std::endl;
    double elapsed =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
            .count();
    cout << elapsed / kNumTests << "ms" << endl;
}

int main()
{
    const uint32_t kNumTests = 10;

    cout << "------ Key Generation ------" << endl;
    SecretKey* sk = new SecretKey();
    GateKey* gk = new GateKey(*sk);
    // MUX Need 3 input
    vector<uint8_t> p(3 * kNumTests);
    vector<TLWE<lvl0param>> cres(kNumTests);
    vector<TLWE<lvl0param>> c(3 * kNumTests);
    bool correct;

    correct = true;

    Test("NOT", HomNOT, NotChegk, p, cres, c, kNumTests, *sk, *gk);
    Test("COPY", HomCOPY, CopyChegk, p, cres, c, kNumTests, *sk, *gk);
    Test("NAND", HomNAND, NandChegk, p, cres, c, kNumTests, *sk, *gk);
    Test("OR", HomOR, OrChegk, p, cres, c, kNumTests, *sk, *gk);
    Test("ORYN", HomORYN, OrYNChegk, p, cres, c, kNumTests, *sk, *gk);
    Test("ORNY", HomORNY, OrNYChegk, p, cres, c, kNumTests, *sk, *gk);
    Test("AND", HomAND, AndChegk, p, cres, c, kNumTests, *sk, *gk);
    Test("ANDYN", HomANDYN, AndYNChegk, p, cres, c, kNumTests, *sk, *gk);
    Test("ANDNY", HomANDNY, AndNYChegk, p, cres, c, kNumTests, *sk, *gk);
    Test("XOR", HomXOR, XorChegk, p, cres, c, kNumTests, *sk, *gk);
    Test("XNOR", HomXNOR, XnorChegk, p, cres, c, kNumTests, *sk, *gk);
    Test("MUX", HomMUX, MuxChegk, p, cres, c, kNumTests, *sk, *gk);
    Test("NMUX", HomNMUX, NMuxChegk, p, cres, c, kNumTests, *sk, *gk);
    Test("ConstantZero", HomCONSTANTZERO, ConstantZeroChegk, p, cres, c,
         kNumTests, *sk, *gk);
    Test("ConstantOne", HomCONSTANTONE, ConstantOneChegk, p, cres, c, kNumTests,
         *sk, *gk);

    return 0;
}