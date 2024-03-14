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

template <class P, class Func, class Chegk>
void Test(string type, Func func, Chegk chegk, vector<uint8_t> p,
          vector<TLWE<P>> cres, vector<TLWE<P>> c, const int kNumTests,
          const SecretKey& sk, const EvalKey& ek)
{
    random_device seed_gen;
    default_random_engine engine(seed_gen());
    uniform_int_distribution<uint32_t> binary(0, 1);
    cout << "------ Test " << type << " Gate ------" << endl;
    cout << "Number of tests:\t" << kNumTests << endl;

    for (uint8_t& i : p) i = binary(engine);
    c = bootsSymEncrypt<P>(p, sk);

    chrono::system_clock::time_point start, end;
    start = chrono::system_clock::now();
    for (int i = 0; i < kNumTests; i++) {
        if constexpr (std::is_invocable_v<Func, TLWE<P>&>) {
            func(cres[i]);
            p[i] = chegk();
        }
        else if constexpr (std::is_invocable_v<Func, TLWE<P>&,
                                               const TLWE<P>&>) {
            func(cres[i], c[i]);
            p[i] = chegk(p[i]);
        }
        else if constexpr (std::is_invocable_v<Func, TLWE<P>&, const TLWE<P>&,
                                               const TLWE<P>&,
                                               const EvalKey&>) {
            func(cres[i], c[i], c[i + kNumTests], ek);
            p[i] = chegk(p[i], p[i + kNumTests]);
        }
        else if constexpr (std::is_invocable_v<Func, TLWE<P>&, const TLWE<P>&,
                                               const TLWE<P>&, const TLWE<P>&,
                                               const EvalKey&>) {
            func(cres[i], c[i], c[i + kNumTests], c[i + kNumTests * 2], ek);
            p[i] = chegk(p[i], p[i + kNumTests], p[i + kNumTests * 2]);
        }
        else if constexpr (std::is_invocable_v<Func, TLWE<P>&, const TLWE<P>&,
                                               const TLWE<P>&, const TLWE<P>&,
                                               const TLWE<P>&,
                                               const EvalKey&>) {
            func(cres[i], c[i], c[i + kNumTests], c[i + kNumTests * 2],
                 c[i + kNumTests * 3], ek);
            p[i] = chegk(p[i], p[i + kNumTests], p[i + kNumTests * 2],
                         p[i + kNumTests * 3]);
        }
        else {
            std::cout << "Invalid Function" << std::endl;
        }
    }
    end = chrono::system_clock::now();
    vector<uint8_t> p2(cres.size());
    p2 = bootsSymDecrypt<P>(cres, sk);
    for (int i = 0; i < kNumTests; i++) assert(p[i] == p2[i]);
    std::cout << "Passed" << std::endl;
    double elapsed =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
            .count();
    cout << elapsed / kNumTests << "ms" << endl;
}

template <class P>
void RunTest()
{
    const uint32_t kNumTests = 10;

    if constexpr (std::is_same_v<P, lvl0param>)
        cout << "------ Test of lvl0param ------" << endl;
    else if constexpr (std::is_same_v<P, lvl1param>)
        cout << "------ Test of lvl1param ------" << endl;

    cout << "------ Key Generation ------" << endl;
    SecretKey* sk = new SecretKey();
    TFHEpp::EvalKey ek;
    ek.emplacebkfft<TFHEpp::lvl01param>(*sk);
    ek.emplaceiksk<TFHEpp::lvl10param>(*sk);
    // MUX Need 3 input
    vector<uint8_t> p(3 * kNumTests);
    vector<TLWE<P>> cres(kNumTests);
    vector<TLWE<P>> c(3 * kNumTests);

    if constexpr (std::is_same_v<P, lvl0param>) {
        cout << "lvl0param" << endl;
        Test<P>("NOT", TFHEpp::HomNOT<P>, NotChegk, p, cres, c, kNumTests, *sk,
                ek);
        Test<P>("COPY", TFHEpp::HomCOPY<P>, CopyChegk, p, cres, c, kNumTests,
                *sk, ek);
        Test<P>("NAND", TFHEpp::HomNAND<lvl01param, lvl1param::μ, lvl10param>,
                NandChegk, p, cres, c, kNumTests, *sk, ek);
        Test<P>("OR", TFHEpp::HomOR<lvl01param, lvl1param::μ, lvl10param>,
                OrChegk, p, cres, c, kNumTests, *sk, ek);
        Test<P>("ORYN", TFHEpp::HomORYN<lvl01param, lvl1param::μ, lvl10param>,
                OrYNChegk, p, cres, c, kNumTests, *sk, ek);
        Test<P>("ORNY", TFHEpp::HomORNY<lvl01param, lvl1param::μ, lvl10param>,
                OrNYChegk, p, cres, c, kNumTests, *sk, ek);
        Test<P>("AND", TFHEpp::HomAND<lvl01param, lvl1param::μ, lvl10param>,
                AndChegk, p, cres, c, kNumTests, *sk, ek);
        Test<P>("ANDYN", TFHEpp::HomANDYN<lvl01param, lvl1param::μ, lvl10param>,
                AndYNChegk, p, cres, c, kNumTests, *sk, ek);
        Test<P>("ANDNY", TFHEpp::HomANDNY<lvl01param, lvl1param::μ, lvl10param>,
                AndNYChegk, p, cres, c, kNumTests, *sk, ek);
        Test<P>("XOR", TFHEpp::HomXOR<lvl01param, lvl1param::μ, lvl10param>,
                XorChegk, p, cres, c, kNumTests, *sk, ek);
        Test<P>("XNOR", TFHEpp::HomXNOR<lvl01param, lvl1param::μ, lvl10param>,
                XnorChegk, p, cres, c, kNumTests, *sk, ek);
        Test<P>("MUX", TFHEpp::HomMUX<P>, MuxChegk, p, cres, c, kNumTests, *sk,
                ek);
        Test<P>("NMUX", TFHEpp::HomNMUX<P>, NMuxChegk, p, cres, c, kNumTests,
                *sk, ek);
        Test<P>("ConstantZero", TFHEpp::HomCONSTANTZERO<P>, ConstantZeroChegk,
                p, cres, c, kNumTests, *sk, ek);
        Test<P>("ConstantOne", TFHEpp::HomCONSTANTONE<P>, ConstantOneChegk, p,
                cres, c, kNumTests, *sk, ek);
    }
    else if constexpr (std::is_same_v<P, lvl1param>) {
        cout << "lvl1param" << endl;
        Test<P>("NOT", TFHEpp::HomNOT<P>, NotChegk, p, cres, c, kNumTests, *sk,
                ek);
        Test<P>("COPY", TFHEpp::HomCOPY<P>, CopyChegk, p, cres, c, kNumTests,
                *sk, ek);
        Test<P>("NAND", TFHEpp::HomNAND<lvl10param, lvl01param, lvl1param::μ>,
                NandChegk, p, cres, c, kNumTests, *sk, ek);
        Test<P>("OR", TFHEpp::HomOR<lvl10param, lvl01param, lvl1param::μ>,
                OrChegk, p, cres, c, kNumTests, *sk, ek);
        Test<P>("ORYN", TFHEpp::HomORYN<lvl10param, lvl01param, lvl1param::μ>,
                OrYNChegk, p, cres, c, kNumTests, *sk, ek);
        Test<P>("ORNY", TFHEpp::HomORNY<lvl10param, lvl01param, lvl1param::μ>,
                OrNYChegk, p, cres, c, kNumTests, *sk, ek);
        Test<P>("AND", TFHEpp::HomAND<lvl10param, lvl01param, lvl1param::μ>,
                AndChegk, p, cres, c, kNumTests, *sk, ek);
        Test<P>("ANDYN", TFHEpp::HomANDYN<lvl10param, lvl01param, lvl1param::μ>,
                AndYNChegk, p, cres, c, kNumTests, *sk, ek);
        Test<P>("ANDNY", TFHEpp::HomANDNY<lvl10param, lvl01param, lvl1param::μ>,
                AndNYChegk, p, cres, c, kNumTests, *sk, ek);
        Test<P>("XOR", TFHEpp::HomXOR<lvl10param, lvl01param, lvl1param::μ>,
                XorChegk, p, cres, c, kNumTests, *sk, ek);
        Test<P>("XNOR", TFHEpp::HomXNOR<lvl10param, lvl01param, lvl1param::μ>,
                XnorChegk, p, cres, c, kNumTests, *sk, ek);
        Test<P>("MUX", TFHEpp::HomMUX<P>, MuxChegk, p, cres, c, kNumTests, *sk,
                ek);
        Test<P>("NMUX", TFHEpp::HomNMUX<P>, NMuxChegk, p, cres, c, kNumTests,
                *sk, ek);
        Test<P>("ConstantZero", TFHEpp::HomCONSTANTZERO<P>, ConstantZeroChegk,
                p, cres, c, kNumTests, *sk, ek);
        Test<P>("ConstantOne", TFHEpp::HomCONSTANTONE<P>, ConstantOneChegk, p,
                cres, c, kNumTests, *sk, ek);
    }
}

int main()
{
    RunTest<lvl1param>();
    RunTest<lvl0param>();
    return 0;
}