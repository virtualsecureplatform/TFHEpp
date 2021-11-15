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

template <class TLWEParam, class Func, class Chegk>
void Test(string type, Func func, Chegk chegk, vector<uint8_t> pres,
          vector<uint8_t> p, vector<TLWE<TLWEParam>> cres,
          vector<TLWE<TLWEParam>> c, const int kNumTests, const SecretKey& sk,
          const EvalKey& ek)
{
    random_device seed_gen;
    default_random_engine engine(seed_gen());
    uniform_int_distribution<uint32_t> binary(0, 1);
    cout << "------ Test " << type << " Gate ------" << endl;
    cout << "Number of tests:\t" << kNumTests << endl;

    for (uint8_t& i : p) i = binary(engine);
    c = bootsSymEncrypt<TLWEParam>(p, sk);

    chrono::system_clock::time_point start, end;
    start = chrono::system_clock::now();
    for (int i = 0; i < kNumTests; i++) {
        if constexpr (std::is_invocable_v<Func, TLWE<TLWEParam>&>) {
            func(cres[i]);
            pres[i] = chegk();
        }
        else if constexpr (std::is_invocable_v<Func, TLWE<TLWEParam>&,
                                               const TLWE<TLWEParam>&>) {
            func(cres[i], c[i]);
            pres[i] = chegk(p[i]);
        }
        else if constexpr (std::is_invocable_v<
                               Func, TLWE<TLWEParam>&, const TLWE<TLWEParam>&,
                               const TLWE<TLWEParam>&, const EvalKey&>) {
            func(cres[i], c[i], c[i + kNumTests], ek);
            pres[i] = chegk(p[i], p[i + kNumTests]);
        }
        else if constexpr (std::is_invocable_v<
                               Func, TLWE<TLWEParam>&, const TLWE<TLWEParam>&,
                               const TLWE<TLWEParam>&, const TLWE<TLWEParam>&,
                               const EvalKey&>) {
            func(cres[i], c[i], c[i + kNumTests], c[i + kNumTests * 2], ek);
            pres[i] = chegk(p[i], p[i + kNumTests], p[i + kNumTests * 2]);
        }
        else if constexpr (std::is_invocable_v<
                               Func, TLWE<TLWEParam>&, const TLWE<TLWEParam>&,
                               const TLWE<TLWEParam>&, const TLWE<TLWEParam>&,
                               const TLWE<TLWEParam>&, const EvalKey&>) {
            func(cres[i], c[i], c[i + kNumTests], c[i + kNumTests * 2],
                 c[i + kNumTests * 3], ek);
            pres[i] = chegk(p[i], p[i + kNumTests], p[i + kNumTests * 2],
                            p[i + kNumTests * 3]);
        }
        else {
            std::cout << "Invalid Function" << std::endl;
        }
    }
    end = chrono::system_clock::now();
    vector<uint8_t> p2(cres.size());
    p2 = bootsSymDecrypt<TLWEParam>(cres, sk);
    for (int i = 0; i < kNumTests; i++) assert(pres[i] == p2[i]);
    std::cout << "Passed" << std::endl;
    double elapsed =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
            .count();
    cout << elapsed / kNumTests << "ms" << endl;
}

template <class T>
void RunTest()
{
    // using F_T1 = void (*)(TLWE<T>&);
    // using F_T2 = void (*)(TLWE<T>&, const TLWE<T>&);
    // using F_T3 =
    //     void (*)(TLWE<T>&, const TLWE<T>&, const TLWE<T>&, const EvalKey&
    //     ek);
    // using F_T4 = void (*)(TLWE<T>&, const TLWE<T>&, const TLWE<T>&,
    //                       const TLWE<T>&, const EvalKey& ek);

    const uint32_t kNumTests = 10;

    if constexpr (std::is_same_v<T, lvl0param>)
        cout << "------ Test of lvl0param ------" << endl;
    else if constexpr (std::is_same_v<T, lvlMparam>)
        cout << "------ Test of lvlMparam ------" << endl;

    cout << "------ Key Generation ------" << endl;
    SecretKey* sk = new SecretKey();
    TFHEpp::EvalKey ek;
    ek.emplacebkfft<TFHEpp::lvl01param>(*sk);
    ek.emplaceiksk<TFHEpp::lvl10param>(*sk);
    // MUX Need 3 input
    vector<uint8_t> pres(kNumTests);
    vector<uint8_t> p(3 * kNumTests);
    vector<TLWE<T>> cres(kNumTests);
    vector<TLWE<T>> c(3 * kNumTests);

    Test<T>("NOT", HomNOT, NotChegk, pres, p, cres, c, kNumTests, *sk, ek);
    Test<T>("COPY", HomCOPY, CopyChegk, pres, p, cres, c, kNumTests, *sk, ek);
    Test<T>("NAND", HomNAND, NandChegk, pres, p, cres, c, kNumTests, *sk, ek);
    Test<T>("OR", HomOR, OrChegk, pres, p, cres, c, kNumTests, *sk, ek);
    Test<T>("ORYN", HomORYN, OrYNChegk, pres, p, cres, c, kNumTests, *sk, ek);
    Test<T>("ORNY", HomORNY, OrNYChegk, pres, p, cres, c, kNumTests, *sk, ek);
    Test<T>("AND", HomAND, AndChegk, pres, p, cres, c, kNumTests, *sk, ek);
    Test<T>("ANDYN", HomANDYN, AndYNChegk, pres, p, cres, c, kNumTests, *sk,
            ek);
    Test<T>("ANDNY", HomANDNY, AndNYChegk, pres, p, cres, c, kNumTests, *sk,
            ek);
    Test<T>("XOR", HomXOR<T>, XorChegk, pres, p, cres, c, kNumTests, *sk, ek);
    Test<T>("XNOR", HomXNOR<T>, XnorChegk, pres, p, cres, c, kNumTests, *sk,
            ek);
    Test<T>("MUX", HomMUX<T>, MuxChegk, pres, p, cres, c, kNumTests, *sk, ek);
    Test<T>("NMUX", HomNMUX<T>, NMuxChegk, pres, p, cres, c, kNumTests, *sk,
            ek);
    Test<T>("ConstantZero", HomCONSTANTZERO, ConstantZeroChegk, pres, p, cres,
            c, kNumTests, *sk, ek);
    Test<T>("ConstantOne", HomCONSTANTONE, ConstantOneChegk, pres, p, cres, c,
            kNumTests, *sk, ek);
}

int main()
{
    RunTest<lvl1param>();
#ifdef ENABLE_AXELL
    RunTest<lvlMparam>();
#endif  // ENABLE_AXELL
    return 0;
}
