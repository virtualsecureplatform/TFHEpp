#include <tfhe++.hpp>
using namespace TFHEpp;

#include <iostream>
using namespace std;

uint8_t ConstantZeroCheck(){
  return 0;
}

uint8_t ConstantOneCheck(){
  return 1;
}

uint8_t NotCheck(const uint8_t in) {
  return (~in)&0x1;
}

uint8_t CopyCheck(const uint8_t in) {
  return in;
}

uint8_t NandCheck(const uint8_t in0, uint8_t in1) {
  return 1 - in0 * in1;
}

uint8_t OrCheck(const uint8_t in0, const uint8_t in1) {
  return (in0 + in1) > 0;
}

uint8_t OrYNCheck(const uint8_t in0, const uint8_t in1) {
  return (in0 + (1-in1)) > 0;
}

uint8_t OrNYCheck(const uint8_t in0, const uint8_t in1) {
  return ((1-in0) + in1) > 0;
}

uint8_t AndCheck(const uint8_t in0, const uint8_t in1) {
  return in0 * in1;
}

uint8_t AndYNCheck(const uint8_t in0, const uint8_t in1) {
  return in0 * (1-in1);
}

uint8_t AndNYCheck(const uint8_t in0, const uint8_t in1) {
  return (1-in0) * in1;
}

uint8_t XorCheck(const uint8_t in0, const uint8_t in1) {
  return (in0 + in1) & 0x1;
}

uint8_t XnorCheck(const uint8_t in0, const uint8_t in1) {
  return (~(in0 ^ in1)) & 0x1;
}

uint8_t MuxCheck(const uint8_t inc, const uint8_t in1, const uint8_t in0){
  return (inc>0)?in1:in0;
}

template <class Func, class Check>
void Test(
        string type, 
        Func func,
        Check check,
        vector<uint8_t> p, 
        vector<TLWElvl0> cres, 
        vector<TLWElvl0> c, 
        const int kNumTests, 
        SecretKey& sk,
        CloudKey& ck)
{
  random_device seed_gen;
  default_random_engine engine(seed_gen());
  uniform_int_distribution<uint32_t> binary(0, 1);
  cout<< "------ Test "<<type<<" Gate ------" <<endl;
  cout<< "Number of tests:\t" << kNumTests <<endl;
  bool correct = true;
  int cnt_failures = 0;

  for (uint8_t &i:p) i = binary(engine);
  c = bootsSymEncrypt(p,sk);


  for (int i = 0; i < kNumTests; i ++){
    if constexpr (std::is_invocable_v<Func, TLWElvl0&>){
        func(cres[i]);
        p[i] = check();
    }else if constexpr (std::is_invocable_v<Func, TLWElvl0&, const TLWElvl0&>){
        func(cres[i], c[i]);
        p[i] = check(p[i]);
    }else if constexpr(std::is_invocable_v<Func, TLWElvl0&, const TLWElvl0&, const TLWElvl0&, const CloudKey&>){
        func(cres[i], c[i], c[i + kNumTests], ck);
        p[i] = check(p[i], p[i + kNumTests]);
    }else if constexpr(std::is_invocable_v<Func, TLWElvl0&, const TLWElvl0&, const TLWElvl0&, const TLWElvl0&, const CloudKey&>){
        func(cres[i], c[i], c[i + kNumTests], c[i + kNumTests*2], ck);
        p[i] = check(p[i], p[i + kNumTests], p[i + kNumTests*2]);
    }else if constexpr(std::is_invocable_v<Func, TLWElvl0&, const TLWElvl0&, const TLWElvl0&, const TLWElvl0&, const TLWElvl0&, const CloudKey&>){
        func(cres[i], c[i], c[i + kNumTests], c[i + kNumTests*2], c[i + kNumTests*3], ck);
        p[i] = check(p[i], p[i + kNumTests], p[i + kNumTests*2], p[i + kNumTests*3]);
    }else{
        std::cout << "Invalid Function" << std::endl;
    }
  }
  vector<uint8_t> p2(cres.size());
  p2 = bootsSymDecrypt(cres, sk);
  for(int i = 0;i<kNumTests;i++) assert(p[i] == p2[i]);
  std::cout<<"Passed"<<std::endl;
}

int main() {
    const uint32_t kNumTests = 10;

    cout<< "------ Key Generation ------" <<endl;
    SecretKey* sk = new SecretKey();
    CloudKey* ck = new CloudKey(*sk);
    //MUX Need 3 input
    vector<uint8_t> p(3 * kNumTests);
    vector<TLWElvl0> cres(kNumTests);
    vector<TLWElvl0> c(3 * kNumTests);
    bool correct;

    correct = true;

    Test("NOT", HomNOT, NotCheck, p, cres, c, kNumTests, *sk, *ck);
    Test("COPY", HomCOPY, CopyCheck, p, cres, c, kNumTests, *sk, *ck);
    Test("NAND", HomNAND, NandCheck, p, cres, c, kNumTests, *sk, *ck);
    Test("OR", HomOR, OrCheck, p, cres, c, kNumTests, *sk, *ck);
    Test("ORYN", HomORYN, OrYNCheck, p, cres, c, kNumTests, *sk, *ck);
    Test("ORNY", HomORNY, OrNYCheck, p, cres, c, kNumTests, *sk, *ck);
    Test("AND", HomAND, AndCheck, p, cres, c, kNumTests, *sk, *ck);
    Test("ANDYN", HomANDYN, AndYNCheck, p, cres, c, kNumTests, *sk, *ck);
    Test("ANDNY", HomANDNY, AndNYCheck, p, cres, c, kNumTests, *sk, *ck);
    Test("XOR", HomXOR, XorCheck, p, cres, c, kNumTests, *sk, *ck);
    Test("XNOR", HomXNOR, XnorCheck, p, cres, c, kNumTests, *sk, *ck);
    Test("MUX", HomMUX, MuxCheck, p, cres, c, kNumTests, *sk, *ck);
    Test("ConstantZero", HomCONSTANTZERO, ConstantZeroCheck, p, cres, c, kNumTests, *sk, *ck);
    Test("ConstantOne", HomCONSTANTONE, ConstantOneCheck, p, cres, c, kNumTests, *sk, *ck);

    return 0;
}