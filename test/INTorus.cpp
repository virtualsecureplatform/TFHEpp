#include<cuhe++.hpp>
#include <iostream>
#include <random>
#include <cassert>
#include <boost/multiprecision/cpp_int.hpp>
#include <gmp.h>
#include <gmpxx.h>

using namespace std;
using namespace cuHEpp;
namespace mp = boost::multiprecision;

int main(){
    constexpr int numTest = 10000;
    random_device seed_gen;
    default_random_engine engine(seed_gen());
    uniform_int_distribution<uint64_t> dist(0, P);

    // Add Test
    for(int i = 0; i<numTest;i++){
        __uint128_t a = dist(engine);
        __uint128_t b = dist(engine);
        INTorus A(a);
        INTorus B(b);
        assert((A+B).value==((a+b)%P));
    }
    cout<<"Add PASS"<<endl;
    //Sub Test
    for(int i = 0; i<numTest;i++){
        __uint128_t a = dist(engine);
        __uint128_t b = dist(engine);
        INTorus A(a);
        INTorus B(b);
        assert((A-B).value==((a+(P-b))%P));
    }
    cout<<"Sub PASS"<<endl;
    //Mul Test
    for(int i = 0; i<numTest;i++){
        __uint128_t a = dist(engine);
        __uint128_t b = dist(engine);
        INTorus A(a);
        INTorus B(b);
        cout<<(A*B).value<<":"<<static_cast<uint64_t>((a*b)%P)<<endl;
        assert((A*B).value==((a*b)%P));
    }
    cout<<"Mul PASS"<<endl;
    //Lsh Test
    for(int upper = 1;upper<=6;upper++){
        for(int l = 32*(upper-1); l<32*upper; l++){
            cout<<l<<endl;
            for(int i = 0; i<numTest;i++){
                uint64_t temp = dist(engine);
                INTorus A(temp);
                // mp::cpp_int a = temp;
                mpz_class a = temp;
                // cout<<static_cast<uint64_t>(a)<<":";
                // if((A<<l).value!=(static_cast<uint64_t>((a<<l)%P))){
                if((A<<l).value!=(a<<l)%P){
                    mpz_class res = (a<<l)%P;
                    cout<<(A<<l).value<<":"<<res.get_str()<<endl;
                    // cout<<(A<<l).value - static_cast<uint64_t>((a<<l)%P)<<endl;
                }
                // assert((A<<l).value==(static_cast<uint64_t>((a<<l)%P)));
                assert((A<<l).value==(a<<l)%P);
            }
            __uint128_t temp = ((1UL<<32)-1)<<32;
            INTorus A(temp);
            mp::cpp_int a = temp;
            if((A<<l).value!=(static_cast<uint64_t>((a<<l)%P))){
                cout<<(A<<l).value<<":"<<static_cast<uint64_t>((a<<l)%P)<<endl;
                cout<<"Here"<<endl;
             }
            assert((A<<l).value==(static_cast<uint64_t>((a<<l)%P)));
        }
        cout<<"Lsh"<<upper*32<<"PASS"<<endl;
    }
}