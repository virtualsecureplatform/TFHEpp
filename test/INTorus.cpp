#include<cuhe++.hpp>
#include <iostream>
#include <random>
#include <cassert>
#include <boost/multiprecision/cpp_int.hpp>

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
        __int128_t a = dist(engine);
        __int128_t b = dist(engine);
        INTorus A(a);
        INTorus B(b);
        assert((A+B).value==((a+b)%P));
    }
    cout<<"Add PASS"<<endl;
    //Sub Test
    for(int i = 0; i<numTest;i++){
        __int128_t a = dist(engine);
        __int128_t b = dist(engine);
        INTorus A(a);
        INTorus B(b);
        assert((A-B).value==((a+(P-b))%P));
    }
    cout<<"Sub PASS"<<endl;
    //Mul Test
    for(int i = 0; i<numTest;i++){
        __int128_t a = dist(engine);
        __int128_t b = dist(engine);
        INTorus A(a);
        INTorus B(b);
        assert((A*B).value==((a*b)%P));
    }
    cout<<"Mul PASS"<<endl;
    //Lsh Test
    for(int upper = 1;upper<=6;upper++){
        for(int l = 32*(upper-1); l<32*upper; l++){
            cout<<l<<endl;
            for(int i = 0; i<numTest;i++){
                __int128_t temp = dist(engine);
                INTorus A(temp);
                mp::cpp_int a = temp;
                // cout<<static_cast<uint64_t>(a)<<":";
                if((A<<l).value!=(static_cast<uint64_t>((a<<l)%P))){
                    cout<<(A<<l).value<<":"<<static_cast<uint64_t>((a<<l)%P)<<endl;
                    cout<<(A<<l).value - static_cast<uint64_t>((a<<l)%P)<<endl;
                }
                assert((A<<l).value==(static_cast<uint64_t>((a<<l)%P)));
            }
            __int128_t temp = ((1UL<<32)-1)<<32;
            INTorus A(temp);
            mp::cpp_int a = temp;
            if((A<<l).value!=(static_cast<uint64_t>((a<<l)%P))){
                cout<<(A<<l).value<<":"<<static_cast<uint64_t>((a<<l)%P)<<endl;
             }
            assert((A<<l).value==(static_cast<uint64_t>((a<<l)%P)));
        }
        cout<<"Lsh"<<upper*32<<"PASS"<<endl;
    }
}