# TFHEpp
TFHEpp is full Scracthed pure C++ Ver. of TFHE. TFHEpp is slightly(about 10%) faster than original TFHE and supports Circuit Bootstrapping.
TFHEpp depends on AVX2 becauise we use SPQLIOS FMA. If you want run TFHEpp without AVX2, see spqlios++ branch. It include pure C++ implementation of SPQLIOS as header only library, but slow.

# Speed Test

Following Code measure how many time homomorphic NAND takes on your computer with TFHEpp. 
```
git clone https://github.com/virtualsecureplatform/TFHEpp
cd TFHEpp
mkdir build
cd build
cmake ..
make
./test/nand 
```

If you want to run semantically equivalent test on original TFHE, run below code.
```
git clone https://github.com/tfhe/tfhe.git --recursive
cd tfhe
mkdir build
cd build
cmake ../src -DENABLE_TESTS=on -DENABLE_NAYUKI_PORTABLE=off -DENABLE_NAYUKI_AVX=off -DENABLE_SPQLIOS_AVX=off -DENABLE_SPQLIOS_FMA=on -DCMAKE_BUILD_TYPE=optim
make
./test/test-gate-bootstrapping-spqlios-fma
```