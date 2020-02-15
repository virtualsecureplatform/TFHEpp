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
sed -i '79c\    inline void add4(double * res, const double * a, const double * b) {' ../src/libtfhe/fft_processors/spqlios/spqlios-fft-impl.cpp && sed -i '83c\    inline void sub4(double * res, const double * a, const double * b) {' ../src/libtfhe/fft_processors/spqlios/spqlios-fft-impl.cpp
cmake ../src -DENABLE_TESTS=on -DENABLE_NAYUKI_PORTABLE=off -DENABLE_NAYUKI_AVX=off -DENABLE_SPQLIOS_AVX=off -DENABLE_SPQLIOS_FMA=on -DCMAKE_BUILD_TYPE=optim
make
./test/test-gate-bootstrapping-spqlios-fma
```

This is for TFHE-10ms.
```
git clone https://github.com/virtualsecureplatform/tfhe-10ms.git --recursive
cd tfhe-10ms
mkdir build
cd build
cmake ../src -DENABLE_TESTS=on -DENABLE_NAYUKI_PORTABLE=off -DENABLE_NAYUKI_AVX=off -DENABLE_SPQLIOS_AVX=off -DENABLE_SPQLIOS_FMA=on -DCMAKE_BUILD_TYPE=optim
make
./test/test-bootstrapping-fft-spqlios-fma 
```

If you have Docker on your system, this will do above on docker.

```
git clone https://github.com/virtualsecureplatform/TFHEpp
cd TFHEpp
docker build -t tfheppbench .
```