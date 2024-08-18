#!/bin/bash
# This script is used to benchmark the FFT performance of TFHEpp
echo "Benchmarking HomNAND using different FFTs" | tee /tmp/log.txt
cmake . -G Ninja -B /tmp/build -DENABLE_TEST=ON
cd /tmp/build
ninja
echo "Benchmarking SPQLIOS" | tee -a /tmp/log.txt
./test/nand | tee -a /tmp/log.txt

cmake . -G Ninja -B /tmp/build -DENABLE_TEST=ON -DUSE_AVX512=ON
cd /tmp/build
ninja
echo "Benchmarking SPQLIOS AVX512" | tee -a /tmp/log.txt
./test/nand | tee -a /tmp/log.txt

cd /TFHEpp
rm -rf /tmp/build
cmake . -G Ninja -B /tmp/build -DENABLE_TEST=ON  -DUSE_CONCRETE_FFT=ON
cd /tmp/build
ninja
ninja
echo "Benchmarking concrete-fft" | tee -a /tmp/log.txt
./test/nand | tee -a /tmp/log.txt

cd /TFHEpp
rm -rf /tmp/build
cmake . -G Ninja -B /tmp/build -DENABLE_TEST=ON  -DUSE_CONCRETE_FFT=ON -DUSE_AVX512=ON
cd /tmp/build
ninja
ninja
echo "Benchmarking concrete-fft AVX512" | tee -a /tmp/log.txt
./test/nand | tee -a /tmp/log.txt

cd /TFHEpp
rm -rf /tmp/build
source /opt/intel/oneapi/setvars.sh
cmake . -G Ninja -B /tmp/build -DENABLE_TEST=ON -DUSE_MKL=ON
cd /tmp/build
ninja
echo "Benchmarking MKL" | tee -a /tmp/log.txt
./test/nand | tee -a /tmp/log.txt

cd /TFHEpp
rm -rf /tmp/build
cmake . -G Ninja -B /tmp/build -DENABLE_TEST=ON -DUSE_FFTW3=ON
cd /tmp/build
ninja
echo "Benchmarking FFTW3" | tee -a /tmp/log.txt
./test/nand | tee -a /tmp/log.txt

cat /tmp/log.txt
