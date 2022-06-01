FROM ubuntu:22.04

LABEL maintainer="naoki9911(Naoki MATSUMOTO) <m.naoki9911@gmail.com>"

# install build dependencies
RUN apt-get update && apt-get upgrade -y
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y build-essential g++ libomp-dev cmake git libgoogle-perftools-dev

# build TFHE
RUN git clone --recursive --depth 1 https://github.com/tfhe/tfhe.git  && mkdir tfhe/build
WORKDIR /tfhe/build
RUN cmake ../src -DENABLE_TESTS=on -DENABLE_NAYUKI_PORTABLE=off -DENABLE_NAYUKI_AVX=off -DENABLE_SPQLIOS_AVX=off -DENABLE_SPQLIOS_FMA=on -DCMAKE_BUILD_TYPE=optim && make -j$(nproc)

# build TFHEpp
COPY . /TFHEpp
RUN mkdir /TFHEpp/build
WORKDIR /TFHEpp/build
RUN cmake .. -DENABLE_TEST=ON -DCMAKE_CXX_COMPILER=g++ -DCMAKE_BUILD_TYPE=Release && make -j$(nproc)

# To run benchmark test
#$ echo "TFHEpp" && /TFHEpp/build/test/nand && echo "original TFHE" && /tfhe/build/test/test-gate-bootstrapping-spqlios-fma 
