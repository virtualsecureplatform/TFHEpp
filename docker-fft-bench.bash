#!/bin/bash
docker build -t tfhepp -f Dockerfile-ubuntu2404 .
docker run -v $PWD:/TFHEpp tfhepp bash -c "cd /TFHEpp && bash ./fft-bench.bash"