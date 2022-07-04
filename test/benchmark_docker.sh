#!/bin/bash
set -x


SCRIPT_DIR=$(cd $(dirname $0); pwd)
cd $SCRIPT_DIR

./show_info.sh

echo "TFHEpp"
/TFHEpp/build/test/nand

echo "original TFHE"
/tfhe/build/test/test-gate-bootstrapping-spqlios-fma
